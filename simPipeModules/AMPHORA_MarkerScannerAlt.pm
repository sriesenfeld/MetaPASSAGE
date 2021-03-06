#!/usr/bin/perl
# S.J. Riesenfeld, Martin Wu
# MetaPASSAGE module: AMPHORA_MarkerScannerAlt
# Last updated: Nov 2010
#
# This module was created by Samantha J. Riesenfeld from the
# MarkerScanner.pl script in AMPHORA (more information below).
# Changes are noted in comments below.  The module does *not* exactly
# mimic the functionality of the MarkerScanner.pl script in
# AMPHORA. (For example, it skips the blast against E.Coli.)  It was
# edited for use with MetaPASSAGE.pl.
#
# NOTE: This version of AMPHORA uses HMMER 2. The path for HMMER 2 is
# set in simPipeVars.pm.  AMPHORA also requires BioPerl, so make sure
# the BioPerl packages are correctly loaded below.

# NOTE from AMPHORA:
#
# AMPHORA (version 1.0) An Automated Phylogenomic Inference Pipeline
# for bacterial sequences.  Copyright 2008 by Martin Wu
#
# AMPHORA is free software: you may redistribute it and/or modify its
# under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# any later version.
#
# AMPHORA is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details
# (http://www.gnu.org/licenses/).  For any other inquiries send an
# Email to Martin Wu martinwu@ucdavis.edu

# Riesenfeld: making this into a package so it can be loaded with
# other modules in simPipeModules

package AMPHORA_MarkerScannerAlt;

use strict;
# Riesenfeld: Added these use lines, for use with MetaPASSAGE.pl
use warnings;
use lib 'simPipeModules';
use simPipeVars qw(:amphora);

require Exporter;
our @ISA = qw(Exporter);
our @EXPORT_OK = qw(markerScannerAlt);

# Riesenfeld: Added in case BioPerl is not installed
my $bioperl_missing=0;
eval {
    require Bio::SeqIO;
    Bio::SecIO->import();
};
if ($@) {
    $bioperl_missing=1;
}
eval {
    require Bio::SearchIO;   
    Bio::SearchIO->import();
};
if ($@) {
    $bioperl_missing=1;
    warn("Note: AMPHORA functionality is not available without installation of BioPerl.\n");
}

my (%markerlist, %seq, %Hmmlength) = ();
my $AMPHORA_home = $amphora_path;
my $in_file;

sub markerScannerAlt($;$$) {

    # Riesenfeld: Check for BioPerl:
    if ($bioperl_missing) {
	die "AMPHORA_MarkerScannerAlt: Cannot use AMPHORA functionality without installing BioPerl!\n";
    }

    # Riesenfeld: Slightly changed input message and added input paths
    $in_file = shift(@_);    
    my $usage = qq~ Usage: markerScannerAlt <fasta-file> <AMPHORA-path> <hmmer2-path>~;
    
    die $usage unless (defined ($in_file) and (-f $in_file));
    my $path = shift(@_);
    if (defined ($path)) {
	$AMPHORA_home = $path;
    }
    if (!(-d $AMPHORA_home)) {
	die "Cannot find AMPHORA at $AMPHORA_home!\n";
    }
    $path = shift(@_);
    if (defined($path)) {
	$hmmer2_path = $path;
    }
    if (! (-d $hmmer2_path)) {
	die "Cannot find HMMER2 binaries at $hmmer2_path!\n";
    }
    
    get_marker_list();
    read_marker_hmms();

    # Riesenfeld: Skipping this blast step; we don't want to only identify E. coli types.
    #   (Previously converted blast step from using WU-BLAST to using NCBI BLAST.)
    # system ("cp $ARGV[0] $$.query");
    # Blastp search
    # my $cmd = $ncbi_blast_path."formatdb -i $$.query >&/dev/null";
    # system ($cmd);
    # $cmd = $ncbi_blast_path."blastall -p blastp -d $$.query -i $AMPHORA_home/Marker/markers.fas -e 0.1 -v 50000 -b 50000 -z 5000 > $$.blastp";
    # system($cmd);
    # get_blast_hits();
    
    # Riesenfeld: Need to initialize %seq since get_blast_hits() is not called
    #   Lines below added (as well as sub definition)
    system ("cp $in_file $$.candidate");
    get_seq_ids();
    
    # HMM search
    # Riesenfeld: using File::Spec->catfile instead of string concatenation
    my $hmmpfam_path = File::Spec->catfile($hmmer2_path, 'hmmpfam');
    if (! (-e $hmmpfam_path)) {
	die "Cannot find hmmpfam!\n";
    }
    my $cmd = "$hmmpfam_path -Z 5000 -E 1e-3 ".
	File::Spec->catfile($AMPHORA_home, "Marker/markers.swhmm")." $$.candidate > $$.hmmsearch";	
    system ($cmd);
    
    # fix the number of sequences in the database for E-value calculation
    get_hmm_hits();
    
    # clean up
    system ("rm $$.*");
}

####################################################################################################################
sub get_marker_list {
    # Riesenfeld: using File::Spec->catfile instead of string concatenation
    my $marker_list = File::Spec->catfile($AMPHORA_home, "Marker/marker.list");
    open (IN, $marker_list)
	|| die "Can't open $marker_list: $!\n";
    while (<IN>) {
	chop;
	/^(\S+)/;
	$markerlist{$1} = 1;
    }
    close IN;
}

sub read_marker_hmms {
    # Riesenfeld: using File::Spec->catfile instead of string concatenation
    my $markers = File::Spec->catfile($AMPHORA_home, "Marker/markers.swhmm");
    open (IN, $markers) || die "Can't open $markers: $!\n";
    my $name;
    while (<IN>) {
        chop;
	if (/^NAME\s+(\S+)/) {
	    $name = $1;
	}
	elsif (/^LENG\s+(\d+)/) {
	    $Hmmlength{$name} = $1;
	}
    }
    close IN;
}

sub get_blast_hits {
    my %hits = ();
    
    my $in = new Bio::SearchIO('-file' => "$$.blastp");
    while (my $result = $in->next_result) {	
	while( my $hit = $result->next_hit ) {
	    $hits{$hit->name()} = 1;
	}
    }
    unless (%hits) {
	system("rm $$.*");
	exit(1) 
    }
    my $seqin = new Bio::SeqIO('-file'=>$in_file);
    my $seqout = new Bio::SeqIO('-file'=>">$$.candidate",'-format'=>'fasta');
    
    while (my $seq = $seqin->next_seq) {
	if ($hits{$seq->id}) {
	    $seq{$seq->id} = $seq;
	    $seqout->write_seq($seq);
	}
    }
}

# Riesenfeld: Added sub below, need to initialize %seq if get_blast_hits does not get called
sub get_seq_ids {
    my $seqin = new Bio::SeqIO('-file'=>$in_file);
    while (my $seq = $seqin->next_seq) {
	$seq{$seq->id} = $seq;
    }
}

sub get_hmm_hits {
    my $hmmsearch = new Bio::SearchIO ('-file'=>"$$.hmmsearch", '-format' => 'hmmer');
    my $marker;
    
    my %hits =();
    
  RESULT:while ( my $result = $hmmsearch->next_result ){
      my $query = $result->query_name();
      my ($query_match, $hit_match, $best_evalue, $hitname) = ();
      if (my $hit = $result->next_hit() ){
	  if ( !(defined $best_evalue) or $best_evalue > $hit->significance) {		
	      $best_evalue = $hit->significance;
	      $hitname = $hit->name;
	      $query_match = 0;
	      $hit_match = 0;
	      while (my $hsp = $hit->next_hsp) {
		  $query_match += ($hsp->end('query')-$hsp->start('query') );
		  $hit_match += ($hsp->end('hit')-$hsp->start('hit')); 
	      }
	      $query_match  /= $seq{$query}->length;
	      $hit_match /= $Hmmlength{$hit->name()};
	  }	
      }
      ## Riesenfeld: Added the check to see if $hitname is defined
      next RESULT unless (defined($hitname) and $markerlist{$hitname});
      next RESULT if ( ($query_match < 0.7) and ($hit_match < 0.7) );		#ignore the hit if the match is partial 
      ## Riesenfeld: this hit is ignored if the match is partial relative to BOTH the query and the hit      
      $hits{$hitname}{$query} = 1;	
  }
    unless (%hits) {
	system("rm $$.*");
	exit(1) 
    }
    for my $marker (keys %hits) {
	my $filename="$marker.pep";	    
	my $seqout = new Bio::SeqIO('-file'=>">$filename",'-format'=>'fasta');
	for my $seqid (keys %{$hits{$marker}}) {
	    $seqout->write_seq($seq{$seqid});
	}
    }
}

1;
