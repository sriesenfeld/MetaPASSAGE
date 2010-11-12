#!/usr/bin/perl
# S.J. Riesenfeld
# Updated: Oct 2010
# MetaPASSAGE module: alnQC
#
# Module that may eventually contain code for doing quality control on
# the alignment that MetaPASSAGE.pl outputs.  For now, it is very simple
# and consists only of a function for getting rid of short sequences.

package alnQC;

use strict;
use warnings;

use IO::File;

use simPipeVars qw(:alnqc :misc_file_ext);

require Exporter;
our @ISA = qw(Exporter);
our @EXPORT_OK = qw(drop_short_seqs);

=head2 drop_short_seqs

  Usage : drop_short_seqs(@args)    
  Function : Takes an input fasta file of sequences and writes a new
    fasta file of sequences such that the new file only contains the
    sequences in the input file that are longer than a given length
    threshold. The default length threshold is 25 characters (e.g.,
    base pairs if the sequence is DNA). If the threshold is <=0, then
    the input file is returned (with 0 as the number of sequences
    dropped).
  Example : 
    # Assume 'seqs.fna' is a fasta file of DNA sequences
    $file_out = drop_short_sequences('seqs.fna', 100);
    # $file_out contains only those sequences in 'seqs.fna' that are
    # longer than 100 bp.
  Returns : an array containing the new filename and the number of
    sequences dropped.
  Args : An array whose first element must be the input filename (a
    fasta file containing sequences). The second argument (optional)
    is a length threshold, which has a default of 25 that is set in
    simPipeVars.pm. The output file will contain sequences in the
    input file that have at least this length (in characters). The
    remaining optional argument is the desired name of the output
    file.

=cut

sub drop_short_seqs ($;$$) {
    
    my ($file_in, $threshold, $file_out) = @_;
    if (!defined($file_in)) {
	die "drop_short_seqs:  Require a file of sequences!\n";
    }
    if (!defined($threshold)) {
	$threshold = $seq_len_threshold_default;
    }
    if ($threshold <= 0) {
	return ($file_in, 0);
    }
    if (!defined($file_out)) {
	if ($file_in =~ /(.+)(\.[^\.]+)/) {
	    $file_out = $1.$qc_len_filename_part.$2;
	} else {
	    $file_out = $file_in.$qc_len_filename_part;
	}
    }
    my $dropped_ctr=0;

    # avoid dependence on BioPerl for now
    open (IN, $file_in) or die "Cannot open $file_in for read: $!\n";
    open (OUT, ">$file_out") or die "Cannot open $file_out for write: $!\n";

    my ($hdr, $seq);
    my $first_hdr = 1;
    my $len = 0;
    while (my $line = <IN>) {
	if ($line =~ /^>/) {
	    if (! $first_hdr) {
		# print "Length of current seq $hdr is: $len\n";
		if ($len >= $threshold) {
		    print OUT $hdr;
		    print OUT $seq;
		} else {
		    $dropped_ctr++;
		}
	    } else {
		$first_hdr=0;
	    }
	    $seq='';
	    $len=0;
	    $hdr = $line;
	} else {   
	    $seq = $seq.$line;
	    my $seq_part = $line;
	    chomp($seq_part);
	    $len = $len + length($seq_part);
	}
    }
    # if (defined($hdr)) {
	# print "Length of current seq $hdr is: $len\n";    
    # }
    if ($len >= $threshold) {
	print OUT $hdr;
	print OUT $seq;
    } elsif (! $first_hdr) {
	$dropped_ctr++;
    }
    # print "$dropped_ctr sequences dropped.\n";
    close(OUT);
    close(IN);
    return ($file_out, $dropped_ctr);
		     }

1;
