#!/usr/bin/perl
# S.J. Riesenfeld
# MetaPASSAGE module: retrieveSeqs
#
# Last updated: Nov 2010 (comments only); Code last updated in
# Sept. 2009
#
# Module for retrieving gene sequences (DNA) from the NCBI database
# for the MetaPASSAGE pipeline, designed particularly for use with
# AMPHORA protein families. May need to be adapted to work with
# non-AMPHORA sequences.  The retrieve_seqs sub is designed to retry
# calls to RefSeq and GenBank that generate errors because sometimes
# an error saying that an ID does not exist in the database is in fact
# due to a time-out or some other unrelated problem.  There are
# regular pauses after a bunch of retrieval calls in order not to get
# shut out of either database.

package retrieveSeqs;

use strict;
# use warnings;
use warnings FATAL => qw( all );  # to keep retrieval from going nuts
use Getopt::Std;

require Exporter;
our @ISA = qw(Exporter);
our @EXPORT_OK = qw(retrieve_seqs extract_gene_ids );
my $bioperl_missing=0;                                                                                            
eval {
    require Bio::DB::GenBank;
    Bio::DB::GenBank->import();
};
if ($@) {
    $bioperl_missing=1;
}
eval {
    require Bio::DB::RefSeq;
    Bio::DB::RefSeq->import();
};
if ($@) {
    $bioperl_missing=1;
}
eval {    
    require Bio::Factory::FTLocationFactory;
    Bio::Factory::FTLocationFactory->import();
};
if ($@) {
    $bioperl_missing=1;
}

=head2 retrieve_seqs

 Usage: retrieve_seqs(@args);

 Function: Writes the sequences for the genes as retrieved from NCBI
   to the output file.
 
 Example: 
   my $num_retrieved = 
      retrieve_seqs ('seq_hdrs.pep', 'seqs.fna', 1);

 Returns: the number of DNA sequences successfully retrieved.

 Args: An array of 3 arguments; last 1 is optional.
   arg 1: Name of a fasta file containing the sequence ids in the
          headers in the form in which they appear in AMPHORA
          Reference Sequence files, e.g., so that the gene id follows
          the '>' and is immediately followed by a '-'.
   arg 2: File name for output.
   arg 3: (optional) a flag for appending; if it is defined and
     non-zero, then <arg2> is checked to see what sequences it
     contains; any sequences with ids that match sequences in the
     input file are skipped and not retrieved; retrieve sequences are
     appended to the output file (rather than re-writing the output
     file).

=cut

sub retrieve_seqs($$;$) {
    if ($bioperl_missing) {
	die "Cannot do sequence retrieval without BioPerl installed!\n";
    }
    my ($infile, $outfile, $append) = @_;
    my $max_calls = 30;
    my $max_fails = 45;
    my $max_repeats = 3;
    my $sleep_time=30;
    my $quick_pause = 2;
    my @gene_ids = extract_gene_ids($infile);
    my @retrieved_ids;
    my @prev_gene_ids;
    if ($append) {
	@prev_gene_ids = extract_gene_ids($outfile);
    }
    print ''.scalar(@gene_ids). " gene ids for retrieval.\n "; 
       # .join(' ', @gene_ids)."\n";
    if ($append) {
	print ''.scalar(@prev_gene_ids). " gene ids already retrieved and written in file $outfile.\n";
    } 
    # unless ($type and ($type eq 'aa')) { $type = 'dna'; }
    if ($append) {
	open(OUT, ">>$outfile") or die "Cannot open $outfile for appending retrieved sequences: $!\n";
    } else {
	open(OUT,">$outfile") or die "Cannot open $outfile for writing retrieved sequences: $!\n";
    }
    ###### Creating bioperl modules's objects for sequence retrieval
    my $db = Bio::DB::RefSeq->new();
    # $db->verbose(2);
    my $gb = new Bio::DB::GenBank;
    # $gb->verbose(2);
    my $loc_factory = Bio::Factory::FTLocationFactory->new;
    
    my $num_refseq_fails=0;
    my $num_genbank_fails=0;
    my @failed_geneids;
    my @failed_genbank_ids;
    ###### Fetching sequences
    my $num_refseq_calls=0;
    my $num_genbank_calls=0;
    for my $geneId (@gene_ids) {
	if ($append and grep($_ eq $geneId, @prev_gene_ids)) {
	    # print "Gene id $geneId already retrieved.\n";
	    next;
	}
	print "Trying to retrieve sequence from RefSeq by gene id: $geneId... ";
	my $ptn_obj;
	my $err_flag=0;
	my $trial;
	for ($trial=0; $trial < $max_repeats; $trial++) {
	    if ($num_refseq_fails > $max_fails) {
		die "Too many ($num_refseq_fails) failed retrieval calls to RefSeq!  Quitting!\n";
	    }
	    if ($err_flag) {
		sleep ($quick_pause);
		print "Trying again:  re-trial $trial.\n";
		$err_flag=0;
	    }
	    $num_refseq_calls++;
	    eval {
		$ptn_obj = $db->get_Seq_by_id($geneId);
	    };
	    if ($@) {
		$num_refseq_fails++;
		print "A Bioperl exception occurred from Bio::DB::RefSeq::get_Seq_by_id: $@\n";
		$err_flag = 1;
	    } else {
		last;
	    }
	}
	if ($err_flag) {
	    push(@failed_geneids, $geneId);
	    print "Skipping geneId: $geneId.\n\n";
	    $err_flag=0;
	    next;
	}
	print "done.\n";

	my $nt_seq; my $gi = ''; my $orient = 1;my $range= '';my ($nt_acc,$loc_str);
	my $feat;
	my @feats = $ptn_obj->top_SeqFeatures;
        # foreach $feat ($ptn_obj->top_SeqFeatures) {
	while (@feats) {
	    $feat = shift(@feats);
	    if($feat->primary_tag eq 'CDS') {
		last;
	    }
	}
	if($feat->primary_tag ne 'CDS') {
	    next;
	}
	# get DNA sequence
	my @coded_by = $feat->each_tag_value('coded_by');
	($nt_acc,$loc_str) = split(/\:/,$coded_by[0]);my $comp_check = 0;
	if($nt_acc =~ s/complement\(//) {
	    $comp_check = 1;
	    $orient = -1;
	    $loc_str =~ s/\)//;
	}
	$range = $loc_str;
	my $nt_obj;
	print "Trying to retrieve sequence from GenBank by id: $nt_acc... ";
	for ($trial=0; $trial < $max_repeats; $trial++) {
	    if ($num_genbank_fails > $max_fails) {
		die "Too many ($num_genbank_fails) failed retrieval calls to GenBank!  Quitting!\n";
	    }
	    if ($err_flag) {
		sleep($quick_pause);
		print "Trying again: trial $trial.\n";
		$err_flag=0;
	    }
	    $num_genbank_calls++;
	    eval {
		$nt_obj = $gb->get_Seq_by_id($nt_acc);
	    };
	    if ($@) {
		$num_genbank_fails++;
		print "A Bioperl exception occurred from Bio::DB::GenBank::get_Seq_by_id: $@\n";
		$err_flag = 1;
	    } else {
		last;
	    }
	}
	if ($err_flag) {
	    push(@failed_geneids, $geneId);
	    push(@failed_genbank_ids, $nt_acc);
	    print "Skipping acc id: $nt_acc.\n\n";
	    next;
	}
	print "done.\n";	
	my $loc_obj = $loc_factory->from_string($loc_str);
	my $feat_obj = Bio::SeqFeature::Generic->new(-location=>$loc_obj);
	$nt_obj->add_SeqFeature($feat_obj);
	my $cds_obj = $feat_obj->spliced_seq;
	$nt_seq = $cds_obj->seq;
	if($comp_check == 1) {
	    $nt_seq =~ tr/ACGTacgt/TGCAtgca/;
	    $nt_seq = reverse($nt_seq);
	}
	my $description = $ptn_obj->desc;
	$description =~ s/\]./\]/;
	my $accession = $ptn_obj->accession_number();
	my $version = $ptn_obj->version();
	# my $fasta_header = '>gi|'.$gi.'|reft|'.$accession.'.'.$version.'| '.$description.'|'.$nt_acc.'|'.$range.'|'.$orient;	    
	# my $fasta_header = '>reft|'.$accession.'.'.$version.'| '.$description.'|'.$nt_acc.'|'.$range.'|'.$orient;
	my $fasta_header = '>'.$accession.'-'.$nt_acc.' '.$description.' (Ver. '.$version.')';
	my $formatted_seq = '';
	for(my $i =0; $i < length($nt_seq);$i+=60){
	    $formatted_seq = $formatted_seq . substr($nt_seq,$i,60) . "\n";
	}
	print OUT "$fasta_header\n$formatted_seq";	
	push(@retrieved_ids, $geneId);
	if ( (($num_refseq_calls % $max_calls) == 0) or (($num_genbank_calls % $max_calls) == 0) ) {
	    print "\n". scalar(@retrieved_ids)." sequences retrieved, \n".
		"$num_refseq_calls retrieval calls made to RefSeq,\n".
		"$num_genbank_calls retrieval calls made to GenBank.\n".
		"Pausing now for $sleep_time seconds...\n\n";
	    sleep ($sleep_time);
	}
    }
    print "Total of ".scalar(@retrieved_ids)." sequences retrieved.\n";
    print "Total of ".scalar(@failed_geneids)." sequences NOT retrieved due to BioPerl Exceptions:\n".
	'   '.join(', ', @failed_geneids)."\n";
    print ''.scalar(@failed_genbank_ids)." of these are due to failed retrievals from GenBank of the DNA sequence ids: \n".
	'   '.join(', ', @failed_genbank_ids)."\n";
    
    remove_sublist(\@failed_geneids, \@gene_ids);
    print "Total of ".scalar(@gene_ids) ." sequences in output file.\n";
    close OUT;
    my $bad_flag = check_output($outfile, \@gene_ids);
    if ($bad_flag) {
	die "Something went wrong retrieving or writing the gene sequences!\n";
    }
    return scalar(@retrieved_ids);
}


# first arg: reference to a fasta file containing the sequence ids in the headers
#            in the form in which they appear in AMPHORA Reference Sequence files, 
#            e.g., so that the gene id follows the '>' and is immediately followed by a '-'.
# Returns an array of the gene id parts of the sequence ids
sub extract_gene_ids ($) {
    my $infile = shift(@_);
    open (INF, "$infile") or die "Cannot open $infile: $!\n";
    my @gene_ids;
    my $header = qr{^>(\w+)-};
    while (my $line = <INF>) {
	if ($line =~ /$header/) {
	    my $gene_id = $1;
	    push (@gene_ids, $gene_id);
	}
    }
    close(INF);
    return @gene_ids;
}		

sub check_output($$) {
    my ($file, $expected_ids_ar) = @_;
    print "Checking output file.\n";
    my @seen_ids;
    my @unseen_ids = @{$expected_ids_ar};
    open (CHECK, "$file") or die "Cannot open $file: $!\n";
    my $bad_flag=0;
    my $header_flag = 0;
    my $header = qr{^>(\w+)-};
    # print ''.scalar(@unseen_ids). " expected ids: ".join(' ', @unseen_ids)."\n";
    while (my $line = <CHECK>) {
	if ($line =~ /$header/) {
	    my $gene_id = $1;
	    # print "Id: $gene_id.\n";
	    if ($header_flag) {
		print "Bad formatting in file $file: a header line appears two lines in a row!  (Missing sequence data?)\n";
		$bad_flag =1;
	    } 
	    $header_flag = 1;
	    if (grep($_ eq $gene_id, @seen_ids)) {
		print "Id $gene_id appears twice in file $file!\n";
		$bad_flag=1;
	    } 
	    if (! grep ($_ eq $gene_id, @unseen_ids) ) {
		print "Id $gene_id not expected in file $file!\n";
		$bad_flag=1;
	    } 
	    push(@seen_ids, $gene_id);
	    @unseen_ids = grep ($_ ne $gene_id, @unseen_ids);
	} else {
	    $header_flag = 0;
	}
    }
    if (scalar(@seen_ids) != scalar(@{$expected_ids_ar}) ) {
	print "Unexpected number of ids read: ".scalar(@seen_ids)." read, but ".scalar(@{$expected_ids_ar}).
	    " expected.\n";
	$bad_flag =1;
    }
    return $bad_flag;
}

sub in_list {
    my ($elt, $list_ref) =@_;
    if (grep($_ eq $elt, @$list_ref)) {
        return 1; }
    else { return 0;}
} 
        

sub remove_sublist {
    my ($r_list_ref, $list_ref) = @_;
    @$list_ref = grep( (!in_list($_, $r_list_ref)), @$list_ref);
}




1;
