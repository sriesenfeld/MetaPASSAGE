#!/usr/bin/perl 
# S.J. Riesenfeld
# Jan 2011
#
# Used to sample sequences (with replacement) from a fasta file (-s
# option), according to frequencies given by a different file (-f
# option), until there are a specified total number of unique
# sequences that have been sampled. Two files are output: one is the
# file of sequences sampled (each sequence is only printed once
# there); the other is a file with extension ".mprf" for use with
# MetaSim that gives the number of times each sequence in the first
# output file was sampled.
#
# This can be used to simulate a subsample from a larger population.
#
# The functionality of this script is expected to be added to
# MetaPASSAGE.

use strict;
# use warnings;

use Getopt::Long;
use File::Spec;

my ($seqfile, $freqfile, $out_seqfile, $out_mprffile);

my $num = 200;
GetOptions ( "s=s" => \$seqfile,
	     "f=s" => \$freqfile,
	     "o=s" => \$out_seqfile,
	     "m=s" => \$out_mprffile,
	     "n:i" => \$num,
    );

open(FREQ, "$freqfile") or die "Cannot open $freqfile: $!";

open (OUTS, ">$out_seqfile") or die "Cannot open $out_seqfile: $!";
open(OUTM, ">$out_mprffile") or die "Cannot open $out_mprffile: $!";

my $id_pattern = qr/([^_]*)_(\S*)/;
my $freq_file_pattern = qr/$id_pattern\s*(\d*)/;

my @pop_ids;
my @unique_ids;
while (my $line=<FREQ>) {
    if ($line =~ $freq_file_pattern) {
	my $id = $1.' '.$2;
	if (grep($id eq $_, @unique_ids)) {
	    print "ID $id is not unique in file.\n";
	} else {
	    push(@unique_ids, $id);
	}
	my $n_copies = $3;
	foreach my $i (1.. $n_copies) {
	    # 1 elt for every copy of sequence in pop
	    push(@pop_ids, $id);
	}
    } else {
	print "Unexpected formatting this line: $line";
    }
}
print scalar(@pop_ids). " sequences in population.\n";
print scalar(@unique_ids). " unique ids in population.\n";

my ($community_ar, $com_unique_ids_ar)
    = sample_with_replacement_until_unique_num($num,\@pop_ids);
print scalar(@$community_ar). " sequences in community sampled.\n";
print scalar(@$com_unique_ids_ar). " unique ids in community sampled.\n";

# print mprf file
print OUTM "# MetaSim profile for community corresponding to sequences in $out_seqfile.\n";
foreach my $id (@$com_unique_ids_ar) {
    my $freq = scalar(grep($_ eq $id, @$community_ar));
    print OUTM "$freq name \"$id\"\n";
}

# print seq file
open(SEQ, "$seqfile") or die "Cannot open $seqfile: $!";
my @printed_seqs=();
my $printing = 0;
while(my $line = <SEQ>) {
    if ($line =~ /^>(.*)/) {
	my $header = $1;
	if ($header =~ $id_pattern) {
	    $header = $1.' '.$2;
	} else {
	    print "Could not reformat header: $header\n";
	}
	# id is present in community so make sure sequence is printed, but only once
	if (scalar(grep( $_ eq $header, @$com_unique_ids_ar))
	    and (! scalar(grep ($_ eq $header, @printed_seqs)))) {
	    push(@printed_seqs, $header);
	    print OUTS ">$header\n";
	    $printing=1;
	} else {
	    $printing=0;
	}
    } elsif ($printing) {
	print OUTS $line;
    }
}

close SEQ;
close FREQ;
close OUTM;
close OUTS;


sub sample_with_replacement_until_unique_num ($$) {
    my ($num_unique, $array_ref) = @_;
    if (($num_unique > scalar(@$array_ref)) or ($num_unique <= 0)) {
	die "First argument should be a positive integer smaller than the size of the array.\n";
    }
    my @sample;
    my @unique_elts;
    while(1) {
	my $random_index = int(rand(scalar(@$array_ref)));	
	my $elt = $array_ref->[$random_index];
	# print "Sampled index $random_index, corresponding to element $elt.\n";
	push(@sample, $elt);
	if (! scalar(grep($_ eq $elt, @unique_elts))) {
	    # print "New element sampled.\n";
	    push(@unique_elts, $elt);
	    if (scalar(@unique_elts) >= $num_unique) {
		last;
	    }
	}
    }
    return (\@sample, \@unique_elts);;
}

