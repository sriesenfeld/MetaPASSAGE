#!/usr/bin/perl
# S.J. Riesenfeld
# Updated: Oct 2010

uexse strict;
use warnings;
use File::Spec;

# This a script example to show the creation of reference databases
# using a gene family for which peptide and DNA sequences are
# available. 
# 
# Usage: perl example1_alt_refdb_create.pl
#
# This script assumes that a file 'rpoB.pep' of peptide sequences and
# a file 'rpoB.fna' of corresponding DNA sequences are present in the
# current directory. It also assumes that an HMM model for the rpoB
# gene family is present in a file called 'rpoB.hmm' (built so that
# partial sequences can be aligned to it).

my $file_basename = 'rpoB';  # The name of the files (without
			     # extensions) containing sequences for a
			     # protein family of 
my $refdb_partialname = 'RefDB_rpoB';  # Basename (including path, but
				       # not extensions) for output
				       # reference database files.
my @num_seqs = (20, 200);  # Total number of sequences to choose.
my @num_seqs_maxPD = (5, 50);  # Choose 25 sequences to maximize
			       # phylogenetic diversity and then
			       # populate the rest randomly.

# Create two reference databases: one of size 20 sequences and one of
# size 200 sequences; each one has a quarter of their sequences being
# chosen by maximizing phylogenetic diversity (the rest are chosen
# uniformly at random).
foreach my $i (0,1) {

    my $num = $num_seqs[$i];
    my $num_maxPD = $num_seqs_maxPD[$i];
    my $refdb_basename = $refdb_partialname.'_'.$num;

    my $cmd = "perl ../createRefDB.pl -f $file_basename ".
	" -n $num -m $num_maxPD -b $refdb_basename -t Genome.tree";

    system($cmd);
    print "\n";
}
