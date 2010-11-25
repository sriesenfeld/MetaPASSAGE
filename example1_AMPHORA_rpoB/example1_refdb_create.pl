#!/usr/bin/perl
# S.J. Riesenfeld
# Updated: Oct 2010

use strict;
use warnings;
use File::Spec;

# This a script example to show the creation of reference
# databases. It highlights the integration of MetaPASSAGE.pl with
# AMPHORA.
# 
# Argument: The path to the AMPHORA directory, or you can set this by
# changing the first line of code below, which sets the default path
# of /usr/local/AMPHORA.
# 
# Usage: perl example1_refdb_create.pl <AMPHORA_path>
#
# This script assumes that a file of DNA sequences corresponding to
# the AMPHORA gene family of interest already exists. To permit user
# testing, a file called 'rpoB.fna' of DNA sequences corresponding to
# the peptide sequences in 'AMPHORA/Reference/Sequences/rpoB.pep' is
# included. The DNA sequence file should be moved so it is in the same
# AMPHORA directory path as 'rpoB.pep'.
# 
# Note: To retrieve sequences for another gene family <gene> in
# AMPHORA, you can use the command:
#
# perl createRefDB.pl -g <gene> -r -d <AMPHORA_path>
#
# where <path_to_AMPHORA> is the full path to the AMPHORA
# directory. The retrieval takes several hours sometimes, depending on
# NCBI and GenBank. If it is interrupted, you can restart where you
# left off by running the command:
#
# perl createRefDB.pl -g <gene> -r -a -d <AMPHORA_path>
#
# The script outputs a summary of the retrieval process as well as any
# errors. For more details, see the comments in the source code of
# createRefDB.pl and the relevant modules in simPipeModules.
#
# Now we assume that 'rpoB.fna' exists and is stored in the same
# location as 'rpoB.pep' in the AMPHORA directory structure.

my $amphora_path = '/usr/local/AMPHORA/';  
# Customize this via an optional input argument so it contains the
# correct path to the AMPHORA directory.

my $alt_path = shift(@ARGV);
if (defined($alt_path)) {
    $amphora_path = $alt_path;
}

if (!(-d $amphora_path)) {
    die "Cannot find AMPHORA here: $amphora_path.". 
	" Please give the correct path as an argument.\n";
}
my $gene_family = 'rpoB';  # Protein family symbol from AMPHORA of
			   # interest.
my $refdb_partialname = 'RefDB_rpoB';  # Basename (including path, but
				       # not extensions) for output
				       # reference database files.
my @num_seqs = (20, 200);  # Total number of sequences to choose.
			   # AMPHORA's rpoB.pep file as reference
			   # sequences
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

    my $cmd = "perl ../createRefDB.pl -g $gene_family -d $amphora_path".
	" -n $num -m $num_maxPD -b $refdb_basename";

    system($cmd);
    print "\n";
}
