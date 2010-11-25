#!/usr/bin/perl
# S.J. Riesenfeld
# Updated: Oct 2010

use strict;
use warnings;
use File::Spec;

# This a script example to show how to create a reference database
# with non-AMPHORA sequences using the MetaPASSAGE package.
# 
# No input required; the script assumes there is a fasta file called
# '16SrRNA.fna' in the current directory. It contains 16s rRNA
# sequences that have been converted to DNA for use with the
# simulation pipeline.  This file is included with the distribution to
# help the user test the software.

my $seqs_file_basename = '16SrRNA';

my $num_seqs = 200;  # Total number of sequences to choose.
		     # AMPHORA's rpoB.pep file as reference sequences

my $refdb_basename = 'RefDB_16SrRNA'.'_'.$num_seqs;  
  # Basename (including path, but not extensions) for output reference
  # database files

my $cmd = "perl ../createRefDB.pl -n $num_seqs -b $refdb_basename".
    " -f $seqs_file_basename";

system($cmd);
