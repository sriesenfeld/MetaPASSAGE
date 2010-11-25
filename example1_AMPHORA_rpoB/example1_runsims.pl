#!/usr/bin/perl
# S.J. Riesenfeld
# Updated: Nov. 2010
#
# Assumes that the script example1_refdb_create.pl has been run.
#
# Usage: 
#   perl example1_runsims.pl <Options>
# <Options>:
#   --amphora_path <AMPHORA_PATH>
#       where <AMPHORA_PATH> is the path to (and including) the
#       AMPHORA top-level directory; default is '/usr/local/AMPHORA/';
#
#   --hmmer2_path <HMMER2_PATH>
#       where <HMMER2_PATH> is the path to the directory containing
#       HMMER 2 binaries, including hmmpfam; necessary for alignment
#       using AMPHORA; default is '/usr/local/hmmer2/bin/';
#
#   --metasim_path <METASIM_PATH>
#       where <METASIM_PATH> is the path to the metasim installation
#       directory; default is '/usr/local/metasim/'; optional to use
#       this; if path exists, just enables a log of MetaSim's run
#       information to be copied and stored for each simulation.
#
# Note: We are expecting this script to be run from the
# example1_AMPHORA_rpoB directory, which is expected to be a
# subdirectory of the MetaPASSAGE directory. The MetaSim error
# configuration file errormodel-80bp.mconf is expected to be in the
# example1_AMPHORA_rpoB_directory. Change the paths below if these
# assumptions are not valid.
#
#
# Function: Example of using the simulation pipeline MetaPASSAGE with
# AMPHORA's repository of rpoB protein sequences. It demonstrates the
# tight integration of the simulation pipeline with AMPHORA. This
# analysis uses tools built into AMPHORA to find and align sequences
# to the hmm model it stores. The alignment is trimmed by AMPHORA with
# a hand-curated mask. 
#
# See the example1_alt_rpoB directory for a similar analysis that does
# not use AMPHORA.
#
# This example suggests how the effect of limited prior knowledge
# about the universe of sequences for a gene family may be studied by
# using different reference databases. Moreover, it shows how a novel
# sequencing error model can be incorporated easily in the simulation
# of reads by MetaSim. It incorporates a filtering of the reads after
# they are generated in order to retain only one per distinct
# full-length sequence. Finally, it gives an idea of how the pipeline
# may be run in batch.

use strict;
use warnings;
use File::Spec;
use Getopt::Long;

my ( $amphora_path, $gene, $dirbasename, $filebasename, $hmmer2_path,
     $refdb_partialname, $refdb_basename, $n_src, $n_target, $n_reads,
     $file_basename_from_refdb, 
     $error_conf_file, $seqs_file_basename, $filter_setting, 
     $metasim_path, $cmd );

$amphora_path = '/usr/local/AMPHORA/';  
$metasim_path = '/usr/local/metasim/'; 
$hmmer2_path = '/usr/local/hmmer2/bin/';

GetOptions(
    '--amphora_path:s' => \$amphora_path,
    '--hmmer2_path:s' => \$hmmer2_path,
    '--metasim_path:s' => \$metasim_path,
    );
if (! (-d $amphora_path)) {
    die "Cannot find AMPHORA here: $amphora_path!". 
	" Please give the correct path using option '--amphora_path'.\n";
}
if (!(-d $hmmer2_path) or !(-e File::Spec->catfile($hmmer2_path, 'hmmpfam'))) {
    die "Cannot find hmmpfam of HMMER2 here: $hmmer2_path!".
	" Please give the correct path using option '--hmmer2_path'.\n";
}
$gene = 'rpoB';  # The gene symbol AMPHORA uses for rpoB sequences
$error_conf_file = 'errormodel-80bp.mconf'; 
  # We make use of an Empirical error model with MetaSim rather than
  # the pipeline default of Sanger with no error. This model
  # configuration is for 80bp Illumina reads and was downloaded from
  # the MetaSim website.

my @refdb_sizes = (20, 200);  
$refdb_partialname = 'RefDB_rpoB_'; # The refdb is a subset of the
				    # universe of AMPHORA rpoB
				    # sequences which simulates
				    # limited prior knowledge of
				    # full-length rpoB sequences;
				    # created with the script
				    # example1_refdb_create.pl.
$n_src = 100;  # Number of different rpoB sequences that are present in
	       # the population; these full-length sequences are used
	       # to generate the reads.
$n_reads = 500;  # 500 reads gives an average of 5 reads per unique
		 # source sequence, if the distribution over the
		 # selected full-length sequences is uniform.
$filter_setting = -1;  # Filters the final read set so there is at
		       # most one read per original source
		       # sequence. This is useful if you want to
		       # compare two analyses where one involves
		       # full-length sequences and the other reads.

$dirbasename= 'sims_rpoB_';
$filebasename='sim';
$file_basename_from_refdb = 'RefDB_';

my $num_trials = 2;  # Would typically make the number of trials for
		     # each setting of the parameters bigger for
		     # statistical purposes.
my $drop_threshold = 10; # Only drop reads shorter than 10bp

foreach my $size (@refdb_sizes) {
    $refdb_basename = $refdb_partialname.$size;

    my $big_dir = $dirbasename.$file_basename_from_refdb.$size;
    if (-d $big_dir) {
	die "Directory $big_dir already exists!\n";
    }
    mkdir($big_dir);

    foreach my $i (1..$num_trials) { # can be made parallel easily
	
	print "\nPerforming simulation $i of $num_trials with $size-sequence Reference DB.\n";
	my $out_basename = $filebasename.$i.'_'.$gene.'_'.$file_basename_from_refdb.$size;
	my $dir = File::Spec->catdir($big_dir, $out_basename);

	$cmd = "perl ../MetaPASSAGE.pl --out_dir $dir ". 
	    " --out_basename $out_basename ".
	    " --amphora_path $amphora_path ". 
	    " --gene $gene ".
	    " --refdb_basename $refdb_basename ".
	    " --num_seqs $n_src --sim".
	    " --ms_error_conf_file $error_conf_file ".
	    " --metasim_path $metasim_path".
	    " --num_reads $n_reads --drop_len $drop_threshold ".
	    " --blast --target_num_reads $filter_setting ".
	    " --align --align_source ".
	    " --hmmer2_path $hmmer2_path";
	# AMPHORA stores sequences for the gene family and has its own
	# models for alignment so there is no need to provide a model
	# or peptide sequences. It is important that an additional
	# file of DNA sequences corresponding to the peptide sequences
	# for the gene family is created and stored in AMPHORA. See
	# the README and the notes in example1_refdb_create.pl.
	
	print "Command: $cmd\n";
	eval {
	    system($cmd);
	};
	if ($@) {
	    die "Error testing MetaPASSAGE.pl: $@!\n";
	}
	print "\n";
    }
}
print "\nAll simulations complete.\n";
