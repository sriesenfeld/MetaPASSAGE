#!/usr/bin/perl
# S.J. Riesenfeld
# Updated: Nov. 2010
#
# Assumes that the script example1_alt_refdb_create.pl has been run.
#
# Usage:
#   perl example1_alt_runsims.pl <Options>
# <Options>:
#
#   --metasim_path <METASIM_PATH>
#       where <METASIM_PATH> is the path to the metasim installation
#       directory; default is '/usr/local/metasim/'; optional to use
#       this; if path exists, just enables a log of MetaSim's run
#       information to be copied and stored for each simulation.
#
# Note: We are expecting this script to be run from the
# example1_alt_rpoB directory, which is expected to be a subdirectory
# of the MetaPASSAGE directory. The MetaSim error configuration file
# errormodel-80bp.mconf is expected to be in the
# example1_alt_rpoB_directory. Change the paths below if these
# assumptions are not valid.
#
# Function: Example of using the simulation pipeline MetaPASSAGE with
# the rpoB family of protein sequences. This example does a similar
# analysis to the one in the example1_AMPHORA_rpoB but without going
# through AMPHORA. The sequence files and Genome tree file are
# borrowed from AMPHORA, simply as examples. Any protein family would
# work. The HMM for rpoB was built using HMMER3 from the AMPHORA
# rpoB.aln file, converted to Stockholm format.
#
# When MetaPASSAGE is used outside AMPHORA, it calls whichever HMMER
# functions are in the user's path. Since the model included was built
# with HMMER3 hmmbuild, the user should have HMMER 3 installed so the
# alignment is properly done (or the user can substitute a different
# rpoB model built using a different version of HMMER.)
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

my ( $in_seqs_file_basename, $dirbasename, $filebasename,
     $model, $refdb_partialname, $refdb_basename,
     $file_basename_from_gene, $file_basename_from_refdb, 
     $n_src, $n_target, $n_reads, $error_conf_file,
     $filter_setting, $metasim_path, $cmd );

$metasim_path = '/usr/local/metasim/'; 

GetOptions(
    '--metasim_path:s' => \$metasim_path,
    );

$in_seqs_file_basename = 'rpoB'; # rpoB.pep and rpoB.fna containing
				 # peptide and DNA sequences,
				 # respectively, in this directory
$model = 'rpoB.hmm'; # Created from AMPHORA's rpoB.aln using HMMER3
		     # hmmbuild

$error_conf_file = 'errormodel-80bp.mconf'; 
  # We make use of an Empirical error model with MetaSim rather than
  # the pipeline default of Sanger with no error. This model
  # configuration is for 80bp Illumina reads and was downloaded from
  # the MetaSim website.

my @refdb_sizes = (20, 200);  
$refdb_partialname = 'RefDB_rpoB_'; # The refdb is a subset of the
				    # universe of rpoB sequences which
				    # simulates limited prior
				    # knowledge of full-length rpoB
				    # sequences; created with the
				    # script
				    # example1_alt_refdb_create.pl.
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
$file_basename_from_gene='rpoB_';
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
	my $out_basename = $filebasename.$i.'_'.$file_basename_from_gene.$file_basename_from_refdb.$size;
	my $dir = File::Spec->catdir($big_dir, $out_basename);
	
	# Simulates reads from the rpoB sequences, blasts them to
	# figure out what frame to translate them in, and aligns them
	# and the reference database to the model.
	$cmd = "perl ../MetaPASSAGE.pl --out_dir $dir ". 
	    " --out_basename $out_basename ".
	    " --seqs_file_basename $in_seqs_file_basename ".
	    " --refdb_basename $refdb_basename ".
	    " --num_seqs $n_src --sim".
	    " --ms_error_conf_file $error_conf_file ".
	    " --metasim_path $metasim_path".
	    " --num_reads $n_reads --drop_len $drop_threshold ".
	    " --blast --target_num_reads $filter_setting ".
	    " --align --align_source ".
	    " --model $model";
	
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
