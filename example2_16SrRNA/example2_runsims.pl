#!/usr/bin/perl
# S.J. Riesenfeld
# Updated: Nov. 2010
#
# Assumes that the script example2_refdb_create.pl has been run.
# 
# Usage: perl example2_runsims.pl <Options>
# <Options>:
#   --model <16S_rRNA_model> 
#       where <16S_rRNA_model> is an INFERNAL CM model for aligning
#       16S rRNA sequences; such a model can be obtained by running
#       INFERNAL's cmbuild on an alignment of 16S rRNA sequences. (The
#       model we used to develop this example was built from a
#       hand-curated alignment of 16S rRNA sequences that is the
#       intellectual property of the Ribosomal Database Project (RDP)
#       and is not publicly available; we suggest asking RDP for it.)
#
#   --metasim_path <METASIM_PATH>
#       where <METASIM_PATH> is the path to the metasim installation
#       directory; default is '/usr/local/metasim/'; optional to use
#       this; if path exists, just enables a log of MetaSim's run
#       information to be copied and stored for each simulation.
#
# Note: We are expecting this script to be run from the directory
# example2_16SrRNA, which is expected to be a subdirectory of the
# MetaPASSAGE directory. Change the paths below if these assumptions
# are not valid.
#
# This is an example of using the simulation pipeline on 16S rRNA
# sequences (that have been converted first to DNA bases). It
# demonstrates how populations with skewed taxonomic distributions may
# be simulated. It also gives an illustration of how the reference
# database can be used in composing the simulated population, in this
# case to make sure that no exact copy of a sequence from the
# reference database appears in the population simulated. Finally, it
# shows how the pipeline can be run in batch.

use strict;
use warnings;
use File::Spec;
use Getopt::Long;

my ( $dirbasename, $filebasename, $refdb_basename, $n_src, $n_src_ref,
     $file_basename_from_gene, $file_basename_ratio, $n_target,
     $n_reads, $seqs_file_basename, $model, $metasim_path, $cmd );

$model = '16SrRNA_model.mod';  # The 16S model used to develop this
			       # example was built using INFERNAL from
			       # an RDP alignment.  You can ask the
			       # RDP for it.
$metasim_path = '/usr/local/metasim';
    # Optional to use this, just enables a log if MetaSim's run
    # information to be copied and stored for each simulation.						  

GetOptions(
    '--model:s' => \$model,
    '--metasim_path:s' => \$metasim_path,
);

if (! (-e $model)) {
    print "Cannot find model $model! Will skip alignment step.\n";
}

$seqs_file_basename = '16SrRNA';  # The universe (for this simulation)
				  # of 16S rRNA sequences that have
				  # been converted to DNA.
$refdb_basename = 'RefDB_16SrRNA_200'; # Subset of the universe which
				       # simulates limited prior
				       # knowledge of full-length 16S
				       # sequences, created with the
				       # script
				       # example2_refdb_create.pl.
$n_src = 200;  # Number of different 16S sequences that are present in
	       # the population; these full-length sequences are used
	       # to generate the reads.
$n_src_ref = 0;  # 0 intersection with the reference database, i.e.,
		 # only populate the community with sequences from the
		 # 16S set that are NOT in the reference database.
$n_reads = 1000;  # 1000 reads gives an average of 5 reads per unique
		  # source sequence, if the distribution over the
		  # selected full-length sequences is uniform.        

# If you don't have a model, you can also skip the alignment stage by
# leaving off the --align and --model options below.
$dirbasename= 'sims_16S_';
$filebasename='sim';
$file_basename_from_gene = '16S_';
$file_basename_ratio = 'ratio_';

my $num_trials = 2;  # Would typically make the number of trials for
		     # each setting of the parameters bigger for
		     # statistical purposes.

foreach my $skew_factor (1.5, 1) { # The 1.5 value skews the
				   # distribution with which MetaSim
				   # samples sequences heavily in
				   # favor of a few rpoB sequences,
				   # with a long tail; the value 1
				   # signifies a uniform distribution
    my $skew_str = $skew_factor;
    $skew_str =~ s/\./p/;  # MetaSim does not seem to like periods in the file name
    my $big_dir = $dirbasename.$file_basename_ratio.$skew_str;
    if (-d $big_dir) {
	die "Directory $big_dir already exists!\n";
    }
    mkdir($big_dir);
    
    foreach my $i (1..$num_trials) { # can be made parallel easily
	
	print "\nPerforming simulation $i of $num_trials with ratio set to $skew_factor.\n";
	my $out_basename = $filebasename.$i.'_'.$file_basename_from_gene.$file_basename_ratio.$skew_str;
	my $dir = File::Spec->catdir($big_dir, $out_basename);

	# Using the default read length here, which is about a mean of
	# 400bp. Using the pipeline default for sequencing error,
	# which is a Sanger model no error. Pads sequences to create a
	# realistic read length distribution. Then by default drops
	# reads that are shorter than 25bp, (useful for quality
	# control in later stages).
	$cmd = "perl ../MetaPASSAGE.pl --out_dir $dir ".
	    " --out_basename $out_basename ".
	    " --refdb_basename $refdb_basename --num_seqs $n_src ".
	    " --num_ref_seqs $n_src_ref --sim ".
	    " --tax_profile_ratio $skew_factor --num_reads $n_reads ".
	    " --metasim_path $metasim_path".
	    " --blast --type RNA2DNA ".
	    " --seqs_file_basename $seqs_file_basename";
	if (-e $model) {
	    $cmd = $cmd ." --align --model $model";	
	}
	# Reads are correctly oriented using blast, and then the
	# reference database and reads are all aligned to the model.
	
	print "Command: $cmd\n";
	eval {
	    system($cmd);
	};
	if ($@) {
	    die "Error testing simPipe.pl: $@!\n";
	}
	print "\n";
    }
}
print "\nAll simulations complete.\n";
