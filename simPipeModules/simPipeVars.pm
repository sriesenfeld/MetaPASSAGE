#!/usr/bin/perl
# S.J. Riesenfeld
# Last updated: Nov 2010
# MetaPASSAGE module: simPipeVars
#
# Variables and simple functions used by MetaPASSAGE.pl and modules in
# SimPipeModules.
#
# All variables below can be customized except the following: 
#     $metasim_profile_ext, $metasim_fasta_ext, and $metasim_model_name;
#     $amphora_ref_file_ext, $amphora_ref_aln_file_ext;
# since these depend on the way MetaSim works internally.

package simPipeVars;
use strict;
use File::Spec;
require Exporter;
our @ISA = qw(Exporter);
our %EXPORT_TAGS = (simpipe => [qw($default_seqs_filename $ref_db_default_filename_part 
                                   $ref_db_default_size  				                          
                                   $sim_basename_default $src_filename_default
				   $single_reads_name $frames_file_ext $strand_file_ext 
                                   $final_src_seqs_name $orientation_name $align_log_filename 
                                   $pad_filename_part $nopad_filename_part &count_seqs
                                   AMPHORA CMALIGN HMMALIGN @aligners)],
		    seqtypes => [qw( PROTEIN RNA2DNA @types)],
		    amphora => [qw($amphora_path $amphora_ref_file_ext $amphora_ref_aln_file_ext
                                   $amphora_ref_seq_path $amphora_ref_seq_wpads_path 
                                   $amphora_marker_profile_path $amphora_scripts_path
                                   $amphora_genome_tree $amphora_markerscanner_altscript 
                                   $hmmer2_path $ncbi_blast_path $gene_symbol_default
                                   &define_amphora_paths)],
		    metasim => [qw($metasim_path $metasim_profile_ext 
				   $metasim_error_log_name 
				   $metasim_fasta_ext $metasim_filename_default
				   $metasim_sanger_name $metasim_error_conf_name
                                   $metasim_log_filename $metasim_file_name $metasim_hdr_pattern
				   $num_reads_default $mean_read_len_default $stddev_read_len_default
                                   $metasim_db_name_default
		                   $mean_clone_len_default $stddev_clone_len_default
                                   $seq_pad_size_default $sanger_mate_prob_default
                                   $sanger_err_start_default $sanger_err_end_default 
                                   $sanger_insert_default $sanger_delete_default)],
		    misc_file_ext => [qw(@fasta_alt_exts 
					 $stockholm_aln_ext $dna_fasta_ext $gen_fasta_ext $pep_fasta_ext $aln_fasta_ext
					 $blastx_file_ext $blastn_file_ext
					 $log_file_ext 
					 &get_basename &check_ext)],
		    mysql => [qw($mysql_db $mysql_uname $mysql_pword)],
                    maxpd => [qw($maxpd_refdb_default_num $maxpd_default_num $maxpd_margin_num
                                 $maxpd_filename %alt_names @taxa_to_avoid)],
                    blast => [qw(@blast_db_pep_exts @blast_db_na_exts $blast_max_expect 
                                 $blast_max_ratio $num_blastdb_hits_default)],
		    alnqc => [qw($seq_len_threshold_default 
                                 $qc_len_filename_part)],
    );

Exporter::export_ok_tags('amphora', 'seqtypes', 'metasim', 'misc_file_ext', 'mysql', 'simpipe', 
			 'maxpd','blast', 'alnqc');

# set some constants for sequence data types
use constant {
    PROTEIN  => 'Protein',       # input: peptide sequences 
    # RNA   => 'RNA',              # input: RNA sequences
    RNA2DNA  => 'RNA2DNA',       # input: RNA sequences that have been converted to DNA already
#   DNA   => 'DNA',            # DNA sequences given (not currently implemented)
};
our @types = (PROTEIN, RNA2DNA); #, RNA);

use constant {
    AMPHORA => 'AMPHORA',   # use AMPHORA to align sequences (peptides, only AMPHORA gene families)
    CMALIGN => 'cmalign',   # use cmalign from INFERNAL to align sequences (DNA or RNA)
    HMMALIGN => 'hmmalign', # use hmmalign from HMMER to align sequences (peptides)
};
our @aligners = (AMPHORA, CMALIGN, HMMALIGN);

# For MetaPASSAGE.pl specifically
our ($default_seqs_filename, $ref_db_default_filename_part,
     $ref_db_default_size, $sim_basename_default,
     $src_filename_default, $single_reads_name, $frames_file_ext,
     $strand_file_ext, $final_src_seqs_name, $orientation_name,
     $align_log_filename, $pad_filename_part, $nopad_filename_part);

$ref_db_default_filename_part = 'RefDB'; # Part of default file name
					 # for the database of
					 # reference sequences
$ref_db_default_size = 100; # Default size of reference database
$sim_basename_default = 'sim'; # Part of default base filename for
			       # files related to one simulation run
$src_filename_default = 'src'; # Part of label for file containing
			       # sample of full-length sequences used
			       # to generate reads
$single_reads_name = '-sr'; # Part of label for file containing at
			    # most one read per gene sequence
$frames_file_ext = '.frames'; # File extension for file of frame data
			      # written by MetaPASSAGE.pl
$strand_file_ext = '.orient'; # File extension for file of strand data
			      # written by MetaPASSAGE.pl
$final_src_seqs_name = '-final'; # Part of label for file containing
				 # the full-length sequences for the
				 # filtered set of reads
$orientation_name='-or'; # Part of label for file containing the
			 # correctly oriented reads (for RNA or DNA)
$align_log_filename = '-align.log'; # Part of label for file that
				    # contains stderr output from
				    # alignment program
$pad_filename_part = '-pd'; # Part of label for file containing padded
			    # sequences
$nopad_filename_part = '-rp'; # Part of label for file to indicate
			      # pads have been removed

# AMPHORA-related
our ($amphora_path, $amphora_ref_seq_path, $amphora_ref_seq_wpads_path,
     $amphora_ref_file_ext, $amphora_ref_aln_file_ext,
     $amphora_marker_profile_path, $amphora_genome_tree, $amphora_scripts_path, 
     $amphora_markerscanner_altscript, $hmmer2_path, $ncbi_blast_path, $gene_symbol_default);
# Path to AMPHORA $amphora_path variable can also be set by an input
# parameter to MetaPASSAGE.pl
$amphora_path = '/usr/local/AMPHORA/';

$amphora_ref_file_ext = '.pep'; # File extension used by AMPHORA for
				# file containing unaligned amino acid
				# reference sequences
$amphora_ref_aln_file_ext = '.aln'; # File extension used by AMPHORA
				    # for file containing aligned
				    # amino acid reference sequences
define_amphora_paths(); # Defines needed subdirectories and files for
			# using AMPHORA (see below)
$hmmer2_path = '/usr/local/hmmer2/bin/'; # Path to HMMER Ver 2
					 # executables (incl. hmmpfam)

# path for NCBI-BLAST programs (leave empty if they are in the default $PATH)
$ncbi_blast_path = '';
# default gene family to use if none is specified
$gene_symbol_default = 'rpoB';

# MetaSim-related
our ($metasim_profile_ext, $metasim_fasta_ext, $metasim_path,
     $metasim_sanger_name, $metasim_error_conf_name, $metasim_error_log_name, 
     $metasim_log_filename, $metasim_filename_default, $metasim_hdr_pattern,
     $num_reads_default, $mean_read_len_default, $stddev_read_len_default, 
     $mean_clone_len_default, $stddev_clone_len_default,
     $metasim_db_name_default, $seq_pad_size_default, 
     $sanger_mate_prob_default, $sanger_err_start_default, $sanger_err_end_default, 
     $sanger_insert_default, $sanger_delete_default);

$metasim_path = '/usr/local/metasim/'; # path to directory where MetaSim writes
				       # its log file by default (installation directory);
                                       # optional, but helps record information about MetaSim runs.
                                       # is only used if the path exists;
                                       # can be set here or by an input parameter to MetaPASSAGE.pl; 
$metasim_error_log_name = 'error.log';         # MetaSim output log name -- don't change

$metasim_profile_ext = '.mprf';  # Taxonomic profile file extension used by MetaSim -- don't change
$metasim_fasta_ext = '.fna';     # Fasta file extension used by MetaSim -- don't change
$metasim_sanger_name = '-Sanger';  # Error-free-model name used by MetaSim to label files -- don't change
$metasim_error_conf_name = '-Empirical';  # Empirical error-model name used by MetaSim to label files -- don't change
$metasim_log_filename = 'metasim.log';   # part of MetaPASSAGE local name of log for MetaSim
$metasim_filename_default = '-reads';    # Part of default name for metasim reads file used by MetaPASSAGE.pl
$metasim_hdr_pattern='^>(\S+).*SOURCE_.=.(\w+)';  # the pattern of the headers for the reads generated by MetaSim
$num_reads_default = 500;  # Default number of reads passed as an argument to MetaSim
$sanger_mate_prob_default=0;  # Default for Metasim Sanger model probability of paired reads
$sanger_err_start_default = 0; # Default for Metasim Sanger model initial error rate
$sanger_err_end_default = 0; # Default for Metasim Sanger model final error rate
$sanger_delete_default = 0; # Default for Metasim Sanger model relative deletion rate
$sanger_insert_default = 0; # Default for Metasim Sanger model relative insertion rate
$metasim_db_name_default = 'database'; # Name chosen by MetaSim for the directory containing its database
##### The effect of the default settings below (given other default settings) is to simulate 
##### a read-length distribution with an average of about 400bp and a standard deviation of about 80bp.
$mean_read_len_default=1200;          # Parameter for simulation: average read length; 
                                         # For better control of read length, especially for sampling fragments of 
                                         # gene sequences rather whole genomes, use the clone-length parameter to
                                         # control read length instead (i.e., set the read length to be much 
                                         # bigger than the clone length).
$stddev_read_len_default=60;          # Parameter for simulation: std dev of read length distribution
$mean_clone_len_default=400;          # Parameter for simulation: average length of clone fragment:
                                         # can be used to control read length directly if read-length parameter is
                                         # set to be much larger than clone length (i.e. the whole clone becomes a read).
$stddev_clone_len_default =80;        # Parameter for simulation: std dev of clone fragment length distribution
$seq_pad_size_default = 220;          # Sequences given to MetaSim can be padded by this 
                                         # number of 'N's on either end, to achieve boundary effects

# Miscellaneous file extensions
our (@fasta_alt_exts, $stockholm_aln_ext, $dna_fasta_ext, $gen_fasta_ext, $pep_fasta_ext, $aln_fasta_ext,
     $blastx_file_ext, $blastn_file_ext, $log_file_ext);
@fasta_alt_exts = qw(.fna .fa .fasta .pep);    # List of possible extensions for fasta files of unaligned sequences
$stockholm_aln_ext = '.sto';         # File extension for stockholm alignment, used by default by hmmalign
$dna_fasta_ext = '.fna';             # File extension to label files of unaligned nucleic acid sequences
$gen_fasta_ext = '.fa';              # File extension to label files of unaligned sequences in fasta format
$pep_fasta_ext = '.pep';             # File extesnion to label files of unaligned amino acid sequences
                                     #    (used by AMPHORA Reference Sequences files)
$aln_fasta_ext = '.aln';              # File extension to label fasta file of aligned peptide sequences
$blastx_file_ext = '.blx';            # File extension for output of blastx 
$blastn_file_ext = '.bln';            # File extension for output of blastn
# $blastp_file_ext = '.bp';            # File extension for output of blastp
$log_file_ext = '.log';              # File extension used to label general log file


# For maxPD.pl related stuff
our ($maxpd_refdb_default_num, $maxpd_default_num, $maxpd_margin_num, 
     $maxpd_filename, @former_names, %alt_names, @taxa_to_avoid);
$maxpd_margin_num = 25;     # Number of species to get from AMPHORA Genome tree via maxPD extra
                            #   so that the total can be whittled down to the number desired (removing taxa to avoid, etc).
# $maxpd_refdb_default_num = 75;    # Number of species gotten via maxPD to be used in reference database.
$maxpd_default_num = 578;   # Max number of species to calculate maxPD for by default
$maxpd_filename = 'maxPD.out';  # Part of the name of the output file if maxPD is called on a tree
# Changes from taxon names in the AMPHORA Reference Sequence headers (the keys) to the current NCBI name (value); 
#  also, alternative names (some are actually used in the AMPHORA Genome.tree though nowhere else in AMPHORA).
%alt_names = ('Silicibacter pomeroyi DSS-3' => 'Ruegeria pomeroyi DSS-3', 
	      'Thiomicrospira denitrificans ATCC 33889' => 'Sulfurimonas denitrificans DSM 1251',
	      'Lactobacillus salivarius UCC118' => 'Lactobacillus salivarius subsp salivarius UCC118',
              'Geobacter uraniumreducens Rf4' => 'Geobacter uraniireducens Rf4',
	      'Burkholderia cepacia AMMD' => 'Burkholderia ambifaria AMMD',
	      'Coxiella burnetii Dugway 7E9-12' => 'Coxiella burnetii Dugway 5J108-111',
	      'Salmonella enterica subsp. enterica serovar Typhi Ty2' => 'Salmonella enterica subsp enterica serovar Typhi str Ty2');
@taxa_to_avoid = qw(Mycoplasma Candidatus);  # Species having any of these words in their name will not be included
                                             # in the part of the ref db built via maxPD



# For blast specifically
our(@blast_db_pep_exts, @blast_db_na_exts, $blast_max_expect, $blast_max_ratio, $num_blastdb_hits_default);
@blast_db_pep_exts = ('.phr', '.pin', '.psq');   # The extensions used by blastall for the peptide database
@blast_db_na_exts = ('.nhr', '.nin', '.nsq');    # The extensions used by blastall for the nucleotide database
$num_blastdb_hits_default = 50;  # Gives the default number of blast hits found.
$blast_max_expect = 1;          # Assume all blastx or blastn expect values are smaller than $blast_max_expect 
                                # (used for parsing blast output); ignore reads that have only larger blast expect values;
                                # should not usually require customization.
$blast_max_ratio = 1e5;         # Translation is done in a frame only if the lowest blastx expect value for 
                                # that frame is 0 or is less than the product of blast_max_ratio and the 
                                # lowest blastx expect value over all frames; 
                                # similarly, a read is only reverse complemented if the lowest blastn expect 
                                # value for that strand is 0 or is less than the product of the blast_max_ratio
                                # and the lowest blastn expect value for both strands.


# alnqc: For quality control on the simulated reads and the alignment, specifically
our($seq_len_threshold_default, $qc_len_filename_part);
$seq_len_threshold_default = 25;   # The minimum length of a DNA sequence read to be retained 
$qc_len_filename_part = '-thr';  # Default part of filename added after short sequences are dropped


sub count_seqs ($) {
    
    my $file = shift(@_);
    my $header = '^>';
    open(IN, '<', $file) or die "Cannot open $file to count sequences.\n";
    my $count=0;
    while (my $line = <IN>) {
	if ($line=~ /$header/) {
	    $count++;
	}
    }
    close(IN);
    return $count;
}


# Returns the given file name -- without the path and without the extension, 
# if the extension is in the given array;
# if no array of extensions is given, @fasta_alt_exts is used by default.
sub get_basename ($@) {
    my $filename = shift(@_);
    my ($volume,$directories,$file) = File::Spec->splitpath( $filename );
    my @exts = @_;
    if (! @exts ) { @exts = @fasta_alt_exts; }
    foreach my $ext (@exts) {
	if ($file =~ ('(.+)'.$ext .'$')) {
	    return ($1);
	}
    }
    return $file;
}

sub check_ext ($;$) {
    my ($filename, $exts_ar) = @_;
    if (!$filename) {
	return;
    }
    my @exts = @fasta_alt_exts;
    if ($exts_ar and scalar(@{$exts_ar})) {
	push(@exts, @{$exts_ar});
    }
    my $real_ext;
    foreach my $ext (@exts) {
	if ($filename =~ ($ext .'$')) {
	    $real_ext = $ext;
	    last;
	}
    }
    return $real_ext;
}

# This function assumes the directory organization created by a default installation of AMPHORA.
sub define_amphora_paths(;$) {
    my $path = shift (@_);
    if (defined ($path)) {
	$amphora_path = $path;
    } elsif (!defined ($amphora_path)) {
	die "The amphora path variable is not yet defined.\n";
    }
    # path for AMPHORA Reference Sequences
    $amphora_ref_seq_path = File::Spec->catdir($amphora_path,'Reference/Sequence/');  
    # path for AMPHORA profile HMMs
    $amphora_marker_profile_path = File::Spec->catdir($amphora_path, 'Marker/');  
    $amphora_genome_tree = File::Spec->catfile($amphora_path,'Reference/Tree/Genome.tree');
    # directory containing the scripts MarkerScanner.pl and MarkerAlignTrim.pl
    $amphora_scripts_path = File::Spec->catdir($amphora_path,'Scripts/');  			 
}

1;
