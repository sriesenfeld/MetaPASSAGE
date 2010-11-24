#!/usr/bin/perl

# Main script for metagenomic simulation pipeline MetaPASSAGE.
# S.J. Riesenfeld
# Last updated:  Nov 2010
# Createded: April 2009 (roughly)
# See accompanying README.

=head2 MetaPASSAGE.pl

    MetaPASSAGE.pl is the main script for the metagenomic simulation
      pipeline.  See the accompanying README for more information and
      more usage examples.
    
  Usage: 

    To get information on all options:

      perl MetaPASSAGE.pl -h

    An example of how to run a complete workflow on a protein family
      of sequences:

      perl MetaPASSAGE.pl --seqs_file_basename <seqs> --num_seqs 100
        --sim --num_reads 1000 --blast --align --model <model>

    where <seqs>.fna and <seqs>.pep are existing files that contain
      corresponding sequences, DNA sequences in the first and amino
      acid sequences in the second; <model> is the name of a file
      containing an HMM profile for the sequences is <seqs>.pep.

    To run a similar workflow on a protein family in AMPHORA:
   
      perl MetaPASSAGE.pl --gene <gene_symbol> --num_seqs 100 --sim 
        --num_reads 1000 --blast --align

    where <gene_symbol> is the name used to identify a marker gene
      family in AMPHORA. 

    To run a similar workflow on a family of 16S rRNA sequences, first
      convert the sequences to DNA by running the script
      util_scripts/RNA2DNA.pl on the sequence file and then:

      perl MetaPASSAGE.pl --seqs_file_basename <seqs> --type RNA2DNA 
        --num_seqs 100 --sim --num_reads 1000 --blast --align 
        --model <model>

    where <seqs>.fna is an existing file that contains 16S DNA
      sequences (converted from RNA) and <model> is an existing file
      containing an INFERNAL profile model for these sequences.

  Reference Database Creation: 

    To create a reference database for use with option
      '--refdb_basename', see the README and the script
      createRefDB.pl.

=cut

use strict;
use warnings;

use File::Spec;
use File::Copy "mv";
use Getopt::Long;
use POSIX qw(INT_MAX);
use lib 'simPipeModules';

use simPipeVars qw(:seqtypes :amphora :metasim :misc_file_ext :simpipe :blast :alnqc);
use refSeqs qw(sample_seqs ); 
use seqsPad qw(seqs_pad seqs_remove_pad);
use runMetaSim qw(run_metasim create_profile filter_reads);
use blast2trans qw(run_blast_transeq);
use runScanAlignTrim qw(run_markerScannerAlt run_alignTrim);
use alnQC qw(drop_short_seqs);

require File::Temp;
use File::Temp;

my $usage = qq{
Usage: $0 <options (see below)>

  MetaPASSAGE.pl is the main script for the metagenomic simulation
    pipeline.  See the documenation in the source code and/or the
    accompanying README for more information and more usage examples.

  Minimally, options '--gene' or '--seqs_file_basename', '--num_seqs'
    and '--sim', (see below) would normally be chosen.  Options
    '--blast' and '--align' would also be standard.  Option
    '--blastdb' may be set if option '--blast' is used.

Options:

    -h, --help 
       [no value taken; if set, this help message is printed and
        program ends]

    --type 
       <string, giving data type: 'Protein' or 'RNA2DNA' (i.e., RNA
        converted to DNA); defaults to 'Protein'; (other types,
        e.g. 'DNA', could be added)>

    -g, --gene 
       <gene, family symbol, e.g., 'rpoB', specfying which protein
        family of AMPHORA is to be used; assumes a standard
        installation of AMPHORA; if not set and option '-y' is used
        then that filename is used instead of the gene symbol>

    --amphora_path 
       <the path to the AMPHORA directory; must be set if '--gene'
        option is used; the default, set in simPipeVars.pm, is
        '/usr/local/AMPHORA/'>

    -y, --seqs_file_basename 
       <optional base name (without the extension, but it should
        include the path) of fasta sequence files to use instead of
        AMPHORA Reference Sequences; there may be one file (DNA only)
        or two (DNA and amino acid) sequence files with the same
        basename: if sequences are DNA, filename extension assumed to
        be '.fna'; if sequences are peptides, filename extension
        assumed to be '.pep'; (set by $dna_fasta_ext and
        $pep_fasta_ext in simPipeVars.pm)>

     -c, --out_basename 
       <optional base name (without extensions) used in the naming
        scheme for most output files>

    -o, --out_dir 
       <optional name of output directory; excluding the blast
        database files and the reference database, output files will
        be written in this directory.>

    -v, --refdb_basename 
       <optional base name (without extensions, but it should include
        the path) of reference database files, which should already
        exist; if type is 'protein', peptide sequence file with
        extension '.pep' must exist; if type is 'RNA2DNA', DNA
        sequence file with extension '.fna' must exist; (these default
        extensions are in the module simPipeVars.pm); if no file is
        found with the given base name, then the reads are aligned
        without the reference sequences; used with option '-a', and
        also with option '-m' and '-d'>

    -n, --num_seqs 
       <total number of distinct sequences to be sampled for the
        simulated population from the whole set of full-length
        sequences (except Candidatus and Mycoplasma sequences), which
        is specified by option '-g' or '-y'; see option '-m', which
        affects the extent to which these sequences are sampled from
        the reference database; the sampled sequences are written to a
        fasta file; if option '-s' is set, the current MetaSim
        database is deleted and these sequences are automatically
        added to the MetaSim database; by default (if neither option
        '-t' nor option '--taxon_profile_ratio' is set), a uniform
        MetaSim taxonomic profile (with extension '.mprf') for these
        sequences is generated>

    -m, --num_ref_seqs
        <number of sequences among the total specified by option '-n'
        that are sampled uniformly at random from the reference
        database; the remainder are sampled uniformly at random from
        the non-reference sequences in the set of all sequences
        (specified by option '-g' or option '-y'); if this option is
        set to a nonnegative number (even to 0), then the reference
        database must exist; if this option is set to a negative
        number, then the reference database is ignored during
        sampling, i.e., the number of sequeces specified by option
        '-n' is sampled from among all input sequences; by default,
        the option is set to -1>

    -i, --src_file 
       <(alternate to option '-n') fasta file of DNA sequences to be
        used to simulate metagenomic reads via option '-s' or for
        parsing via option '-q' a file of metagenomic reads specified
        by option '-f'; use option '-j' to add the sequences to the
        MetaSim database if they have not already been added; option
        '-i' is ignored if option '-n' is set; if option '-s' is set
        and '-t' is not, then a uniform MetaSim taxonomic profile for
        these sequences is generated.>

    --pad
       <integer, size of pad of 'N's that will be concatenated to each
        end of each full-length sequence before the full-length
        sequences are added to the metasim database and used to
        generate reads; this is used to create a more realistic length
        distribution over the reads for a gene family and to get rid
        of boundary effects; by default, it is set to the default
        value of the mean clone length (in simPipeVars.pm); if it is
        set to 0, no padding is done>

    -j, --add_seqs 
       [no value taken; if set, add the sequences specified by option
	'-i' to the MetaSim database (after deleting the currently
	existing database); ignored if option '-i' is not set]

    --tax_profile_ratio
       <a number that is used to generate a skewed distribution for
        the relative abundances of sequences in the simulated
        population; the output is a taxonomic profile with extension
        '.mprf' for MetaSim; technically, the value of this option is
        the ratio between the second-to-smallest relative abundance
        and the minimum relative abundance, which is assumed to be 1,
        of distinct sequences in the simulated community; the
        abundances are proportional to a geometric distribution for a
        subset of the sequences (either half the sequences or until
        the maximum relative abundance is reached); the remaining
        sequences have relative abundance 1, to elongate the tail of
        the whole distribution; by default, the ratio is set to 1, so
        the whole distribution is uniform>

    --tax_profile_max_abundance
       <an integer designating the maximum relative abundance in the
        taxonomic profile generated for use with MetaSim, where the
        minimum relative abundance is assumed to be 1; by default, the
        maximum integer value is used>

    -t, --tax_profile 
       <optional .mprf taxonomic profile file (described in
        documentation for MetaSim) to be used in simulation via option
        '-s'; sequences referenced in this profile file should be
        added previously or in this command to the MetaSim database.>

    --metasim_path 
       <optional string giving the location of the metasim
        installation directory or wherever MetaSim writes its
        error.log file by default; used to record the information
        MetaSim outputs about each simulation run; default, set in
        simPipeVars.pm, is '/usr/local/metasim'; the functionality is
        only available if the path exists>

    -s, --sim 
       [no value taken; if set, runs MetaSim simulation with either an
        automatically generated profile file, based on sequences
        chosen according to option '-n' or '-i', or a profile file
        specified by option '-t'; by default, uses a Sanger-based
        sequencing model with 454-like read-length distribution and no
        sequencing error (Sanger is used because it is easy to control
        from the command line and offers a default of no sequencing
        error); defaults set in simPipeVars.pm; see options below for
        modifying aspects of this model; see option
        '--ms_error_conf_file' for specifying a different empirical
        model]

    -r, --num_reads 
       <optional number of reads to simulate with MetaSim; default is
        $num_reads_default (set in simPipeVars.pm); ignored if option
        '-s' is not set>

    -p, --mean_read_length 
       <optional mean length of read to simulate with MetaSim; default
        is $mean_read_len_default (please see simPipeVars.pm, where
        this is set); passed to MetaSim. Note: it may make sense in
        the simulations to control the read length via the clone
        length, i.e., to set the read length to be very large and the
        clone length to be what you want the read length to be>

    -u, --stddev_read_length 
       <optional standard deviation of read length distribution;
        default is $stddev_read_len_default (please see
        simPipeVars.pm, where this is set); passed to MetaSim>

    --cl, --mean_clone_length
       <optional mean length of clone fragment for simulation with
        MetaSim; default is $mean_clone_len_default (please see
        simPipeVars.pm, where this is set); passed to MetaSim. Note:
        it may make sense in the simulations to control the read
        length via the clone length, i.e., to set the read length to
        be very large and the clone length to be what you want the
        read length to be>

    --cs, --stddev_clone_len
       <optional standard/max deviation of clone fragment length
        distribution; default is $stddev_clone_len_default (please see
        simPipeVars.pm, where this is set); passed to MetaSim>

    --sanger_err_start
       <value is the initial error rate for the MetaSim Sanger error
        model; passed to MetaSim; default is 0>

    --sanger_err_end
       <value the final error rate for the MetaSim Sanger error
        model; passed to MetaSim; default is 0>

    --sanger_deletions
       <value is the relative deletion rate for the MetaSim Sanger
        error model; passed to MetaSim; default is 0>

    --sanger_insertions
       <value is the relative insertion rate for the MetaSim Sanger
        error model; passed to MetaSim; default is 0>

    --sanger_mate_probability 
       <value is the probability of paired reads for the MetaSim
        Sanger error model; passed to MetaSim; default is 0>

   --ms_error_conf_file
       <a file that can be used with MetaSim's Empirical sequencing
       model to simulate reads with various kinds of sequencing error;
       the default sequencing model is Sanger without any error>

   --mn, --no_delete_metasim_db
       [no value taken; if set, the metasim database is not deleted
        after reads have been generated; by default, the database is
        deleted because it takes space to store and, if a single
        database is stored for multiple runs, different calls to
        metasim may interfere with each other]

   --mg, --global_metasim_db
       [no value taken, if set, the metasim is called from the
        calling directory rather than from the output directory (and
        hence, a single database is stored in the calling directory
        rather than one in each output directory, for multiple runs);
        by default, metasim is run in the output directory so that
        different calls to metasim do not interfere with each other]

    -f, --reads_file 
       <(alternate to option '-s') fasta file of metagenomic reads to
        be blasted and translated as per option '-x'; ignored if '-s'
        is set>

    -q, --target_num_reads 
      <optional number of reads written into a new file of simulated
       metagenomic reads, which are a subset of those reads that have
       been created via option '-s' or specified via option '-f'; if
       set to a negative value, chooses the subset so that there is at
       most one read per sequence for each sequence in the fasta file
       of complete gene sequences created by option '-n' or specified
       via option '-i'; if set to a positive value, chooses a subset
       of the specified size uniformly at random from the full set of
       reads (which should preserve the skew of the population
       distribution, if there is one); if set to 0 or not set, no
       filtering is done; the reason one might use this sampling
       filter is to generate more reads than will be used, as not all
       reads will necessarily be recognized by BLAST (or AMPHORA) as
       belonging to the gene family; if '-x' is set, filtering is done
       after reads have been translated/oriented (and possibly scanned
       via AMPHORA), which means one should be careful about how many
       reads are initially created, since the BLAST step and AMPHORA
       scan are slow>

    -d, --blastdb 
       <the base name (including the path) of an existing blast
        database, for use with option '-x'; if the files are not
        found, ncbi-blast formatdb is run to create a blast database
        of that name in the current directory (to prevent errors from
        formatdb, make sure the path of directories exists even if the
        blast database files do not yet exist); this db is created
        either from the reference database, if it is specified by
        option '-v' and exists, or from AMPHORA Reference Sequences
        (peptides) for the gene specified via the '-g' option, or from
        the sequences in the file specified via the '-y' option; if
        option '-x' is used without option '-d', the database is
        assumed to have the default base name of the reference
        database, if it is specified by option '-v' and exists, or of
        the given gene symbol specified by option '-g', or the name of
        the file specified by option '-y', and if it does not exist,
        it is created>

    -x, --blast, --blastx, --blastn
       [no value taken; if set, and if type is PROTEIN, run blastx on
	metagenomic reads against blast database to find frame for
	translation and translate sequences into peptides; if set, and
	if type is RNA2DNA, run blastn on metagenomic reads against
	DNA blast database to find orientation; blast database may
	be previously built or getting built via option '-d'; by
	default (if option '-f' is not set), the metagenomic reads
	used are those generated according to option '-s'; option '-s'
	or '-f' is mandatory if option '-x' is used; if option '-q' is
	also set, the file of single metagenomic reads written by that
	option setting is used instead.]

    -k, --save_blast_out
       [no value taken; if set, the (large) blast files used for
        translation/orientation are saved; by default, they are
        discarded]

    -e, --blast_expect
       <optional expectation value cutoff for blastx; default chosen
        by blastall is 10.0>

    -b, --num_blastdb_hits
       <optional number of database sequences to show scores for from
        blastx; default is 50>

    -a, --align
       [no value taken; if set and an AMPHORA gene family is specified
        with option '-g', the peptide reads (created via option '-x'
        or input via option '-l') are scanned with a modified version
        of the AMPHORA script MarkerScanner.pl and aligned, along with
        the reference database sequences (if given), to the AMPHORA
        HMM 2 profile for the specified gene, using the AMPHORA script
        MarkerAlignTrim.pl (see option '--hmmer2_path'); if set and
        other sequences have been input via option '-y', then the
        reads, which have been oriented or translated via option '-x'
        or input via option '-l', are aligned, along with the
        reference database sequences (if specified via option '-v'),
        to the model in the file provided via option '--model'; the
        alignment method used is INFERNAL cmalign if the sequence
        type, specified via option '--type', is DNA converted from RNA
        ('RNA2DNA'), and the method is HMMER hmmalign if the sequence
        type is amino acid ('Protein'); the INFERNAL or HMMER commands
        are expected to be in the user's path]

    --align_source
       [no value taken; if set, then the full-length sequences that
        correspond to the final set of reads are aligned, along with
        the reference database (if given), to the model in the file
        provided via option '--model'; see option '-a' for details on
        how the alignment method is decided; this option is ignored if
        option '-a' is not also selected]

    --save_align_log
       [no value taken; if set, then a log of the alignment process is
        saved; this can be useful when aligning with AMPHORA; default
        is not to save the file]

    --hmmer2_path
       <the path to HMMER 2 binaries, in particular, hmmpfam; used
        only with AMPHORA; default set in simPipeVars.pm; see also
        option '-a'>

    -l, --tr_reads_file 
       <optional fasta file of peptide translations of simulated
        metagenomic reads, if type is PROTEIN, or of correctly
        oriented DNA reads, to be aligned via option '-a'>

    --drop_len
       <an integer giving the minimum length (in amino acids, if type
        is PROTEIN, and in base pairs otherwise) of a read to be
        included in the output alignment; reads generated that are
        shorter than this length will be excluded (this only applies
        if an alignment is being constructed); default value is 25
        amino acids / 75 base pairs (set in simPipeVars.pm); if set to
        0, no sequences are dropped>

    --model
       <filename for a file containing the model used for alignment;
        this option is ignored if AMPHORA is being used to align
        sequences; if option '-y' is used and not option '-g', then
        either cmalign or hmmalign will be used to do the alignment
        (depending on the type), and a model must be provided; see
        option '-a'>

    -z, --noscan 
       [no value taken; if set, scanning will be skipped during the
        alignment via option '-a']

    -w, --log
       [no value taken; if this option is used, the output that is by
        default written to the terminal is instead written to a log
        file.]

};

my $error_msg = qq{
    Check input options, as well as locally customized variables in simPipeVars.pm.
    };

my ($help_flag, $add_seqs_flag, $sim_flag, $blast_flag, $align_flag,
    $align_src_flag, $noscan_flag, $save_blast_out_flag, $log_flag,
    $log_fh, $log_file, $metasim_log, $sanger_err_start,
    $sanger_err_end, $sanger_insertions, $sanger_deletions,
    $sanger_mate_probability, $metasim_error_conf_file,
    $gene_symbol, $amphora_path_nondefault, $alt_gene_symbol,
    $num_seqs, $num_ref_seqs, $type, $cur_dir, $out_dir,
    $full_path_cur_dir, $basename, $fullpath_basename, $use_refdb,
    $ref_db_file_basename, $ref_db_file_pep, $ref_db_file_dna,
    $alt_seqs_basename, $sample_file_basename, $sample_file_dna,
    $sample_file_dna_padded,$sample_file_pep, $pad_size,$src_seqs_dna,
    $sample_final_src_seqs, $final_num_src_seqs, $tax_profile_ratio,
    $tax_profile_max_abundance, $metasim_profile, $metasim_reads_file,
    $dna_reads_file, $thr_metasim_reads_file, $scanned_pep_reads_file,
    $filtered_sim_reads_file, $reads_to_align, $num_reads,
    $num_filtered_reads, $final_num_reads, $no_del_metasim_db_flag,
    $global_metasim_db_flag, $mean_read_len, $stddev_read_len,
    $mean_clone_len, $stddev_clone_len, $num_blastdb_hits,
    $blast_expect_cutoff, $blastdb_name, $blastdb_ne_flag,
    $frames_file, $trans_file, $tr_reads_file, $hmm_profile, $model,
    $len_threshold, $align_log_flag, @args);

####### SET-UP DEFAULT VALUES FOR SOME VARIABLES

my $default_type = PROTEIN;  # default type is PROTEIN
$basename = $sim_basename_default;
$cur_dir = File::Spec->curdir();
$out_dir = $cur_dir;
$full_path_cur_dir = File::Spec->rel2abs($cur_dir);
$num_reads=$num_reads_default; 
$global_metasim_db_flag = 0;
$no_del_metasim_db_flag = 0;
$num_ref_seqs = -1;  #ignore ref db when sampling, i.e., sample from any of the sequences
$pad_size = $seq_pad_size_default;
$tax_profile_ratio = 1;
$tax_profile_max_abundance = INT_MAX;
$num_filtered_reads = 0;
$num_blastdb_hits = $num_blastdb_hits_default;
$sanger_mate_probability = $sanger_mate_prob_default;
$sanger_insertions = $sanger_insert_default;
$sanger_deletions = $sanger_delete_default;
$sanger_err_start = $sanger_err_start_default;
$sanger_err_end = $sanger_err_end_default;
$align_log_flag = 0;

####### GET OPTION SETTINGS
# using single letters as alternate optional names so that this version is compatible with previous versions
GetOptions(
    'help|h|?' => \$help_flag,
    'type:s' => \$type,
    'gene|g:s' => \$gene_symbol,
    'amphora_path:s' =>\$amphora_path_nondefault,
    'out_basename|c:s' => \$basename,
    'out_dir|o:s' => \$out_dir,
    'seqs_file_basename|y:s' => \$alt_seqs_basename,
    'refdb_basename|v:s' => \$ref_db_file_basename,
    'log|w' => \$log_flag,
    'num_seqs|n:i' => \$num_seqs,
    'num_ref_seqs|m:i' => \$num_ref_seqs,
    'src_file|i:s' => \$sample_file_dna,
    'pad:i' => \$pad_size,
    'add_seqs|j' => \$add_seqs_flag,
    'tax_profile_ratio:f' => \$tax_profile_ratio,
    'tax_profile_max_abundance:i', => \$tax_profile_max_abundance,
    'tax_profile|t:s' => \$metasim_profile,
    'sim|s' => \$sim_flag,
    'num_reads|r:i' => \$num_reads,
    'mean_read_len|p:i' => \$mean_read_len,
    'stddev_read_len|u:i' => \$stddev_read_len,
    'mean_clone_len|cl:i' => \$mean_clone_len,
    'stddev_clone_len|cs:i' => \$stddev_clone_len,
    'sanger_err_start:f' => \$sanger_err_start,
    'sanger_err_end:f'=> \$sanger_err_end,
    'sanger_insertions:f' => \$sanger_insertions,
    'sanger_deletions:f' => \$sanger_deletions,
    'sanger_mate_probability:f' => \$sanger_mate_probability,
    'reads_file|f:s' =>	\$metasim_reads_file,
    'global_metasim_db|mg' => \$global_metasim_db_flag,
    'no_delete_metasim_db|mn' => \$no_del_metasim_db_flag,
    'target_num_reads|q:i' => \$num_filtered_reads,
    'ms_error_conf_file:s' => \$metasim_error_conf_file,
    'metasim_path:s' => \$metasim_path,
    'blastdb|d:s' => \$blastdb_name,
    'blast|blastx|blastn|x' => \$blast_flag,
    'save_blast_out|k' => \$save_blast_out_flag,
    'num_blastdb_hits|b:i' => \$num_blastdb_hits,
    'blast_expect|e:f' => \$blast_expect_cutoff,
    'tr_reads_file|l:s' => \$tr_reads_file,
    'align|a' => \$align_flag,
    'align_source' => \$align_src_flag,
    'save_align_log' => \$align_log_flag,
    'hmmer2_path:s' => \$hmmer2_path,
    'drop_len:i' => \$len_threshold,
    'model:s' => \$model,
    'noscan|z' => \$noscan_flag    
    );


if ($help_flag) {
    die "$usage\n";
}

if (defined ($amphora_path_nondefault)) {
    $amphora_path = File::Spec->rel2abs($amphora_path_nondefault);
    define_amphora_paths($amphora_path);
}

####### FILENAMES & MISC SET-UP

unless ( (-d $out_dir) and (-w $out_dir) ) {
    mkdir($out_dir) or die "Cannot find, write, or make directory $out_dir.\n";
}

$fullpath_basename = File::Spec->catfile($out_dir, $basename);
$fullpath_basename = File::Spec->rel2abs($fullpath_basename);

if ($log_flag) {
    $log_file = $fullpath_basename.$log_file_ext;
    open(LOG, ">$log_file") or die "Cannot open log file $log_file: $!\n";
    print "\nInformation about this run of MetaPASSAGE.pl is being logged in file ".
	"\n  ".File::Spec->abs2rel($log_file)."\n";
    $log_fh = \*LOG;
} else { $log_fh = \*STDOUT; }

$metasim_log = $fullpath_basename.'-'.$metasim_log_filename;

if (! defined($type)) {
    $type= $default_type;
    print $log_fh "Assuming default sequence data type of $type.\n";
} elsif (!scalar( grep {$_ eq $type} @types )) {
    die "Type $type is not defined! Please use one of the following types: \n".
	join(', ', @types)."\n".
	"(Default is ".$default_type.".)\n";
} 

if (! defined($len_threshold)) {
    $len_threshold = $seq_len_threshold_default;
}

if (! defined($gene_symbol) ) {
    if (defined($alt_seqs_basename)) {
	# use basename of file of alternative input sequences as gene symbol
	my ($volume,$directories,$file) = File::Spec->splitpath( $alt_seqs_basename );
	$alt_gene_symbol = get_basename($file, (@fasta_alt_exts, $stockholm_aln_ext, $aln_fasta_ext));
	# print $log_fh "No gene symbol provided via option '-g'; using $alt_gene_symbol.\n";
    } 	
} 

if (defined ($ref_db_file_basename) ) {
    $ref_db_file_pep = $ref_db_file_basename.$pep_fasta_ext;
    $ref_db_file_dna = $ref_db_file_basename.$dna_fasta_ext;
} 

if ( ( ($type eq PROTEIN) and defined ($ref_db_file_pep) and ((-e $ref_db_file_pep)) )
     or ( ($type ne PROTEIN) and defined ($ref_db_file_dna) and ((-e $ref_db_file_dna)) ) ){
    $use_refdb=1;
} else {
    if (defined ($ref_db_file_basename)) {
	print $log_fh "Cannot find reference database file ". 
	    ( ($type eq PROTEIN) ? $ref_db_file_pep : $ref_db_file_dna ). "!\n". 
	    "Trying to run the pipeline without it.\n";
    }
    $use_refdb=0;
}

####### LOOK FOR PRE-EXISTING BLAST DATABASE
my @blast_db_exts = ($type eq PROTEIN) ? @blast_db_pep_exts : @blast_db_na_exts;
$blastdb_ne_flag=0;    
if ($blast_flag or defined($blastdb_name)) {
    if (! defined ($blastdb_name)) {
	if ($use_refdb) {
	    $blastdb_name = $ref_db_file_basename;
	} elsif (defined($gene_symbol) or defined($alt_gene_symbol)) {
	    $blastdb_name = File::Spec->catfile($cur_dir, defined($gene_symbol) ? $gene_symbol : $alt_gene_symbol);
	} else {
	    die "Please specify the name of a blast database via option '-d' or a gene symbol via option '-g'\n".
		"(for use with AMPHORA) or a file of sequences via option '-y' that can be \n".
		"used to format a blast database.\n";
	}
    }
    foreach my $ext (@blast_db_exts) {
	if (!(-e $blastdb_name.$ext)) {
	    $blastdb_ne_flag=1;
	    last;
	}
    }
    if ( $blastdb_ne_flag && (! defined ($gene_symbol)) && (! defined ($alt_gene_symbol)) ) {
	die "Cannot find all blast database files with basename $blastdb_name.\n".
	    "Please specify the name of an existing blast database via option '-d'!\n";
    }
}

####### SAMPLE FULL-LENGTH SEQUENCES FOR POPULATION

if ($num_reads < 0) {
    die "Please specify a positive number of reads!\n";
}

if (defined ($num_seqs) and ($num_seqs < 0)) {
    die "Please enter a positive total number of full-length sequences to sample.\n";
} elsif ($num_seqs) {
    if (defined($num_ref_seqs) and ($num_ref_seqs >=0)) {
	if (!defined ($ref_db_file_basename)) {
	    die "If option '-m' is set to be nonnegative, then a reference database must be specified via option '-v'.\n";
	}
	if (!$use_refdb) {
	    die "Cannot find reference database ".( ($type eq PROTEIN) ? $ref_db_file_pep : $ref_db_file_dna )."!\n".
		"  If option '-m' is set to be nonnegative, you must specify an existing reference database\n".
		"  via option '-v' for sampling purposes.\n";
	} 
	if ($num_ref_seqs > $num_seqs) {
	    die "The number of sequences to sample from the reference database must be at most the total number of full-length sequences to be sampled.\n";
	}
    }
    # set files names of sampled sequences and related files 
    $sample_file_basename = $fullpath_basename.'-'.$src_filename_default;
    $sample_file_pep = $sample_file_basename.$pep_fasta_ext;
    $sample_file_dna = $sample_file_basename.$dna_fasta_ext;
    print $log_fh "Sampling $num_seqs full-length sequences";
    if (defined($num_ref_seqs) and ($num_ref_seqs >=0)) {
	print $log_fh " ($num_ref_seqs reference sequences).\n";
    } else {
	print $log_fh ".\n";
    }
    @args = ($num_seqs, $sample_file_basename);
    if (defined($alt_seqs_basename)) {
	push(@args, $alt_seqs_basename);
	if ($num_ref_seqs >= 0) {
	    push(@args, $num_ref_seqs, $ref_db_file_dna);
	}
    } else {
	if (! defined ($gene_symbol)) {
	    die "Please use option '-g' to provide a gene symbol (for use with AMPHORA) or option '-y'\n"
		." to provide a file of sequences in order to generate reads.\n";
	}
	my $amphora_ref_seqs_basename = File::Spec->catfile($amphora_ref_seq_path, $gene_symbol);
	push(@args, $amphora_ref_seqs_basename);
	if ($num_ref_seqs >= 0) {
	    push(@args, $num_ref_seqs, $ref_db_file_pep);
	}
    }
    sample_seqs(@args);    
    if (-e $sample_file_pep) { 
	print $log_fh "Peptide sequences sampled from AMPHORA ".
	    "and written to file ".File::Spec->abs2rel($sample_file_pep).
	    ";\n";
    } elsif (! defined($alt_seqs_basename)) {
	die "Trouble sampling sequences and writing to file ". File::Spec->abs2rel($sample_file_pep)."!\n";
    }
    if (-e $sample_file_dna) {
	print $log_fh "DNA sequences have been sampled and written to ".
	    File::Spec->abs2rel($sample_file_dna).".\n";
    } else {
	die "Trouble sampling and writing sequences to file!\n";
    }    
} elsif (defined ($sample_file_dna)) {  
    if (-e $sample_file_dna) {
	foreach my $ext (@fasta_alt_exts) {
	    if ($sample_file_dna =~ ('(.+)'.$ext .'$')) {
		$sample_file_pep = $1.$pep_fasta_ext;
		last;
	    }
	}
	if (!(-e $sample_file_pep)) {
	    $sample_file_pep = undef;
	}
    } else {
	die "Cannot find file $sample_file_dna containing sequences for metagenomic simulation.\n";
    }
}

####### GENERATE SIMULATED READS USING METASIM
if (($sim_flag) and !(-d $metasim_path)) {
    print $log_fh "Note: Cannot find metasim installation directory here $metasim_path for error log copy. Ignoring.\n";
}

if ( (defined ($sample_file_dna) and $add_seqs_flag) or ($num_seqs and $sim_flag) ) {
    # pad sequences
    if ($pad_size) {
	($sample_file_dna_padded, $pad_size)=seqs_pad($sample_file_dna, $pad_size);
    }
    print $log_fh "Full-length sequences padded on each end with $pad_size \'N\'s were written to file".
	"\n  ". File::Spec->abs2rel($sample_file_dna_padded)."\n";
    print $log_fh "Adding sequences from ". File::Spec->abs2rel($sample_file_dna) ." to the MetaSim database.\n";
    $src_seqs_dna = defined($sample_file_dna_padded) ? $sample_file_dna_padded : $sample_file_dna;    

    # add sampled full-length sequences to the MetaSim database
    if (! $global_metasim_db_flag) {
	# default behavior: run metasim from the output directory (because metasim stores its database in calling dir)
	chdir($out_dir) || die "Cannot change directories to output directory: $out_dir: $!\n";
    }
    
    run_metasim('-a', $src_seqs_dna, '-f', $metasim_log);   
    if (! $global_metasim_db_flag) {
	chdir($full_path_cur_dir) || die "Cannot change directories back to original directory: $full_path_cur_dir: $!\n";
    }    
}
if ( !(defined($metasim_profile)) and ($num_seqs or defined($src_seqs_dna) ) 
     and $sim_flag ) {
    # create taxonomic profile
    $metasim_profile = create_profile('-i', $src_seqs_dna, '-o', $out_dir, '-k', $tax_profile_ratio, 
				      '-a', $tax_profile_max_abundance);
    if (! (-e $metasim_profile)) {
	die "Error creating taxonomic MetaSim profile $metasim_profile.\n";
    }
    print $log_fh "Taxonomic profile for MetaSim written to $metasim_profile.\n";    
} elsif (defined($metasim_profile)) {
    # or use given taxonomic profile
    if (! (-e $metasim_profile)) {
	die "Cannot find taxonomic profile $metasim_profile for MetaSim\n";
    }
}
if ($sim_flag) {
    if (! defined ($metasim_profile) ){
	die "A taxonomic profile for MetaSim must be generated or specified via option '-n', '-i', or '-t'.\n";
    }
    # simulate metagenomic reads using determined taxonomic profile
    if (! (-e $metasim_profile)) {
	die "Cannot find taxonomic profile $metasim_profile for MetaSim.\n";
    }
    $metasim_profile = File::Spec->rel2abs($metasim_profile);
    print $log_fh "Running MetaSim to create $num_reads reads.\n";
    $metasim_reads_file = get_basename($metasim_profile, $metasim_profile_ext);
    if ($metasim_reads_file =~ /(.*)-$src_filename_default$/) {
	$metasim_reads_file = $1;
    }    
    $metasim_reads_file = $metasim_reads_file.$metasim_filename_default.$dna_fasta_ext;
    @args = ('-m',  $metasim_profile, '-r', $num_reads, '-f', $metasim_log, 
	     '-o', $metasim_reads_file);
    if (defined ($metasim_error_conf_file)) {
	$metasim_error_conf_file = File::Spec->rel2abs($metasim_error_conf_file);
	push(@args, '-e', $metasim_error_conf_file);
    } else {
	push(@args, '-l', $mean_read_len, '-s', $stddev_read_len, 
	     '-c', $mean_clone_len, '-p', $stddev_clone_len, 
	     '--mate_prob', $sanger_mate_probability,
	     '--err_insert', $sanger_insertions, '--err_delete', $sanger_deletions,
	     '--err_start', $sanger_err_start, '-err_end', $sanger_err_end);
    }
    unless( $no_del_metasim_db_flag) {
	# default behavior: metasim database will be deleted after reads are generated
	push(@args, '-x', 1);
    }
    if (! $global_metasim_db_flag) {
	# default behavior: run metasim from the output directory 
        # (because metasim stores its database in calling directory) and multiple runs interfere with each other
	chdir($out_dir) || die "Cannot change directories to output directory: $out_dir: $!\n";
    } else {
	push(@args, '-d', $out_dir);
    }
    $dna_reads_file=run_metasim(@args);    
    my $empty_seqs;
    if ($pad_size) {
	($metasim_reads_file, $empty_seqs) = seqs_remove_pad($dna_reads_file);
	print $log_fh "$empty_seqs read". (($empty_seqs == 1) ? '' : 's').
	    " dropped from the set of metasim reads due to malformation after padding was removed.\n";
    } else {
	$metasim_reads_file = $dna_reads_file;
    }
    unless ( $global_metasim_db_flag) {
	$metasim_reads_file = File::Spec->rel2abs($metasim_reads_file);
	chdir($full_path_cur_dir) || die "Cannot change directories back to original directory: $full_path_cur_dir: $!\n";
	$metasim_reads_file = File::Spec->abs2rel($metasim_reads_file);
    }
    print $log_fh "MetaSim run information written to ".File::Spec->abs2rel($metasim_log)."\n";
    if (! (-e $metasim_reads_file)) {
	die "Error creating reads file $metasim_reads_file!\n";
    }
    print $log_fh "Simulated reads written to file ". File::Spec->abs2rel($metasim_reads_file).".\n";
    unless ($no_del_metasim_db_flag) {
	print $log_fh "MetaSim database deleted.\n";
    }

} elsif (defined ($metasim_reads_file) and !(-e $metasim_reads_file)) {
    # use specified file of metagenomic reads
    die "Cannot find file $metasim_reads_file containing simulated metagenomic reads.\n";
}

##### DROP SHORT READS NOW

if (defined ($metasim_reads_file)) {
    my $ndrop;
    if ($len_threshold > 0) {
	($thr_metasim_reads_file, $ndrop) = 
	    drop_short_seqs($metasim_reads_file, $len_threshold);
	print $log_fh "Dropped the $ndrop read".(($ndrop ==1) ? '' : 's')." shorter than $len_threshold base pairs;\n".
	    "  the rest have been written to $thr_metasim_reads_file.\n";    
    } else {
	$thr_metasim_reads_file = $metasim_reads_file;
    }
}
####### POSSIBLY FILTER NOW, IF NO TRANSLATION OR ALIGNMENT IS BEING DONE,
#######    ACCORDING TO TARGET NUMBER DESIRED;
####### ALSO CREATE A FILE OF FULL-LENGTH SEQUENCES CORRESPONDING TO THE FINAL SET OF READS 

if ((!$align_flag) and (!$blast_flag) and ($num_seqs or defined($sample_file_dna)) and 
    (defined($thr_metasim_reads_file) or defined($tr_reads_file)) ) {
    my $reads_to_filter = defined($thr_metasim_reads_file) ? $thr_metasim_reads_file : $tr_reads_file;
    if ($num_filtered_reads < 0) {
	print $log_fh "Filtering the reads so that there is at most one read per original sampled sequence.\n";
    } elsif ($num_filtered_reads > 0) {	
	print $log_fh "Filtering the reads:  choosing a subset of size $num_filtered_reads uniformly at random.\n";
    }   
    ($filtered_sim_reads_file, $sample_final_src_seqs, $final_num_reads, $final_num_src_seqs) =
	filter_reads($reads_to_filter, $num_filtered_reads, $sample_file_dna, $out_dir, 
		     ($type eq PROTEIN) ? $sample_file_pep : $sample_file_dna); 
    print $log_fh "Determining the full-length sequences that correspond to the final set of reads\n".
    print $log_fh "Each read id converted to include the id for the corresponding full-length sequence.\n";
    print $log_fh ''.$final_num_reads. " reads written to ". File::Spec->abs2rel($filtered_sim_reads_file).".\n";
    if (-e $sample_final_src_seqs) {
	print $log_fh ''.$final_num_src_seqs.
	    " full-length sequences corresponding to these reads written to file ". $sample_final_src_seqs."\n";
    }
} 

####### FORMAT BLAST DATABASE (IF NECESSARY)

if ($blastdb_ne_flag) {
    print $log_fh "Cannot find at least one of the blast database files with base name $blastdb_name.\n".
	"  Building a new blast database.\n";
    my $blast_seqs_file;
    if ($use_refdb) {
	$blast_seqs_file = ($type eq PROTEIN) ? $ref_db_file_pep : $ref_db_file_dna;
    } else {
	$blast_seqs_file= (defined ($alt_seqs_basename)?  $alt_seqs_basename : 
			   File::Spec->catfile($amphora_ref_seq_path, $gene_symbol)) . 
			   (($type eq PROTEIN) ? $pep_fasta_ext : $dna_fasta_ext);
    }
    print $log_fh "Formatting blast database $blastdb_name from $blast_seqs_file.\n";
    my $cmd = 'formatdb -i '. $blast_seqs_file .' -n '. $blastdb_name . (($type ne PROTEIN) ? ' -p F' : '');
    eval {
	system ($cmd);
    }; 
    if ($@) {
	die "Error running formatdb for blast database: $@\n";
    }
    foreach my $ext (@blast_db_exts) {
	if (!(-e $blastdb_name.$ext)) {
	    die "Error running formatdb for blast database: Database files not found!\n"
		."Make sure path of directories exists for specified blast database basename:\n".
		$blastdb_name."\n";
	}
    }
} 

####### BLAST READS AGAINST BLAST DATABASE (FORMATTED FROM FULL-LENGTH
####### SEQUENCES) TO ORIENT OR TRANSLATE READS

if ($blast_flag) {
    my @blastArgs = ('-d', $blastdb_name);
    if (defined($num_blastdb_hits)) {
	push(@blastArgs, '-b', $num_blastdb_hits);
    }
    if (defined($blast_expect_cutoff)) {
	push(@blastArgs, '-e', $blast_expect_cutoff);
    }
    if ( !defined ($thr_metasim_reads_file) ) {
	die "Option '-s' or '-f' must be used to specify a file of simulated metagenomic reads to be blasted.\n";
    }
    if ($type eq PROTEIN) {
	print $log_fh "Running blastx on reads to figure out which frames to translate in,\n".
	    "  and then translating reads into peptides -- this may take some time.\n";
    } else {
	print $log_fh "Running blastn on reads to figure out how to orient them \n".
	    "  (i.e. if they should be reverse-complemented),\n".
	    "  and then writing correctly oriented reads to file -- this may take some time.\n";
    }
    my $blast_infile = $thr_metasim_reads_file;
    my $blast_outfile;    
    # the $frames_file has orientation information if $type ne PROTEIN;
    # $trans_file has correctly oriented reads if $type ne PROTEIN.
    ($trans_file, $frames_file, $blast_outfile) = run_blast_transeq($blast_infile, $out_dir, \@blastArgs, 
								    $type, $save_blast_out_flag);
    if ( !(-e $trans_file) or !(-e $frames_file) or ( $save_blast_out_flag and !(-e $blast_outfile)) )  {
	die "Error blasting and translating $blast_infile.\n";
    }
    if ($save_blast_out_flag) {
	print $log_fh "Blast output written in file ". File::Spec->abs2rel($blast_outfile).".\n";
    }
    print $log_fh ''.( ($type eq PROTEIN) ? "Frame" : "Orientation")." data written in file ". 
	File::Spec->abs2rel($frames_file).".\n";
    print $log_fh "Reads ". ( ($type eq PROTEIN) ? "translated (possibly in more than one frame)" :
			      "correctly oriented")."\n  and written in file ". 
			      File::Spec->abs2rel($trans_file).".\n";
    
}

####### POSSIBLY SCAN READS IF USING AMPHORA

if (defined ($tr_reads_file) and (! defined ($trans_file)) ){
    $trans_file = $tr_reads_file;
}	    
if ( (!$noscan_flag) and (defined($gene_symbol)) and ($type eq PROTEIN) ) {
    unless (defined ($trans_file) and (-e $trans_file)) {
	die "Cannot find file containing reads to be scanned.\n";
    }
    print $log_fh "Scanning peptide reads with modified AMPHORA code markerScannerAlt ".
	"\n  for recognizable $gene_symbol fragments; this may take some time.\n"; 	  
    @args = ('-g', $gene_symbol, '-r', $trans_file, '-d', $out_dir, '-a', $amphora_path, '-h', $hmmer2_path);
    $scanned_pep_reads_file = run_markerScannerAlt(@args);
    if (!(-e $scanned_pep_reads_file)) {
	die "Error scanning reads!  Expecting file $scanned_pep_reads_file to have been written but cannot find it.\n";
    }
    my $count = count_seqs($scanned_pep_reads_file);
    print $log_fh "$count sequences picked up by modified AMPHORA code markerScannerAlt\n".
	"  as (possibly partial) $gene_symbol sequences.\n";
    print $log_fh "  Results of scan written to ". File::Spec->abs2rel($scanned_pep_reads_file).".\n";    
    $reads_to_align = $scanned_pep_reads_file;
} else {
    $reads_to_align = $trans_file;   
}

####### POSSIBLY FILTER READS, ACCORDING TO TARGET NUMBER DESIRED ALSO
####### CREATE A FILE OF FULL-LENGTH SEQUENCES CORRESPONDING TO THE
####### FINAL SET OF READS

if ( $num_seqs or defined($sample_file_dna) ) {
    if ($num_filtered_reads < 0) {
	print $log_fh "Filtering the reads so that there is at most one read per original sampled sequence.\n";
    } elsif ($num_filtered_reads > 0) {	
	print $log_fh "Filtering the reads:  choosing a subset of size $num_filtered_reads uniformly at random.\n";
    }
    unless (-e $reads_to_align) {
	die "Cannot find file $reads_to_align containing reads to be filtered.\n";
    }
    ($filtered_sim_reads_file, $sample_final_src_seqs, $final_num_reads, $final_num_src_seqs) =
	filter_reads($reads_to_align, $num_filtered_reads, $sample_file_dna, $out_dir, 
		     ($type eq PROTEIN) ? $sample_file_pep : $sample_file_dna ); 
    print $log_fh "Determining the full-length sequences that correspond to the final set of reads\n";
    print $log_fh "Each read id converted to include the id for the corresponding full-length sequence.\n";
    print $log_fh ''.$final_num_reads. " reads written to ". File::Spec->abs2rel($filtered_sim_reads_file).".\n";
    $reads_to_align = $filtered_sim_reads_file;
    if (-e $sample_final_src_seqs) {
	print $log_fh ''.$final_num_src_seqs.
	    " full-length sequences corresponding to these reads written to file ". File::Spec->abs2rel($sample_final_src_seqs)."\n";
    }
} 

####### ALIGN READS, POSSIBLY WITH REFERENCE DATABASE SEQUENCES

if ($align_flag) {
    unless (defined($reads_to_align) and (-e $reads_to_align) ){
	die "Cannot find file containing reads to be aligned.\n";
    }
    my $aligner;
    if ($type ne PROTEIN) {
	$aligner = CMALIGN;
    } else {
	if (!defined($alt_seqs_basename)) {
	    $aligner=AMPHORA;
	} else {
	    $aligner = HMMALIGN;
	}
    }
    if (($aligner ne AMPHORA) and !(-e $model)) {
	die "Please input existing model for alignment with option '--model'!\n";
    }
    if ( ($aligner eq AMPHORA) and (!defined ($gene_symbol)) ) {
	die "Please use option '-g' to provide a gene symbol (for use with AMPHORA) or option '--model'\n"
	    ." to provide a profile for alignment of sequences.\n";
    }
    
    print $log_fh "Aligning the reads in ".  File::Spec->abs2rel($reads_to_align);
    if ($use_refdb) {
	if ($type eq PROTEIN) { 
	    print $log_fh "\n  and reference sequences in $ref_db_file_pep ";
	} else {
	    print $log_fh "\n  and reference sequences in $ref_db_file_dna ";
	}
    } elsif (defined ($ref_db_file_basename)) {
	print $log_fh "\n  (without reference database sequences) ";
    }
    print $log_fh "\n  to ".
	( !defined($alt_seqs_basename) ? "AMPHORA profile HMM for $gene_symbol, ".
	  "using AMPHORA script MarkerAlignTrim.pl;\n"."  also, trimming according to the mask;\n" 
	  : "input model $model, using $aligner;\n" ).	  	      
	  "  this may take some time.\n";	
    
    my ($alignment, $align_log, %align_args);
    
    %align_args = ('-r', $reads_to_align, '-d', $out_dir, '-p', 1, '-t', $type, '-a', $aligner);
    if (defined($gene_symbol)) {
	$align_args{'-g'}= $gene_symbol;
	$align_args{'--amphora'} = $amphora_path;
    } elsif (defined($model)) {
	$align_args{'-m'}=$model;
    } else {
	die "Must provide either a model via option '--model' or a gene symbol via option '-g'\n".
	    "  (for use with AMPHORA) to do alignment!\n";
    }
    if ($use_refdb) {
	$align_args{'-b'}= ($type eq PROTEIN) ? $ref_db_file_pep : $ref_db_file_dna;
    }
    
    ($alignment, $align_log) = 
	run_alignTrim(%align_args);
    
    if (!(-e $alignment)) {
	die "Error in alignment process!  Expecting file $alignment to have been created, but cannot find it.\n";
    }
    if (-e $align_log) {
	if ($align_log_flag) {
	    print $log_fh "Information about or errors in alignment process written to log file ".
		"\n  ".File::Spec->abs2rel( $align_log ).".\n";
	} else {
	    unlink($align_log);
	}
    }
    print $log_fh "Reads and reference sequences have been aligned to profile;".
	(($aligner eq AMPHORA) ? ",\n  alignment has been trimmed;" : '').
	"\n  alignment written to ". File::Spec->abs2rel( $alignment).".\n";
    
    if ($align_src_flag and (defined $sample_final_src_seqs) and (-e $sample_final_src_seqs)) {
	print $log_fh "Aligning the full-length sequences in ".  
	    File::Spec->abs2rel($sample_final_src_seqs) .
	    "\n  that correspond to the final set of reads";
	if ($use_refdb) {
	    if ($type eq PROTEIN) { 
		print $log_fh "\n  and reference sequences in $ref_db_file_pep ";
	    } else {
		print $log_fh "\n  and reference sequences in $ref_db_file_dna ";
	    }
	}
	print $log_fh "\n  to ".
	    ( !defined($alt_seqs_basename) ?  "AMPHORA profile HMM for $gene_symbol, ".
	      "using AMPHORA script MarkerAlignTrim.pl;\n".
	      "  also, trimming the alignment according to the mask;\n" 
	      : "input model $model using $aligner;\n" ).	  	      
	      "  this may take some time.\n";	
	$align_args{'-r'}=$sample_final_src_seqs;
	$align_args{'-p'}=0;
	($alignment, $align_log) = run_alignTrim(%align_args);
	if (!(-e $alignment)) {
	    die "Error in alignment process!  Expecting file $alignment to have been created, but cannot find it.\n";
	}
	if (-e $align_log) {
	    if ($align_log_flag) {
		print $log_fh "Any errors or warnings in alignment process\n".
		    "  written to log file ".
		    File::Spec->abs2rel( $align_log ).".\n";
	    } else {
		unlink($align_log);
	    }
	}
	print $log_fh "Full-length sequences and reference sequences have been aligned".
	    (($aligner eq AMPHORA) ? ",\n  and alignment has been trimmed\n" : "\n").
	    "  and written to file ". File::Spec->abs2rel($alignment).".\n";
    }
    
}



