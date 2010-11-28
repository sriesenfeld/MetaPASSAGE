#!/usr/bin/perl
# S. Riesenfeld
# Last updated: Nov 2010
# Created: May 2009 (roughly)
#
# This script does two things:

# (1) It builds a reference database, which is a set of sequences that
# is meant to be used with MetaPASSAGE.pl to simulate a priori
# knowledge of full-length gene sequences from an input universe of
# sequences. (The sampling MetaPASSAGE.pl does to create a population
# can have a specified intersection size with the reference database,
# e.g., half or 0 or all of the population sequences can be chosen
# from the reference database.) The reference database can be formed
# by sampling uniformly at random from the input sequences or by
# maximizing (for part or all of the sequences in the reference
# database) the sequences' phylogenetic diversity.  This is explained
# further below.  [Coming evenutally: The Reference DB can be created
# by appending to an input smaller Ref DB].
#
# A phylogenetic tree with leaves whose labels contain taxonomic
# information matching that in the sequence headers of the input
# sequences can be used in two ways. (As an example of this type of
# tree, see the AMPHORA file 'Genome.tree' and a file of AMPHORA
# Reference Sequences, such as, rpoB.pep.) First, it can be used to
# control the diversity of the reference database, by choosing the
# sequences contributing the most incremental phylogenetic
# diversity. Second, it can be used to control the sequences that
# appear in the reference database. If a tree is provided, only
# sequences coming from taxa that appear in the tree can be selected
# for the reference database, even if phylogenetic diversity is not
# used as a parameter for the selection of reference database
# sequences.
#
# (2) A prerequisite for running the full simulation pipeline
# MetaPASSAGE.pl with protein sequences is to have a file of the
# corresponding DNA sequences. If you use any of the AMPHORA protein
# families, then this script can be used to retrieve the corresponding
# DNA sequences from NCBI.  For each peptide sequence from the set of
# AMPHORA Reference Sequences for a specified gene family, the
# corresponding DNA sequence is retrieved and written to file.  If you
# are not working with AMPHORA, you can modify this script and/or the
# functions in the retrieveSeqs.pm module to retrieve DNA sequences
# that you need, or you may need to acquire them in some other way.
#
# Known Bugs:
#
# Some of the gene ids used for the AMPHORA Reference Sequences are
# not known to NCBI -- for now, these sequences are just left out.
#
# Also, the retrieval of available sequences from the NCBI databases
# occasionally (seemingly randomly) goes on the fritz, probably due to
# problems with RefSeq and GenBank.  Most errors and warnings are
# trapped and dealt with, but be prepared.  If you don't get back a
# full set of gene sequences, just re-run the script one or two times
# (using the append option), and it should work.
 

use strict;
use warnings;

use File::Spec;

use Getopt::Std;

use lib 'simPipeModules';

use simPipeVars qw(:amphora :misc_file_ext :simpipe :maxpd);
use refSeqs qw(make_reference_db sample_seqs);
use retrieveSeqs qw( retrieve_seqs extract_gene_ids );

my $usage = qq{
Usage: $0 <options (see below)>

    -h [no value taken; 
	if set, this help message is printed and program ends]

    -g <gene symbol from AMPHORA; see options '-d' and '-f'>

    -d <the path to the AMPHORA directory; the default for this value
        is set in simPipeVars.pm>

    -r [no value taken; if set, retrieve DNA sequences corresponding
        to the gene symbol input with option '-g'; automatic retrieval
        is only currently implemented for AMPHORA sequences; also see
        options '-f', '-d', and '-a'; if '-r' is set, no
        reference-database-related options are processed]

    -a [no value taken; if set and option '-r' is set, the file
        specified by option '-f' (or the default file -- see the info
        on option '-f') is appended to rather than rewritten so that
        it contains each DNA sequence once for the AMPHORA taxa and
        the gene specified; if set and option '-r' is not set: THIS
        CASE IS NOT YET IMPLEMENTED: if option '-b' is set and
        reference db files of that basename exist already, then they
        are simply appended to create a reference db with a total
        number of sequences specified by option '-n'; this latter
        option is only implemented for the case when option '-m' is
        set to 0, so the appended sequences are chosen uniformly at
        random from the sequences that are not already in the ref db.]

    -n <total number of sequences to choose to form a reference
        database for a specific family of sequences (specified by
        option '-g' or by option '-f'); default value is set in
        simPipeVars.pm>

    -m <number of sequences to choose from the set using maxPD
        (according to a phylogenetic tree input using option '-t' or
        from maxPD_list() information input using option '-o'; using
        this option requires that the sequence headers have taxonomic
        information; this option works well with AMPHORA but may not
        work with other input sequences and trees; default value is
        set in simPipeVars.pm; the remaining sequences are sampled
        uniformly at random from the unselected sequences and written
        to a file whose base name can be chosen automatically or
        specified by option '-b'; if option '-m' is set to a positive
        integer, then either option '-o' must be set to a file
        containing output of maxPD_list() on a tree, or a tree (see
        option '-t') must be specified (or both)>

    -b <(optional) base name for the output reference database files
         (should include the path, and the directory structure of the
         path must already exist)>

    -t <(optional) file name for Newick tree to be given as input to
         maxPD_list(); default is the AMPHORA Reference Genome tree>

    -c <(optional) maximum number of taxa in tree for which maxPD
         should be computed; the taxa in the output are the taxa used
         to populate the reference database, both by maxPD value and
         by random sampling (can be set to number of taxa in input
         tree); this default is set to work well with AMPHORA protein
         families but should probably be set explicitly for other use>

    -o <(optional) file name for maxPD_list() output; if option '-t'
        is also set, then the output of running maxPD_list() on the
        given tree is written to the file with this name (overwriting
        existing file); if option '-t' is not set, then this file
        should contain the output of running maxPD_list() on a tree
        whose labels contain taxa represented in the sequences being
        used to create the reference database> 

    -f <(optional) base file name (including path, but without the
        extension) for fasta file(s) containing sequences for ref db
        creation. If '-g' <gene> option is set, then it is assumed
        that the peptide sequences for refdb simulation should be
        taken from AMPHORA/Reference/Sequence/<gene>.pep. It is also
        assumed that the corresponding DNA sequences are in the same
        directory with the name <gene>.fna. If '-f' <basename> is
        used, then it is assumed that <basename>.fna contains DNA
        sequences and <basename>.pep, if it exists, contains the
        corresponding peptide sequences. (The latter is optional.)
        Using option '-g' causes option '-f' to be ignored.>

};
 
# hash to contain options
my %opts;
if ( (!getopts('hg:ran:m:b:o:t:c:d:f:', \%opts)) or (! (%opts) ) or ($opts{'h'}) ) {
    die "$usage\n";
}
my ($gene_symbol, $total_num_seqs, $maxpd_num_seqs, $basename, $dir,
$maxpd_outfile, $treefile, $maxpd_compute_num, $all_seqs_basename);
$all_seqs_basename = $opts{'f'};
print "Command given: $0";
while( my ($option, $value) = each(%opts)) {
    print ' '.'-'.$option.' '.$value;
}
print "\n";

$gene_symbol = $opts{'g'};
if (!defined ($gene_symbol) and !defined ($all_seqs_basename)) {
    die "Either an AMPHORA symbol (option '-g') or a file basename (option '-f') must be provided.\n";
    
}
if ($opts{'r'} and !defined ($gene_symbol)) {
    die "Automatic sequence retrieval is only defined for AMPHORA gene families.".
	"  To use option '-r', option '-g' must also be used.\n";
}
# if (!$opts{'r'} and $opts{'a'} and ($opts{'m'}>0)) {
#    die "The reference DB files cannot be formed by appending to input ref DB files if option \'m\' is set to a positive number!\n";
# }
if (defined ($opts{'d'})) {
    define_amphora_paths(File::Spec->rel2abs($opts{'d'}));
}

if (!$opts{'r'}) {
    if (!$opts{'n'}) {
	$total_num_seqs = $ref_db_default_size;
    } else {
	$total_num_seqs = $opts{'n'};
    }
    $maxpd_num_seqs = $opts{'m'};    
    
    if (defined ($maxpd_num_seqs) and ($total_num_seqs < $maxpd_num_seqs)) {
	die "The total number of sequences in the reference database must be at least the number chosen by maxPD.\n";
    }
    $basename = $opts{'b'};    
    if (!defined($basename)) {
	$basename = $ref_db_default_filename_part.(defined($gene_symbol) ? '-'.$gene_symbol : '');
    } 
    $treefile = $opts{'t'};
    $maxpd_outfile = $opts{'o'};
    if ($maxpd_num_seqs and (!$maxpd_outfile or !(-e $maxpd_outfile)) and (!$treefile)) {
	$treefile = $amphora_genome_tree;
	print "No tree file specified; using as default the AMPHORA Genome tree in $amphora_genome_tree.\n";
	if (!(-e $amphora_genome_tree)) {
	    die "Cannot find AMPHORA Genome tree: $amphora_genome_tree!\n";
	}
    }
    
    if ($maxpd_num_seqs and !$maxpd_outfile) {
	my $tree_basename;
	my ($vol, $dirs, $file) = File::Spec->splitpath( $treefile);
	if ($file =~ /(.+)\./) {
	    $tree_basename = $1;
	}
	$maxpd_outfile = ($tree_basename ? $tree_basename.'-' : '') . $maxpd_filename;
    } 
    
    $maxpd_compute_num = $opts{'c'};
    if (!$maxpd_compute_num) {
	$maxpd_compute_num = $maxpd_default_num;
    }
}
$dir = $amphora_ref_seq_path;

my ($all_seqs_file_unaln, $all_seqs_file_dna);
if (defined($gene_symbol)) {
    $all_seqs_basename = File::Spec->catfile($dir, $gene_symbol);
}
$all_seqs_file_dna = $all_seqs_basename.$dna_fasta_ext;
$all_seqs_file_unaln = $all_seqs_basename.$pep_fasta_ext;

if (defined ($gene_symbol) and !(-d $dir)) {
    die "Cannot find directory $dir containing the AMPHORA Reference Sequences.\n";    
}
if (!$opts{'r'} and !(-e $all_seqs_file_dna)) {
    die "Cannot find DNA sequence file $all_seqs_file_dna!\n";
}
if ($opts{'r'} and (defined ($gene_symbol)) and !(-e $all_seqs_file_unaln)) {
    die "Cannot find peptide sequence file $all_seqs_file_unaln!\n";
} 

if ($opts{'r'} and (defined($gene_symbol)) and 
    ( defined($all_seqs_file_dna) and (!(-e $all_seqs_file_dna) or $opts{'a'}) ) ) {
    print "Retrieving DNA sequences... this may take a while...\n";
    my $num_retrieved;
    if ($opts{'a'}) {
	$num_retrieved=retrieve_seqs($all_seqs_file_unaln, $all_seqs_file_dna, 1);
    } else {
	$num_retrieved=retrieve_seqs($all_seqs_file_unaln, $all_seqs_file_dna );
    }
    print "$num_retrieved sequences retrieved from NCBI and written to $all_seqs_file_dna.\n";
} 

if (!$opts{'r'}) {
    if (defined ($gene_symbol)) {
	make_reference_db( '-n', $total_num_seqs, '-m', $maxpd_num_seqs, 
			   '-g', $gene_symbol, '-o', $basename, 
			   '-d', $dir, '-p', $maxpd_outfile, '-t', $treefile, 
			   '-c', $maxpd_default_num ); #, '-a', $opts{'a'}); 
    } else {
	make_reference_db( '-n', $total_num_seqs, '-m', $maxpd_num_seqs, 
			   '-o', $basename, 
			   '-p', $maxpd_outfile, '-t', $treefile, 
			   '-c', $maxpd_default_num, 
			   '-f', $all_seqs_basename );  #, '-a', $opts{'a'});
    } 
}





