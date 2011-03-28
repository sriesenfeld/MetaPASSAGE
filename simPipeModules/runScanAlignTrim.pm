#!/usr/bin/perl
# S.J. Riesenfeld 
# Updated: Nov 2010 
# MetaPASSAGE module: runScanAlignTrim
#
# Module for running AMPHORA scripts and modules adapted from AMPHORA,
# as well as other alignment programs, including INFERNAL and HMMER.

package runScanAlignTrim;

use strict;
use warnings;

use File::Spec;
use File::Copy;

require File::Temp;    # BUG:  newdir() not found right now!!
use File::Temp qw(tempdir);

require Exporter;
our @ISA = qw(Exporter);
our @EXPORT_OK = qw(run_markerScannerAlt run_alignTrim);

use simPipeVars qw(:misc_file_ext :amphora :simpipe :seqtypes);
use AMPHORA_MarkerScannerAlt qw(markerScannerAlt);

=head2 run_markerScannerAlt

  Usage : run_markerScannerAlt(%args);
  Function: Wrapper for running a module that implements some of the
    functionality of the AMPHORA script MarkerScanner.pl. Does an
    hmmsearch for the input file $readsfile of reads against the
    marker gene families in the AMPHORA installation. The reads that
    hit and belong to the gene family corresponding to the input
    $gene_symbol are returned.
  Example :
    # Assume $reads_file is a set of peptide reads, created by
    # simulation with MetaSim from the DNA sequences for the 'rpoB'
    # gene family and then translated in the correct frame.
    my $rpoB_read_hits = run_markerScannerAlt('-g', 'rpoB', '-r',
					      $reads_file, '-d',
					      $working_directory);
  Returns : a filename 
  Args : A hash %args where the keys of %args must include '-r' and
    '-g' and may optionally include '-d', '-f', '-h', and
    '-a'. $args{'-r'} is a filename of peptide reads. $args{'-g'} is a
    string used to name one of the gene families in
    AMPHORA. $args{'-d'} is the directory where output files should be
    written. $args{'-f'} is the desired name of the output
    file. $args{'-a'} is the path to the AMPHORA
    directory. $args{'-h'} is the path to HMMER 2 binaries.

=cut

sub run_markerScannerAlt (%) {
    
    my %options = @_;
    my ($gene, $outdir, $outfile, $readsfile);

    if (!($readsfile = $options{'-r'})) {
	die "You must specify a file containing reads for run_MarkerScannerAlt.\n";
    }
    if (!($gene = $options{'-g'})) {
	die "You must specify a gene for run_MarkerScannerAlt.\n";
    }
    # my $workingdir;
    my $curdir = File::Spec->curdir();
    $curdir = File::Spec->rel2abs($curdir);
    if (!($outdir = $options{'-d'})) {
	$outdir = $curdir; 
    }    
    $outdir = File::Spec->rel2abs( $outdir ) ;    
    if (defined ($options{'-a'})) {
	$amphora_path = $options{'-a'};
	define_amphora_paths($amphora_path);
    }
    if (defined ($options{'-h'})) {
	$hmmer2_path = $options{'-h'};
    }
    my $scanned_reads;
    if (!($outfile = $options{'-f'})) {
	my $basename = get_basename($readsfile, $pep_fasta_ext);
	$outfile = "$basename-$gene".$pep_fasta_ext;	
	$scanned_reads = File::Spec->catfile($outdir, $outfile);
    } else {
	$scanned_reads = File::Spec->rel2abs($outfile);
    }    
    my $tempdir = tempdir ( 'tempdirXXXXX', DIR => $outdir, CLEANUP=>1 );
    my $stuff_to_scan = 'query'.$pep_fasta_ext;
    $stuff_to_scan = File::Spec->catfile ($tempdir, $stuff_to_scan);
    copy($readsfile, $stuff_to_scan) or 
	die "Cannot copy $readsfile to temporary directory $tempdir to run (alternate) AMPHORA scan: $!\n";
    chdir ($tempdir) or die "Cannot move to temporary working directory $tempdir to run AMPHORA scan: $!\n";
    markerScannerAlt($stuff_to_scan, $amphora_path, $hmmer2_path);
    
    my $temp_outfile = $gene.$amphora_ref_file_ext;
    
    if (! (-e $temp_outfile)) {
	die "Expecting markerScannerAlt to have created file $temp_outfile in temporary directory ".$tempdir.
	    ", but cannot find it.\n";
    }
    move($temp_outfile, $scanned_reads) or die "Cannot move temporary file $temp_outfile to $scanned_reads: $!\n";
    chdir ($curdir) or die "Cannot move back to original working directory $curdir: $!\n";
    if (-d $tempdir) {
	system("rm -r $tempdir");
    }
    
    return $scanned_reads;
}

=head2 run_alignTrim

  Usage : run_alignTrim(%args);
  Function: Wrapper to build alignments of the input reads or
    sequences, which are assumed to belong to a single gene family,
    along with reference database full-length sequences (if provided),
    using AMPHORA, HMMER without AMPHORA, or INFERNAL. If either
    non-AMPHORA method is used, a corresponding profile model for the
    gene family must be provided to align the sequences.
  Example: 
    # First example: Assume $reads_file is a set of peptide reads
    # identified as belonging to the 'rpoB' gene family in AMPHORA and
    # let $refdb be a file containing some full-length rpoB sequences
    my ($alignment, $logfile) = run_alignTrim('-r', $reads_file, '-g',
					      'rpoB', '-b', $refdb,
					      '-p', 1);
    # Second example: Assume $reads_file is a set of 16S rRNA reads
    # (previously converted to DNA) that have been correctly oriented,
    # and let $model_file be a file containing an INFERNAL 16S rRNA
    # CFG model.
    my ($alignment, $logfile) = run_alignTrim('-r', $reads_file, '-p',
					      1, '-a', 'cmalign', 
					      '-m', $model_file, 
					      '-t', 'RNA2DNA'); 					      
  Returns : an array of two filenames; the first contains the alignment
    and the second contains a log of output from the alignment processes
  Args : A hash %args where the keys of %args must include '-r' and
    may include '-g', '-a', '-p', '-m', and '-t'; other completely
    optional keys are '-b', '-d', '-l', '-f', and '--amphora'.
    $args{'-r'} is a file of sequences, either reads or full-length,
    from the same gene family.  If the gene family is from AMPHORA,
    then it should be specified with $arg{'-g'}. Otherwise,
    $args{'-a'} should specify the method of alignment to be used,
    either 'cmalign' from INFERNAL or 'hmmalign' from HMMER. (This
    program must be in the user path.)  If the sequences to be aligned
    are not full-length, then $arg{'-p'} should be set to 1. If any
    non-AMPHORA method of alignment is chosen, then $args{'-m'} should
    be a file containing a profile model for the gene family
    corresponding to the alignment method. If non-peptide sequences
    are being aligned $args{'-t'} must be set to 'RNA2DNA', which
    means the sequences are RNA sequences that have been converted to
    DNA (the only non-amino-acid option implemented
    currently). $args{'-b'} is a file containing full-length
    'reference database' sequences from the gene family.  $args{'-d'}
    is the directory where where output files should be
    written. $args{'-l'} specifies the desired log file
    name. $args{'-f'} is the desired name of the output
    file. $args{'--amphora'} is the path to AMPHORA.

## NOTE!!: ADD COMMENTS HERE ABOUT ALIGNING --WITHALI USING ALN ON
## REFDB, also may want to add the facility for hmmalign

=cut

sub run_alignTrim (%) {    
    my %options = @_;
    my ($gene, $outdir, $logfile, $outfile, $seqfile, $refdb, 
	$partial_flag, $aligner, $type, $model, $withali_aln);
    
    if (!($seqfile = $options{'-r'})) {
	die "You must specify a file containing sequences for run_alignTrim.\n";
    }
    $refdb = $options{'-b'};
    $outfile = $options{'-f'};
    $partial_flag = $options{'-p'};   
    $aligner = $options{'-a'};  # 0 is AMPHORA, 1 is cmalign (for now)
    $type = $options{'-t'}; 
    $model = $options{'-m'};
    $withali_aln = $options{'-w'};
    
    if (defined( $options{'--amphora'})) {
	$amphora_path = $options{'--amphora'};
	define_amphora_paths($amphora_path);
    }
    
    if (!$type) { $type = PROTEIN };
    if (!$aligner) { $aligner = AMPHORA };
    if (!scalar(grep($aligner, @aligners))) {
	die "Must choose a method of alignment from the following: ".join(@aligners, ', ')."\n";
    }
    if (($aligner eq AMPHORA) and !($gene = $options{'-g'})) {
	die "You must specify a gene symbol for run_alignTrim to use ". AMPHORA ." for the alignment.\n";
    }
    if ( ($aligner ne AMPHORA) and !$model) {
	die "Please provide a model to align the sequences to!\n";
    }
    my $curdir = File::Spec->curdir();
    $curdir = File::Spec->rel2abs($curdir);
    if (!($outdir = $options{'-d'})) {
	$outdir = $curdir; 
    }    
    $outdir = File::Spec->rel2abs( $outdir ) ;
    my ($basename, $alignment);
    if (!$outfile) {
	$basename = get_basename($seqfile);
	$outfile = ''.$basename. ( defined($refdb) ? "-$ref_db_default_filename_part" : '').
	    (($aligner eq AMPHORA) ? $aln_fasta_ext : $stockholm_aln_ext);	
	$alignment = File::Spec->catfile($outdir, $outfile);
    } else {
	$alignment = File::Spec->rel2abs($outfile);
	$basename = get_basename($outfile, ($aln_fasta_ext, $stockholm_aln_ext));
    }
    if (!($logfile = $options{'-l'})) {
	$logfile = File::Spec->catfile($outdir, $basename.$align_log_filename);
    }
    $logfile = File::Spec->rel2abs($logfile);
    
    my $temp_stuff_to_align;
    my ($temp_fh, $tempdir);
    
    if ($aligner eq AMPHORA) {
	$tempdir = tempdir ( 'tempdirXXXXX', DIR => $outdir, CLEANUP=>1 );    
	$tempdir = File::Spec->rel2abs($tempdir);
	$temp_stuff_to_align = File::Spec->catfile($tempdir, $gene.$amphora_ref_file_ext);	
	
    } else {
	my $temp_fh = File::Temp->new(TEMPLATE => 'tempXXXXX');
	$temp_stuff_to_align = $temp_fh->filename();
    }
    
    if (defined($refdb)) {
	copy ($refdb, $temp_stuff_to_align) or die "Cannot copy $refdb to temporary file $temp_stuff_to_align to run AMPHORA scripts: $!\n";
	open(ALIGN, ">>$temp_stuff_to_align") 
	    or die "Cannot open temporary file $temp_stuff_to_align of sequences: $!\n";
	open(SEQ, "$seqfile") or die "Cannot open file $seqfile to read sequences: $!\n";
	while (my $line = <SEQ>) {
	    print ALIGN $line;
	}
	close(ALIGN);
	close(SEQ);
    } else {
	copy ($seqfile, $temp_stuff_to_align) or die 
	    "Cannot copy $seqfile to temporary file $temp_stuff_to_align to perform alignment: $!\n";
    }    
    # want to align reference sequences and sequences together
    # will cheat a bit and use the profile made from all sequences in our real db
    
    my $cmd;
    if ($aligner eq AMPHORA) {
	chdir ($tempdir) or die "Cannot move to temporary working directory $tempdir to perform alignment: $!\n";
	my $temp_logfile = 'temp.log';
	my $amphora_align_script = File::Spec->catfile($amphora_scripts_path, 'MarkerAlignTrim.pl');
	if (! (-e $amphora_align_script)) {
	    die "Cannot find AMPHORA script MarkerAlignTrim.pl in directory $amphora_scripts_path!\n";
	}
	$cmd='perl '.$amphora_align_script.' -Trim '. ($partial_flag ? '-Partial ' : ''). '2> '. $temp_logfile;
	eval {
	    system ($cmd);
	};
	if ($@) {
	    die "Error running AMPHORA's MarkerAlignTrim.pl: $@\n";
	}
	my $temp_outfile = $gene.$amphora_ref_aln_file_ext;
	
	if (! (-e $temp_outfile)) {
	    die "Expecting AMPHORA script MarkerAlignTrim.pl to have created file $temp_outfile in temporary directory ".$tempdir.
		", but cannot find it.\n";
	}
	move($temp_outfile, $alignment) or die "Cannot move temporary file $temp_outfile to $alignment: $!\n";
	if (-e $temp_logfile) {
	    move($temp_logfile, $logfile) or die "Cannot move temporary log file $temp_logfile to $logfile: $!\n";
	}
	chdir ($curdir) or die "Cannot move back to original working directory $curdir: $!\n";

    } else {
	
 	if ($aligner eq CMALIGN) {	 
	# Would like to first align refdb to model, then fix that alignment, but cannot do this.
	# To do something like this, we would have to actually excise RefDB sequences from 
        # the alignment used to build the model and maintain the RF line (since --rf cannot be used without --withali).
	# To be parallel, should we just align all sequences to model simultaneously, as in AMPHORA? 
	# Doing this for now.
	
	# my $aligned_refdb = File::Spec->catfile($tempdir, $refdb.$stockholm_aln_ext);	
	# my $cmd = "cmalign --rf --hbanded -o $aligned_refdb".
	  #  (($type eq RNA2DNA) ? '--dna':''). " $model $refdb 2> $logfile";  ## DOESN'T WORK

	    $cmd = 'cmalign --hbanded'.($partial_flag ? ' --sub' : '').
		(($type eq RNA2DNA) ? ' --dna':''). 
		' -o '.$alignment.
		(defined ($withali_aln) ? " --withali $withali_aln --rf " : ' ').
		$model.' '. $temp_stuff_to_align
		.' 2>&1 > ' .$logfile; 
	} else {
	    # could try to do this with --withali option and align refdb first
	    $cmd = 'hmmalign -o '.$alignment.' '. $model.' '. $temp_stuff_to_align.
		' 2>&1 > '.$logfile;
	}
	eval {
	    system ($cmd);
	};

	if ($@) {
	    die "Error running ".$aligner." on $seqfile: $@\n";
	}
	if (!(-e $alignment)) {
	    die "Expecting $aligner to have created file $alignment but cannot find it: $!\n";
	}
    }

    if (defined($tempdir) and (-d $tempdir)) {
	system("rm -r $tempdir");
    }    
    if (defined($temp_stuff_to_align) and (-f $temp_stuff_to_align)) {
	system("rm $temp_stuff_to_align");
    }
    
    return ($alignment, $logfile);   
}

1;

