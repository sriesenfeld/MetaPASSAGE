#!/usr/bin/perl
# S.J. Riesenfeld
# Last updated: Nov 2010
# MetaPASSAGE module: runMetaSim
#
# Module for running MetaSim with various default options, for
# creating MetaSim taxonomic profile, and for filtering the simulated
# reads in order to get one read per sequence.
#
# Note: MetaSim documentation says that the database is built by
# default in the installation directory, but in my experience, it ends
# up by default in the directory where MetaSim is called.  However,
# the log file does end up by default in the installation directory.
# This program tries to work around these defaults (which are not
# adjustable from the command line).  There are several other quirks
# of MetaSim (especially the command-line version) which this wrapper
# attempts to handle.

package runMetaSim;

use strict;
use File::Spec;
use File::Copy "mv";
use List::Util "min";

use File::Temp qw(tempdir);
use POSIX qw(ceil floor INT_MAX);

require Exporter;
our @ISA = qw(Exporter);
our @EXPORT_OK = qw(run_metasim filter_reads create_profile);

use refSeqs qw(get_sequence_headers get_random_subset get_gene_ids_from_headers);
use simPipeVars qw(:metasim :misc_file_ext :simpipe);

=head2 run_meta_sim

 Usage: run_metasim( <%args, a hash of options: require at least '-a'
   or '-m' option, as described below> ) Function: Wrapper function
   for running MetaSim from the command line. By default, uses a
   Sanger-based sequencing model, but without any sequencing
   error. (Sanger is used because it is easy to control from the
   command line and offers a default of no sequencing error.) There
   are options for modifying this model or specifying a separate
   empirical model. Defaults set in simPipeVars.pm.
 Notes: Default values of optional arguments set in simPipeVars.pm.
   The metasim database directory is assumed to be in the current
   directory. This wrapper guesses at the name of the file that
   MetaSim creates, which cannot be chosen by command-line input to
   MetaSim, and it will die if it cannot find the file. Regarding read
   length, it may make sense in simulations to control the read length
   via the clone length, i.e., to set the read length to be very large
   and the clone length to be what you want the read length to
   be. This seems to work best for the Sanger model. The default
   values for the clone and read lengths give a mean read length of
   about 400bp. Since command-line MetaSim assumes that the sequences
   are circular, it is useful to pad non-circular sequences before
   adding them to the database.  
 Example:
   # Delete any existing MetaSim database in this directory and add
   # sequences in file $file to a new database. 
   run_metasim('-a', $file, '-f', $log, '-x', 1);
   # The file $profile may be created with the create_profile
   # function.  The command below generates reads without error with a
   # mean read length of 400bp and then deletes the existing MetaSim
   # database.
   my $reads_file = run_matsim('-m',  $profile, '-r', 500, '-f', $log, 
    '-x', 1 );
   # To use a sequencing error configuration file $err_conf_file (like
   # the Illumina 80bp configure file that can be downloaded from the
   # MetaSim website:
   my $reads_file = run_matsim('-m',  $profile, '-r', 500, '-f', $log, 
    '-x', 1 , '-e', $err_conf_file);
 Returns: The name of the file containing reads
 Args : A hash %args where the key options and their values are
   described below. Either $args{'-a'} or $args{'-m'} must be defined;
   the rest are optional.
   Key: '-a'; Value: Fasta file of sequences to be added to the
     MetaSim database.
   Key: '-m'; Value: MetaSim taxonomic profile file, with extension
     '.mprf', giving the distribution of relative abundances of
     sequences to be used to generate reads.
   Key '-r'; Value: Number of reads; default set in simPipeVars.pm.
   Key '-f': Value: Name of log file to which will be appended the
     output MetaSim returns (not the reads) and the contents of the
     metasim error.log. (The error.log can only be reported if the
     metasim path is set correctly, which can be done with an option
     to the pipeline or in simPipeVars.pm.)
   Key '-d'; Value: Output directory; default is current directory.
   Key '-o'; Value: Output file name for the reads; the default is to
     rename the file MetaSim spits out, using a name based on the
     input file.
   Key '-l'; Value: mean read length, parameter passed to MetaSim;
     default set in simPipeVars.pm
   Key: '-s'; Value: standard deviation of read length distribution,
     parameter passed to MetaSim; default set in simPipeVars.pm
   Key: '-c'; Value: mean value of clone 'fragment length'; set by
     default to 2*(mean read length)
   Key: '-p'; Value: 'std/max deviation of fragment length distrib';
     set by default to ceil(.1*(clones mean))
   Key: '-x'; Value: 1, or 0/undefined: if set, remove existing
     database of sequences. If option '-a' is set, then any
     pre-existing database is deleted, before sequences are added to a
     new one (which is not deleted); if option '-m' is set, then the
     existing database is used to generate reads and then deleted. It
     is useful to delete the database for reasons beyond storage
     space; MetaSim complains if sequences that are the same but have
     different identifiers are stored, which can happen with gene
     sequences. Deleting the database between rounds helps avoid this
     situation.
   Key: '-e'; Value: set to the name of an input error model file,
     which is then used instead of the Sanger (without error) model
     used by default. 
   Key: '--mate_prob'; Value: the probability of paired reads for the
        MetaSim Sanger error model; default is 0
   Key: '--err_insert'; Value: the relative insertion rate for the
        MetaSim Sanger error model; default is 0
   Key: '--err_delete'; Value: the relative deletion rate for the
        MetaSim Sanger error model; default is 0
   Key: '-err_start'; Value: the initial error rate for the MetaSim
     Sanger error model; default is 0
   Key: '-err_end'; Value: the final error rate for the MetaSim Sanger
        error model; default is 0

=cut

sub run_metasim(%) {
    # Name of MetaSim output file is chosen by MetaSim based on input file;
    # name of error log is also chosen by MetaSim.    

    my ($errorLog, $cmd);
    if (-e $metasim_path) {
	$errorLog = File::Spec->catfile($metasim_path, $metasim_error_log_name); 
    }
    
    # my $cur_dir = File::Spec->rel2abs(File::Spec->curdir());
    # print "Running metasim from directory $cur_dir.\n";

    my %options = @_;
    my ($returnVal, $fastafile, $profilefile);
    
    my ($workingdir, $outdir, $readsfile);
    my ($mean_read_len, $stddev_read_len, $num_reads, $clones_mean, $clones_param2, 
	$mate_prob, $err_start, $err_end, $err_insert, $err_delete);
    if (! ($mean_read_len = $options{'-l'}) ) {
	$mean_read_len = $mean_read_len_default;
    }
    if (! ($stddev_read_len = $options{'-s'}) ) {
	$stddev_read_len = $stddev_read_len_default;
    }
    if (!($num_reads = $options{'-r'})) { 
	$num_reads = $num_reads_default;
    }
    unless ($num_reads > 0) { die "Require a positive number of reads."; }
            
    if (! ($clones_mean = $options{'-c'})) {
	$clones_mean = $mean_clone_len_default;
    }
    if (! ($clones_param2 = $options{'-p'})) {
	$clones_param2 = $stddev_clone_len_default;
    }
    if (! ($mate_prob = $options{'--mate_prob'})) {
	$mate_prob = $sanger_mate_prob_default;
    }
    if (! ($err_start = $options{'--err_start'})) {
	$err_start = $sanger_err_start_default;
    }
    if (! ($err_end = $options{'--err_end'})) {
	$err_end = $sanger_err_end_default;
    }
    if (! ($err_insert = $options{'--err_insert'})) {
	$err_insert = $sanger_insert_default;
    }
    if (! ($err_delete = $options{'--err_delete'})) {
	$err_delete = $sanger_delete_default;
    }
    if (! ($outdir = $options{'-d'})) {
	$outdir = File::Spec->curdir();
    }
    my $err_conf_file = $options{'-e'};
    if (defined($err_conf_file) and ! (-e $err_conf_file)) {
	die "Cannot find the error config model input $err_conf_file.\n";
    }
    my ($out_file, $out_fh, $tempdir_flag);
    unless ($options{'-a'} or $options{'-m'} or $options{'-x'}) { 
	die "Missing option values in run_metasim.\n";
    }
    if ($options{'-f'}) {
	$out_file = $options{'-f'};
	open(OUT_MS, ">>$out_file") or die "Cannot open log $out_file: $!\n";
	$out_fh = \*OUT_MS;
    } else { $out_fh = \*STDOUT; }
    
    if ($options{'-x'} and (-d $metasim_db_name_default) and (defined($options{'-a'}) || ! defined($options{'-m'}))) {
	# print "(Removing previously existing MetaSim database.)\n";
	eval {
	    system("rm -r database");
	};
	if ($@) {
	    die "Error removing old MetaSim database.\n";
	}
    }	
    if ($fastafile = $options{'-a'}) { 
	unless (-e $fastafile) { die "Cannot find file $fastafile of sequences to add to MetaSim database.\n";}
	# remove old metasim database so there are no conflicts in the names of matching gene sequences
        # print "Adding sequences from $fastafile to MetaSim database.\n";	
	$cmd = "MetaSim cmd --add-files $fastafile";
	print $out_fh "\nCommand given: \n$cmd\n";
	eval {
	    $returnVal = qx($cmd);
	    print $out_fh "\nMetaSim returns: \n$returnVal\n";
	    if (defined ($errorLog)) {
		print $out_fh "\nMetaSim error log: \n";
		printErrorLog($errorLog, $out_fh);
	    }
	};
	if ($@) {
	    die "Error adding sequences in MetaSim: $@\n";
	}
    }
    my $workingdir;
    if ($profilefile = $options{'-m'}) {
	unless (-e $profilefile) { die "Cannot find MetaSim profile file $profilefile.\n"; }
	if ($options{'-o'}) {
	    $readsfile = File::Spec->catfile($outdir, $options{'-o'});
	    $workingdir= tempdir ( 'tempdirXXXXX', DIR => $outdir, CLEANUP => 1);
	    $tempdir_flag =1;
	}
	if (! (defined $workingdir) or !(-d $workingdir)) {
	    $workingdir = $outdir;
	}	
	my ($volume, $directories, $file) = File::Spec->splitpath( $profilefile );
	my $basename = get_basename($file, $metasim_profile_ext);
	# The name of the metasim output file is chosen by metasim based on the input file;
	# below is a guess at the name it has chosen.
	my $metasim_reads_file = File::Spec->catfile($workingdir, $basename.
						     (defined($err_conf_file)? $metasim_error_conf_name :
						      $metasim_sanger_name).$metasim_fasta_ext); 
	if (-e $metasim_reads_file) {
	    unlink ($metasim_reads_file);
	}
	# Added the clones parameter settings; it seems to work well
	# to set the read length to be much larger than the clone
	# length and let the clone length actually determine the read
	# length; Added "uniform-weights" as a parameter to try to
	# keep MetaSim from using info about the genome size or copy
	# number.  Note that it seems that the default setting is to
	# assume sequences are ***circular***.  Not much I seem to be
	# able to do anything about this from the command line, so
	# adding padding and generating extra reads is recommended.
	if (!defined($err_conf_file)) {
	    $cmd = "MetaSim cmd --reads $num_reads --clones-mean $clones_mean --clones-param2 $clones_param2 -d $workingdir --uniform-weights --sanger --sanger-mean $mean_read_len --sanger-param2 $stddev_read_len --sanger-err-start $err_start --sanger-err-end $err_end --sanger-deletions $err_delete --sanger-insertions $err_insert --sanger-mate-probability $mate_prob $profilefile";
	} else {
	    $cmd = "MetaSim cmd -mg $err_conf_file --reads $num_reads -d $workingdir --uniform-weights $profilefile";
	}
	print $out_fh "\nCommand given: \n$cmd\n";
	eval {
	    $returnVal = qx($cmd);
	    print $out_fh "\nMetaSim returns: \n$returnVal\n";
	    if (defined ($errorLog)) {	    
		print $out_fh "\nMetaSim error log: \n";
		printErrorLog($errorLog, $out_fh);
	    }
	};
	if ($@) {
	    die "Error simulating reads in MetaSim: $@\n";
	}
	if (! (-e $metasim_reads_file)) {
	    die "Error:  Expecting file $metasim_reads_file to be written by MetaSim, but it has not been written.\n";	
	}
	if (defined $readsfile) {
	    mv($metasim_reads_file, $readsfile);	    
	} else {
	    $readsfile = $metasim_reads_file;
	}
	if ($options{'-x'} and (-d $metasim_db_name_default)) {
	    # print "Removing existing MetaSim database.\n";
	    eval {
		system("rm -r database");
	    };
	    if ($@) {
		die "Error removing old MetaSim database.\n";
	    }
	}
	# CLEANUP doesn't work consistently, so I'm deleting this directory explicitly.
	if ($tempdir_flag and (-d $workingdir)) {
	    system("rm -r $workingdir");
	}
	return $readsfile;
    } 
}

=head2 create_profile

  Usage : create_profile(%args);

  Function: Creates a file containing relative abundances for the
    input sequences (or sequence headers), in a format accepted by
    MetaSim. (MetaSim calls it a taxonomic profile and uses the
    extension '.mprf'.) By default, the abundances are all 1
    (signifying a uniform distribution for MetaSim sampling and
    simulation). An input parameter can specify a ratio to create a
    skewed distribution that more similarly mimics real microbial
    populations. By default, the skewed distribution is composed of a
    geometric distribution over half of the sequences and the
    remaining have minimum abundance value 1 (elongating the tail of
    the distribution). See the arguments for more information.

  Example :
    # Assume $file_in a fasta file containing sequences or just
    # sequence ids.
    my $file_out = create_profile('-i', $file_in); 
    # $file_out contains a uniform MetaSim taxonomic profile
    my $file_out= create_profile('-i', $file_in, '-k', 1.5); 
    # $file_out contains a skewed MetaSim taxonomic profile

  Returns : a filename 

  Args : A hash %args where the keys of %args must include '-i' and
    may optionally include '-k', '-a', '-n', and '-o'.  $args{'-k'} is
    the ratio used in the geometric distribution (it is the ratio
    between the second smallest relative abundance value and the
    minimum relative abundance value, which is 1). $args{'-a'}
    specifies the maximum abundance value permitted; by default, it is
    the maximum integer value. $args{'-n'} specifies the maximum
    number of sequences that should have relative abundance greater
    than 1; by default, this is half the number of sequences in the
    input file. $args{'-o'} specifies either an output directory where
    the output file should be written or, if that directory does not
    exist, the name of the output file.

=cut

sub create_profile (@) {
    my $outdir;
    my %options = @_;
    my $infile = $options{'-i'};
    if (!defined($infile)) {
	die "Provide an input file to create_profile using option '-i'.\n";
    }
    my $outfile = $options{'-o'};
    my $ratio = $options{'-k'};
    my $max_abundance = $options{'-a'};
    my $max_num_seqs_abund_gt1 = $options{'-n'};
    if (!defined($ratio)) {
	$ratio=1;
    }    
    if ($ratio < 1) {
	die "Input parameter for option '-k' must be at least 1 for create_profile!\n";
    }
    if (! defined ($max_abundance)) {
	$max_abundance = INT_MAX;
    }
    my $n_seqs;
    if ($ratio > 1) {
	$n_seqs = count_seqs($infile);
    }
    if ($ratio>1) {
	$max_num_seqs_abund_gt1 = floor($n_seqs/2);
    }
    my $n_seqs_abund_gt1=0;
    if ($ratio > 1) {
	$n_seqs_abund_gt1 = min($n_seqs-1, $max_num_seqs_abund_gt1);
    }
    my @abund_vals=();
    my $val = $ratio;
    my $max_abund_div_k=0;
    if ($ratio>1) {
	$max_abund_div_k= floor($max_abundance / $ratio);
    }
    if ($ratio > 1) {
	foreach my $i (1..$n_seqs_abund_gt1) {
	    push(@abund_vals, $val);
	    if ($val > $max_abund_div_k) {
		last;
	    } else {
		$val *= $ratio;
	    }
	}
    }
    open(IN, $infile) or die "Cannot open input file $infile. $!\n";
    if (-d $outfile) { $outdir = $outfile; }
    if ( (!$outfile) or ($outdir) ) {
	$outfile = get_basename($infile).$metasim_profile_ext;
    }
    if ($outdir) {
	$outfile = File::Spec->catfile($outdir, $outfile);
    }
    open(OUT_PROF, ">$outfile") or die "Cannot open output file $outfile. $!\n";
    
    print OUT_PROF "\# MetaSim taxon profile created from $infile;\n";
    print OUT_PROF "\# ratio of $ratio.\n";
    
    my $val = 1;
    if (scalar(@abund_vals)>0) {
	$val = pop(@abund_vals);
    }    
    my $count=1;
    while(my $line = <IN>) {
	if ($line =~ '>(.+)') {
	    print OUT_PROF ''.ceil($val)." name \"$1\"\n"; 
	    if (scalar(@abund_vals) > 0) {
		$val = pop(@abund_vals);
	    } else {
		$val = 1;
	    }
	}
    }
    close(IN);
    close(OUT_PROF);
    return $outfile;
}

=head2 filter_reads

 Usage: filter_reads (@args)
 Function: Writes a file containing a subset of the reads, with the
   ids are changed to have the form <read id given by MetaSim>-<gene
   id of source sequence>.  The subset of reads is chosen according to
   the value of N (see arg 2).  If N is undef or negative, then for
   each source sequence, one read is chosen from the set (if it is
   nonempty) of reads that correspond to that source.  If N is
   positive, then a subset of size N is chosen uniformly at random
   from all reads.  If N == 0, then no sub-sampling is done but ids
   are converted. If a file containing the corresponding full-length
   source sequences (i.e., the sequences given to MetaSim to produce
   the reads) is also given as input, then the subset of these
   sequences corresponding to the output set of reads is written to a
   separate output file.  Can be used more generally on files of
   peptide reads that are translations of MetaSim DNA reads, as long
   as the sequence headers have not changed.
 Example: 
   ($filtered_reads_file, $final_src_seqs_pep, $final_num_reads,
        $final_num_src_seqs) = 
     filter_reads('reads.fna', -1, 'src_seqs.fna', 'outdir',
       'src_seqs.pep');
 Returns: A list containing the name of the output reads file, the
   name of the output source sequence file, the number of reads
   selected and written to file, and the number of source sequences
   written to file.
 Args: A list of up to 6 arguments; first 2 are non-optional. 
  arg 1: file of reads generated by MetaSim (or their translations into
    peptides);
  arg 2: an integer N: N<0 if output should contain at most one
    read per sequence in file named by arg 3; if N==0, then all reads
    are maintained (can be useful with the optional args to get a file
    of source sequences corresponding to the reads); otherwise N is
    the number of reads that should be sampled uniformly at random
    from the file named by arg 1;
 arg 3: (optional, required if arg 2 is < 0) a fasta file
   containing headers (with or without sequences) for sequences used
   by MetaSim to generate the reads;
 arg 4: (optional) an output file name for the filtered reads or a
   directory in which all output files should be written;
 arg 5: (optional) a file containing sequences (DNA or peptides) that
   correspond to the sequences given to MetaSim as sources for the
   reads; if arg 5 is specified, then the full sequences corresponding
   to those reads which have been filtered are written to a separate
   output file;
 arg 6: (optional) the name of the output file for use with arg 5.

=cut

sub filter_reads ($$;$$$$) {
    my $outdir;
     my ($fragfile, $N, $seqfile_hdrs, $outfile, $seqfile, $seq_outfile) = @_;
    my @frag_hdrs = get_sequence_headers($fragfile);
    my (@chosen_frag_hdrs, $chosen_frag_metasim_ids_ar, @src_hdrs, $src_gene_ids_ar, @chosen_src_gene_ids);

    if (! defined($N)) {
	return;
    }
    my @exts=(@fasta_alt_exts, $pep_fasta_ext);
    my $pattern = qr{$pep_fasta_ext$};    
    my $ext;
    if (-d $outfile) { $outdir = $outfile; }
    if ( (!$outfile) or ($outdir)) {
	if ($fragfile =~ /$pattern/) {	    
	    $ext = $pep_fasta_ext;
	} else {
	    $ext = $metasim_fasta_ext;
	}	
	if ($N > 0) {
	    $outfile = get_basename($fragfile, @exts).'-'.min(scalar(@frag_hdrs), $N).'r'.$ext;
	} elsif ($N < 0) {
	    $outfile = get_basename($fragfile, @exts).$single_reads_name.$ext;
	} else {
	    $outfile = get_basename($fragfile, @exts).'-'.scalar(@frag_hdrs).'r'.$ext;
	}	    
    }
    if ($outdir) {
	$outfile = File::Spec->catfile($outdir, $outfile);
    }
    my $num_found=0;
    if ($N>0) {
	# pick random subset of size $N of the reads
	@chosen_frag_hdrs = get_random_subset($N, \@frag_hdrs);
    } elsif ($N < 0) {
	# pick one read uniformly at random from each set of reads that come from the same source sequence
	@src_hdrs = get_sequence_headers($seqfile_hdrs);
	$src_gene_ids_ar = get_gene_ids_from_headers(\@src_hdrs);
	foreach my $gene_id (@{$src_gene_ids_ar}) {
	    # print "Considering gene id: $gene_id\n";
	    my @reads_matched = grep { $_ =~ $gene_id } @frag_hdrs;
	    # print ''.scalar(@reads_matched)." read headers matched to that gene id:\n";
	    # print map {"  $_\n"} @reads_matched;
	    if (scalar(@reads_matched)) {
		my @chosen_match = get_random_subset(1, \@reads_matched);
		# print "Chosen match: ".$chosen_match[0]."\n";
		push(@chosen_frag_hdrs, @chosen_match);
		# $num_found++;
	    }
	}
    } else {
	@chosen_frag_hdrs = @frag_hdrs;
    }
    
    $chosen_frag_metasim_ids_ar = cut_to_ids(\@chosen_frag_hdrs);
    open (INFRAG, $fragfile) or die "Cannot open $fragfile containing reads generated by MetaSim: $!\n";
    open(OUT, ">$outfile") or die "Cannot open $outfile for writing single reads: $!\n";
    
    my $basic_header = qr{^>};
    my $hdr_remainder = qr{^>\S+\s(.*)};
    my $metasim_header = qr{$metasim_hdr_pattern};  # the pattern that can identify the gene id of the source sequence for a read
    my $copy=0;
    my $num_frags_written=0;
    while(my $linefrag = <INFRAG>) {
	# only do something if we find a header, in which case, extract the read id and the gene id
	if ($copy) {
	    if ($linefrag =~ $basic_header) {
		$copy=0;
	    } else {
		print OUT $linefrag;
		next;
	    }
	}
	
	my ($frag_id, $gene_id);
	if ($linefrag =~ /$metasim_header/) {
	    $frag_id = $1; 
	    $gene_id = $2;	    
	} else {
	    next;
	}
	if (! scalar( grep {$_ =~ $frag_id} @{$chosen_frag_metasim_ids_ar}) ) {
	    next;
	}
	if (! scalar( grep {$_ =~ $gene_id} @chosen_src_gene_ids) ) {
	    push(@chosen_src_gene_ids, $gene_id);
	}
	print OUT ">$frag_id-$gene_id ";
	if ($linefrag =~ /$hdr_remainder/) {
	    print OUT "$1\n";
	} else {
	    print OUT "\n";
	}
	$num_frags_written++;	
	$copy=1;
    }
    close(OUT);
    close(INFRAG);
    
    if ( (! $seq_outfile) and $seqfile ) {
	if ($seqfile =~ /$pattern/) {	    
	    $ext = $pep_fasta_ext;
	} else {
	    $ext = $metasim_fasta_ext;
	}	
	$seq_outfile = get_basename($seqfile, @exts).$final_src_seqs_name.scalar(@chosen_src_gene_ids).$ext;
    }
    if ($outdir and $seq_outfile) {
	$seq_outfile = File::Spec->catfile($outdir, $seq_outfile);	
    }    
    my $num_src_seqs=0;
    if ($seqfile and (-e $seqfile)) {
	open (OUTSEQ, ">$seq_outfile");
	open(INSEQ, $seqfile);
	$copy=0;
	while (my $line=<INSEQ>) {
	    if ($copy) {
		if ($line =~ $basic_header) {
		    $copy=0;
		} else {
		    print OUTSEQ $line;
		    next;
		}
	    }
	    if (scalar( grep {$line =~ $_} @chosen_src_gene_ids)) {
		$num_src_seqs++;
		print OUTSEQ $line;
		$copy=1;
	    }
	}
    }
    close(INSEQ);
    close(OUTSEQ);
    if ($num_src_seqs != scalar(@chosen_src_gene_ids)) {
	die "In filter_reads: problem finding all source sequences in file $seqfile!\n";
    }
    
    return ($outfile, $seq_outfile, $num_frags_written, $num_src_seqs);
}


sub cut_to_ids($) {
    my $seq_hdrs_ref = shift(@_);
    my $pattern = qr{^>(\S+)\s};
    my @ids;
    foreach my $seq_hdr (@$seq_hdrs_ref) {
	if ($seq_hdr =~ /$pattern/) {
	    my $id = $1;
	    push(@ids, $id);
	}
    }
    return \@ids;
}


sub printReturnVal($) {
    my $val = shift(@_);

}

sub printErrorLog ($;$) {
    my $errorLog = shift(@_);
    my $fh = shift( @_);
    if (! defined($fh)) {
	$fh = \*STDOUT;
    }
    print $fh qq{MetaSim 'error.log':\n};
    open(IN, $errorLog) or die "Cannot open MetaSim error log $errorLog : $!\n";
    while (my $line = <IN>)
    {   print $fh $line;
    }
    print $fh "\n";
}

1;
