#!/usr/bin/perl
# S. Riesenfeld
# Last updated: 26 May 2009
# Updated: 8 April 2009
# MetaPASSAGE module: blast2trans
#
# Module for running blastx or blastn on reads to figure out what
# frame to translate them in or whether to reverse complement them
# (depending on whether the blast database is peptides or DNA).

package blast2trans;

use strict;
use warnings;
require File::Temp;
use File::Temp;
use List::Util qw(min);

require Exporter;
our @ISA = qw(Exporter);
our @EXPORT_OK = qw(run_blast_transeq); # copyNextSeq);

# use File::Temp qw(:mktemp);

use simPipeVars qw(:seqtypes :misc_file_ext :simpipe :blast);

=head2 run_blast_transeq

 Usage: run_blast_transeq(@args);

 Function: Runs NCBI blastx or blastn on all the sequences in the
   given file with the given input parameters.  Then parses the output
   file to find the frame or strand information for each sequence,
   compiles that data, and writes it to a file.  Then calls EMBOSS
   transeq to translate each sequence in the frames that have the best
   expect values, or reverse complements sequences according to the
   strand with the best expect value.  Saves those translations or
   orientations in a single file.  If the sequence is from Mycoplasma,
   then the genetic code for Mycoplasma is used with transeq.
   Otherwise, the bacterial genetic code is used for translation.

 Example: 
   ($trans_file, $frames_file) = run_blast_transeq('reads.fna',
      'outdir', ['-d', 'rpoB'], 'Protein');
      							 
 Returns: An array of two or three file names; the first is the file
   of translated or oriented sequences, the second is the file
   containing the frames or strand data for each query sequence; and
   the last is the blast output file, if it is saved.

 Args: An array of 5 arguments; the last 2 are optional.
   arg 1: name of file contaiining query sequences.
   arg 2: name of output directory.
   arg 3: a reference to an array containing options for blast,
       normally excluding the file output name, but including the
       blast database name (option '-d')
   arg 4: (optional) type of sequences in blast database ('Protein',
       etc, defined in simPipers.pm), default is 'Protein'.
   arg 5: (optional) a flag; if defined and non-zero, then the blast 
       file is saved; otherwise it is deleted.

=cut

sub run_blast_transeq ($$$;$$) {
    my ($in_file, $out_dir, $blast_args_ar, $type, $save) = @_;
    if (!defined($type)) {
	$type = PROTEIN;
    }
    # if $type is not PROTEIN then 'frames' should be thought of as 'strand' for the variable names below

    unless (-d $out_dir) { die "Cannot find directory $out_dir for writing output of blast.\n"; }
    
    open(IN, $in_file) or (die "Cannot open file $in_file of sequence files.\n");
    my $basename = get_basename($in_file);
    my $blast_outfile;
    if (PROTEIN) {
	$blast_outfile=$basename.$blastx_file_ext;
    } else {
	$blast_outfile = $basename.$blastn_file_ext;
    }
    $blast_outfile = File::Spec->catfile($out_dir, $blast_outfile);
    
    my $cmd;
    if ($type eq PROTEIN) {
	$cmd = 'blastall -p blastx -i ' . $in_file . ' -o '. $blast_outfile .' '. join(' ', @{$blast_args_ar});
    } else {
	$cmd = 'blastall -p blastn -i ' . $in_file . ' -o '. $blast_outfile .' '. join(' ', @{$blast_args_ar});
    }
    eval {
	system ($cmd);
    };
    if ($@) {
	die "Error running blastall: $@\n";
    }
    
    open(BL_OUT, "<$blast_outfile") or die "Cannot open $blast_outfile for reading blast output: $!\n";

    my $frames_file;
    if ($type eq PROTEIN) {
	$frames_file = File::Spec->catfile($out_dir, ($basename . $frames_file_ext));
    } else {
	$frames_file= File::Spec->catfile($out_dir, ($basename . $strand_file_ext));
    }
    open(FR_OUT, ">$frames_file") or die "Cannot open $frames_file for writing ".
	( ($type eq PROTEIN) ? "frame" : "strand")." output: $!\n";
    my $trOutfile;
    if ($type eq PROTEIN) {
	$trOutfile = File::Spec->catfile($out_dir, ($basename .$pep_fasta_ext));
    } else {
	$trOutfile = File::Spec->catfile($out_dir, ($basename .$orientation_name.$dna_fasta_ext));
    }
    open(TR_OUT, ">$trOutfile") or die "Cannot open $trOutfile for writing ".
	( ($type eq PROTEIN) ? "translated": "correctly oriented")." reads: $!\n";

    my $seq_count=0;
    my @frames = ('-1','-2','-3','+1','+2','+3');
    my @strand = (1, -1);   # 1 means same orientation, -1 means query and hit have opposite orientations
    
    # my ($temp_seq_dnafile, $temp_seq_pepfile);
    # ($temp_seq_dnafile = mktemp('tempXXXXX')) or die "Cannot create temporary file for blasting and translating sequences: $!\n";
    # ($temp_seq_pepfile = mktemp('tempXXXXX')) or die "Cannot create temporary file for blasting and translating sequences: $!\n";
    
    my ($temp_seq_dna_fh, $temp_seq_dna_file, $temp_seq_tra_fh, $temp_seq_tra_file);
    $temp_seq_dna_fh = File::Temp->new(TEMPLATE => 'tempXXXXX');
    $temp_seq_dna_file = $temp_seq_dna_fh->filename();
    $temp_seq_tra_fh = File::Temp->new(TEMPLATE => 'tempXXXXX');
    $temp_seq_tra_file = $temp_seq_tra_fh->filename();
    my $rc_count=0;
    my $no_rc_count=0;
    while ( copyNextSeq(\*IN, $temp_seq_dna_file,\$seq_count)) {
	
	my (%distinctFrames, %bestExpects, %names4bestExpects);
	
	my ($short_queryname, $long_queryname, $bestExpect);	
	if ($type eq PROTEIN) {
	    ($short_queryname, $long_queryname, $bestExpect) = 
		getFrames($type, \*BL_OUT, \@frames, \%distinctFrames, \%bestExpects, \%names4bestExpects);
	} else {
	    ($short_queryname, $long_queryname, $bestExpect) = 
		getFrames($type, \*BL_OUT, \@strand, \%distinctFrames, \%bestExpects, \%names4bestExpects);
	}
	my $trans_table;
	
	unless($short_queryname) {
	    next; }
	print FR_OUT "Query: ".$long_queryname."\n";

	# print "distinctFrames:\n";
	# print_hash(\%distinctFrames);
	
	my $framesTrans = 0;
	my $framescount = 0;
	for my $key ( keys %distinctFrames ) {
	    my $curframe;
	    if ( $distinctFrames{$key}) {
		$framescount++;
		print FR_OUT ''.(($type eq PROTEIN) ? "Frame: ": "Orientation: "). "$key\t Min expect: $bestExpects{$key}\n".
		    "\t\tSequence hit: ".$names4bestExpects{$key}."\n";	 
		# only translate in this frame if its expect value 
		# is < blast_max_ratio*(best expect over all frames) and also less than $best_max_expect		
		if ( ($bestExpects{$key}==0) or 
		     ( ($bestExpects{$key} < $blast_max_ratio*$bestExpect ) and ($bestExpects{$key} < $blast_max_expect) ) ) {
		    $framesTrans++;
		    if ($type eq PROTEIN) {
			if ($key =~ '^\+(.)') {
			    $curframe = $1;
			} else { $curframe = $key; }
			
			# figure out which table to use with transeq
			if (($names4bestExpects{$key} =~ /Mycoplasma/) or ($names4bestExpects{$key} =~ /Spiroplasma/)) {
			    $trans_table=4;	    
			} else {
			    $trans_table=11;
			}
			my $cmd2 = 'transeq '.$temp_seq_dna_file.' '.$temp_seq_tra_file.' -frame='.$curframe.
			    ' -table='.$trans_table.' -trim=Y -auto=Y -alternative=Y -clean=Y';
			eval{
			    system ($cmd2);
			};
			if ($@) {
			    die "Error running transeq: $@\n";
			}
			print FR_OUT "Translated with transeq using table $trans_table.\n";
			# FIX!!: If we translate more than once then, there will be multiple peptide reads with the same id
			# but different sequences! Need to come back to figure out how to deal with this. 
			# Append file of current translated sequence to file of translated sequences
			append_file (\*TR_OUT, $temp_seq_tra_file);
			# print "Read $short_queryname:  blastx identified $framescount frame" . 
			#    ($framescount ==1 ? '' : 's') .'; translated' .($framescount==1 ? '': (" in $framesTrans frame". 
			#							      	   ($framesTrans==1 ? '': 's'))). ".\n";
		    } else {
			if ($key < 0) {
			    open(TMP_IN, "$temp_seq_dna_file");
			    copyNextSeq(\*TMP_IN, $temp_seq_tra_file, \$rc_count, 1);
			    print FR_OUT "Sequence reverse complemented.\n";
			    close(TMP_IN);
			    # Append file of currently oriented sequence to file of previously oriented sequences
			    append_file (\*TR_OUT, $temp_seq_tra_file);
			} else {			    
			    # Append file of currently oriented sequence to file of previously oriented sequences
			    append_file (\*TR_OUT, $temp_seq_dna_file);
			    print FR_OUT "Sequence already correctly oriented.\n";
			}
		    }
		}		
		print FR_OUT "\n";
	    }
	}
	print FR_OUT "\n";

    }

    close(BL_OUT);
    
    if (!$save) {
	unlink ($blast_outfile);
    }
    # unlinking should happen automatically
    close(FR_OUT);
    close(TR_OUT);
    return ($trOutfile, $frames_file, $blast_outfile);
    
}


# Scans for the next set of frames (BLASTX) or strand (BLASTN) data
# and saves the information in given hashes; returns the short query
# name, the long query name, the best expect value over all frames or
# strands; returns 0 if fails.

sub getFrames ($$$$$$) {
    
    my ($type, $blOut_fhr, $framesr, $distinctFramesr, $bestExpectsr, $names4bestExpectsr) = @_;

    my $short_queryname=0;
    my $long_queryname;
    my $bestExpect = $blast_max_expect;

    my $copying_queryname=0;
    my $copying_hitname=0;
    my $curkey;
    my $curval;
    my $cur_hitname;

    undef ($curkey);
    undef ($curval);
    undef ($cur_hitname);

    # initialize hashes
    foreach my $key (@$framesr) {
	$distinctFramesr->{$key}=0;
	$bestExpectsr->{$key}=$blast_max_expect;
	$names4bestExpectsr->{$key}='';
    }
    my $process_key_flag=0;
    while(my $line= <$blOut_fhr>) {
	chomp($line);
	if ($copying_queryname or $copying_hitname) {
	    if ($line =~ /^[\s]*$/) {
		# found blank line, stop copying
		if ($copying_queryname) {
		    $long_queryname=join(' ', split(' ', $long_queryname));
		    $copying_queryname=0;
		} else {
		    $cur_hitname=join(' ', split(' ', $cur_hitname));
		    $copying_hitname=0;
		}
	    } else {
                # keep copying whole read name 
		if ($copying_queryname) {
		    $long_queryname.= (' '.$line);  
		} else {
		    $cur_hitname.=(' '.$line);
		}		    
	    }
	} elsif ($line =~ /Query= (.*)/) {
	    $copying_queryname = 1;
	    $long_queryname=$1;
	    undef ($cur_hitname);
            # get short name	    
	    if ($line =~ /Query= (.+)? [|].*/) {
		$short_queryname = "$1";
		# print "Looking at query $short_queryname...\n";
	    } elsif ($line =~ 'Query= (.+)') {
		$short_queryname = "$1";
		# print "Looking at query $short_queryname...\n";
	    } else {	
		print ("line: $line\n");
		$short_queryname = ''; 
		print "Query name not formatted as expected.\n";
	    }
	} elsif ($line =~ />(.*)/) {
	    $copying_hitname=1;
	    $cur_hitname = $1;
	} elsif ($line =~ /Expect = (\S+)/) {
	    $curval = $1;
	    $curval =~ s/^(e)(.*)/1e$2/;
	} elsif (($type eq PROTEIN) and ($line =~ /Frame = (\S+)/)) {
	    $curkey = $1;
	    $process_key_flag = 1;
	} elsif (($type ne PROTEIN) and ($line =~ /Strand = (\w+) \/ (\w+)/)) {
	    $process_key_flag = 1;
	    if ($1 eq $2) {
		$curkey = $framesr->[0];
	    } else {
		$curkey = $framesr->[1];
	    }
	} elsif ( (( $line =~ /BLASTX/ ) or ($line =~ /BLASTN/)) and ($short_queryname)) {
	    # stop processing the file
	    last;
	}	

	if ($process_key_flag) {
	    $process_key_flag=0;
	    if ( (!defined($curkey)) or (!defined($curval)) or (!defined($cur_hitname)) ) {
		print "Blast output not formatted as expected, in processing query $short_queryname.\n";
		print "  ".(($type eq PROTEIN) ? "Frame" : "Strand"). " data may have errors.\n";
	    } else {
		$distinctFramesr->{$curkey}=1;
		if ($bestExpectsr->{$curkey} > $curval) {
		    # keep track of best expect value and the name of the associated sequence in the database
		    # for the frame given by $curkey
		    $bestExpectsr->{$curkey} = $curval;
		    $names4bestExpectsr->{$curkey} = $cur_hitname;
		    # keep track of best expect value over all frames
		    if ($bestExpect > $curval) {
			$bestExpect = $curval;
		    }
		}
	    }
	    undef ($curval);
	    undef ($curkey);	    
	} 
    }
    return ($short_queryname, $long_queryname, $bestExpect);
}


# Look for the next fasta header in the file corresponding to the
# given filehandle (must be readable); copy the sequence (with the
# header) to the specified file (second argument); update the number
# of sequences (referenced by third argument); if the optional fourth
# argument $rc_seq_flag is defined and non-zero (by default it is not
# defined), then the sequence is reverse complemented during the copy;
# the fifth argument gives the sequence type (assumed by default not
# to be RNA).  Return 1 if success; returns 0 if reach end of first
# file without finding a sequence to copy.
sub copyNextSeq ($$$;$$) {
    
    my ($seqfile_fhr, $tempfile, $countr, $rc_seq_flag, $type) = @_;
    if (! defined ($type)) {
	$type = RNA2DNA;
    }
    open(TMP, ">$tempfile") or die "In copyNextSeq:  cannot open temporary file $tempfile: $!\n";
    my $copying =0;
    my $seq='';
    my $rc_seq;
    my $line_len=0;
    while(my $line = <$seqfile_fhr>) {
	if ($line =~ '>.*' ) { # found a fasta header
	    if ($copying) {  # stop copying
		seek($seqfile_fhr, -length($line), 1) # place the same line back onto the filehandle
		    or die "Cannot reset position in filehandle: $!\n";	     
		last;
	    } else { # start copying (should be header)
		$$countr++;
		$copying=1;		
		print TMP $line;	       
	    }	    
	} elsif ($copying) {	    
	    if ($rc_seq_flag) {
		my $seq_part = $line;
		chomp($seq_part);
		if (!$line_len) {
		    $line_len = length($seq_part);
		}
		$seq = $seq.$seq_part;
	    } else {
		print TMP $line;
	    }	    
	}
    }
    
    if ($rc_seq_flag) {
	# now rev comp whole seq and print it nicely
	# if ($type eq RNA) {
	    # $seq =~ tr/ACGUacgu/UGCAugca/;
	# } else { # DNA
	    $seq =~ tr/ACGTacgt/TGCAtgca/;
	# }
	my $rc_seq = reverse($seq);
	while (length($rc_seq) > 0) {
	    my $len=min(length($rc_seq), $line_len);
	    print TMP (substr($rc_seq, 0, $len)."\n");
	    $rc_seq = substr($rc_seq, $len);
	}
    }
    close(TMP);
    return $copying;
}


sub append_file ($$) {
    my ($bigfile_fhr, $file_to_append) = @_;
    open (IN_APPEND, $file_to_append) or die "In append_file:  cannot open $file_to_append to read: $!\n";
    while (my $line = <IN_APPEND>) {
	print $bigfile_fhr $line;
    }
    close IN_APPEND;
}


1;
