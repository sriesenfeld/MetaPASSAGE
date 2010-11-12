#!/usr/bin/perl
#
# S.J. Riesenfeld
# Last updated: Nov. 2010
# MetaPASSAGE module: refSeqs
# 
# Module for, among other things: 1) sampling a given number of gene
# sequences from an input set of sequences in order to create a
# Reference Database, possibly with the objective of maximizing the
# phylogenetic diversity represented by a subset of the Reference
# Database; 2) sampling a given number of peptide sequences from a
# database; writing the corresponding DNA sequences if desired
# (assuming a version of the database that contains the corresponding
# DNA sequences has been created, probably via
# retrieveSeqs::retrieve_seqs()); these functions may also have other
# applications.
#
# When using with AMPHORA, for now, these functions ignores all
# AMPHORA Reference Sequences whose ids begin with the word 'MERGED'
# because the corresponding DNA sequences cannot be retrieved easily
# from the NCBI database.
#
# The functions are designed not to select or sample sequences whose
# names contain words listed in @taxa_to_avoid (set in
# simPipeVars.pm), which has by default Mycoplasma and Candidatus.

package refSeqs;

use strict;
use warnings;

require Exporter;
our @ISA = qw(Exporter);
our @EXPORT_OK = qw(make_reference_db sample_seqs get_random_subset get_gene_ids_from_headers get_sequence_headers);

use File::Spec;

use simPipeVars qw{:amphora :misc_file_ext :simpipe :maxpd};
use maxPD qw(maxPD_list);


=head2 process_maxpd_output

 Usage: process_maxpd_output(<arg1>, <arg2>, <arg3>);

 Function: Chooses top <arg2> (if possible) species names from file
   <arg3> (which should be the output of a call to maxPD_list, e.g.,
   the results of calling maxPD_list on, say, AMPHORA Genome tree),
   such that the species names do not contain words in @taxa_to_avoid
   (set in simPipeVars.pm to contain 'Mycoplasma' and
   'Candidatus'). Returns a reference to an array of the corresponding
   sequence headers -- only species names with non-zero maxPD values
   may be returned -- and also a reference to an array of the
   remaining sequence headers that contain only taxa that appear in
   the maxPD output list and do not have names with words in
   @taxa_to_avoid.  (This gives a means of filtering which species
   appear in the output at all.)  Converts the species names so they
   have spaces instead of underscores in them.

 Example: 
  # 'maxpd_out_file' already contains the output of a run of
  # maxPD_list(); @seq_headers is an array of fasta sequence headers;
  ($selected_hdrs_ar, $remaining_hdrs_ar) =
  process_maxpd_output(\@seq_headers, 100, 'maxpd_out_file'); 

 Returns: An array containing two references to arrays; the first is a
   reference to an array containing the selected headers; the second
   is a reference to an array containing the remaining headers (that
   do not have names in @taxa_to_avoid and also appear in <arg3>).

 Args: The input is an array of three arguments:
   arg 1: a reference to an array of fasta sequence headers;
   arg 2: a number headers to select; (may be left undefined or
     negative, in which case the whole list of sequence headers is
     returned, simply without the taxa that have names in
     @taxa_to_avoid)
   arg 3: a file containing the output of a run of maxPD_list(), must
      be defined if <arg2> is non-negative; 

=cut

sub process_maxpd_output($$$) {
    my ($seq_headers_ar, $maxpd_num, $maxpd_output) = @_;
    my @taxa;
    my @skipped_taxa;
    my @maxpd_seq_headers;
    my @zero_maxpd;
    my @rem_headers;
    my @avoided_hdrs;    
    my $pos_maxpd_taxa=0;
    my @all_taxa;
    
    if (defined($maxpd_num) and ($maxpd_num >=0)) {
	open(MAXPD, $maxpd_output) or die "Cannot open file $maxpd_output to read output of maxPD: $!\n";    
	while (my $line =<MAXPD>) {	
	    my @line_parts = split(' ', $line);
	    my $taxon = $line_parts[0];
	    $taxon = join(' ', split('_', $taxon));
	    if (!(defined($line_parts[1])) or ($line_parts[1] > 0)) {
		push(@taxa, $taxon);	
	    } else {
		push(@zero_maxpd, $taxon);
	    }
	}
	print "Total of ". scalar(@taxa). " taxa with non-zero maxPD value listed in the file $maxpd_output;\n".
	    "total of ".scalar(@zero_maxpd). " other taxa listed.\n";
	if ($maxpd_num > scalar(@taxa)) {
	    print "Can choose at most ".scalar(@taxa). " taxa with non-zero maxPD value!\n";
	    $maxpd_num = scalar(@taxa);
	}
	$pos_maxpd_taxa = scalar(@taxa);        
	@all_taxa = (@taxa, @zero_maxpd);
    }
    print "Filtering the ".scalar(@{$seq_headers_ar})." given sequence headers now.\n";
    if (!scalar(@all_taxa)) {	
	my @headers = @{$seq_headers_ar};
	my $num_avoided=0;
	foreach my $hdr (@headers) {
	    if ( grep {$hdr =~ $_} @taxa_to_avoid ) {
		next;
	    } else {
		push (@rem_headers, $hdr);
	    }
	}
	$num_avoided++;
	print "($num_avoided headers excluded because contain one of: ".join(', ', @taxa_to_avoid).".)\n";
    } else {
	# $maxpd_num defined and >= 0;
	# pick the top $maxpd_num of the taxa output by maxPD_list and find the matching sequence headers for them
	# return all remaining headers, without hits to @taxa_to_avoid
	foreach my $taxon (@all_taxa) {
	    if ($pos_maxpd_taxa >= 0) {
		$pos_maxpd_taxa--;
	    }
	    # my @matches = grep ($_ =~ $taxon, @{$seq_headers_ar});
	    my $name_pattern = join('.*', split(' ', $taxon));
	    # my @hdr_matches = grep {$_ =~ /[\[,\{]$name_pattern[\]\}]/ } @{$seq_headers_ar};
	    my @hdr_matches = grep {$_ =~ /$name_pattern/ } @{$seq_headers_ar};
	    if (! scalar(@hdr_matches))  {
		print "No matches for taxon \'$taxon\'!\n";
		while (my ($key, $value) = each(%alt_names)) {
		    my $key_flag=0;
		    my $name_pattern_alt;
		    if ($key =~ $name_pattern) {
			$name_pattern_alt = join('.*', split(' ', $value));
		    } elsif ($value =~ $name_pattern) {
			$key_flag=1;
			$name_pattern_alt = join('.*', split(' ', $key));
		    }
		    if ($name_pattern_alt) {
			$name_pattern = $name_pattern_alt;
			@hdr_matches = grep {$_ =~ /$name_pattern/ } @{$seq_headers_ar};
		    }
		    if (scalar(@hdr_matches)) {
			print "  But matched taxon under alternative name \'".($key_flag ? $key : $value)."\'.\n";
		    last;
		    }
		}
	    }
	    if (! scalar(@hdr_matches))  {	
		push(@skipped_taxa, $taxon);
		next;
	    } elsif (scalar(@hdr_matches) > 1) {
		print "More than one (".scalar(@hdr_matches).") sequence header matches for $taxon in $maxpd_output:\n".	    
		    "  Trying to reduce the number of matches.\n";
		my @alt_hdr_matches = grep {$_ =~ /[\[,\{]$name_pattern[\]\}]/ } @{$seq_headers_ar};	    
		if (scalar(@alt_hdr_matches)) {
		    @hdr_matches = @alt_hdr_matches;
		}
		if (scalar(@hdr_matches > 1)) {
		    print "  Still multiple matches -- using at most one for the reference database.\n";
		} else {
		    print "  Succeeded in getting one best match.\n";
		}
	    }
	    for my $i (0..$#hdr_matches) {
		my $hdr = $hdr_matches[$i];
		if ( grep {$hdr =~ $_} @taxa_to_avoid ) {
		    if (! grep {$hdr eq $_} @avoided_hdrs) {
			push(@avoided_hdrs, $hdr);
		    }
		} elsif ( ($i==0) and ( scalar(@maxpd_seq_headers) < $maxpd_num ) and ($pos_maxpd_taxa >= 0)) {
		    if (! grep {$hdr eq $_} @maxpd_seq_headers ) {		
			push(@maxpd_seq_headers, $hdr);
		    }
		} elsif (! grep {$hdr eq $_} @rem_headers) {
		    push(@rem_headers, $hdr);
		}
	    }
	}
	print "A total of ".scalar(@maxpd_seq_headers)." sequence headers chosen for reference database via (positive) maxPD.\n"; 
	
	print "A total of ".scalar(@rem_headers). " remaining headers that contain a taxon from the list in $maxpd_output.\n";
	print "(".scalar(@avoided_hdrs)." headers not included in either count because they are one of:  "
	    .join(', ', @taxa_to_avoid).".)\n";
	print "A total of ".scalar(@skipped_taxa)." taxa had no matching sequence headers.\n";
    }
    return (\@maxpd_seq_headers, \@rem_headers);
}

=head2 make_reference_db

 Usage: make_reference_db(%args);
 
 Function: Creates a database (single file, or multiple files if both
   DNA and peptide sequence files are present) of sequences, based on
   the input sequences and options.

 Example: 
   # Use with AMPHORA:
   make_reference_db( '-n', 200, '-m', 25, '-g', 'rpoB', '-g',
       $gene_symbol, '-o', 'RefDB_rpoB', '-d',
       '/usr/local/AMPHORA/Reference/Sequence', '-p',
       'maxPD_out_file', '-t', 'AMPHORA_Genome_tree.tre', '-c', 570 );
   # Generic use:
   make_reference_db( '-n', 200, '-o', 'RefDB',
       '-f', $all_seqs_basename );

 Returns: Nothing.

 Args: Hash with the keys and values as follows. The values for keys
     '-n' and '-o' are required.  Aslo, the value for key '-f' or key
     '-g' must be defined.  The rest are optional.

   Key: '-n'; Value: Positive integer, total number of sequences to
     put in reference database;
   Key: '-g': Value: String, gene symbol use to name protein family in
     AMPHORA.  
   Key: '-f'; Value: <basename> (including path, but not extensions)
     of fasta file(s) containing source sequences for the Reference
     Database; at least one of <basename>.pep (peptide sequences) and
     <basename>.fna (corresponding DNA sequences) must exist;
     extensions for files (unaligned peptides and unaligned DNA
     sequences) are given in simPipeVars.pm; if option '-f' is
     specified, all AMPHORA-related options are ignored ('-s', '-d').
   Key: '-d'; Value: Directory contaning the AMPHORA Reference
     sequences; this option is ignored if option '-g' is not used; if
     not specified, it is assumed to be given by the simPipeVars.pm
     variable $amphora_ref_seq_path or to be the current directory;
   Key: '-c'; Value: Number of taxa in the given tree for which maxPD
     should be computed (the taxa in the output are the taxa used to
     populate the reference database, both by maxPD value and by
     random sampling);
   Key: '-m'; Value: Positive integer, number of sequences that should
     come in order of maxPD contribution; requires key '-t' to be set
     or key '-p' to be set with a pre-existing file;
   Key: '-t'; Value: Name of file containing a Newick tree, to be fed
     to maxPD_list;
   Key: '-p'; Value: Name of file (pre-existing or to be created) that
     contains (or will contain) the output of maxPD_list on a tree
     where each leaf label contains a taxon appearing in a header of
     the sequences input via option '-g' (for AMPHORA) or '-f'; if
     '-t' is set, then the output of maxPD_list is dumped in a file
     with this name; otherwise, this file is assumed to already
     contain the correct output;

   ### STILL TO BE IMPLEMENTED (BELOW)
   Key: '-a'; Value: Integer; if defined and nonzero, then the program
     checks for pre-existing files of the name specified by key '-o';
     if they exist, they are read and new sequences are only appended
     to the old db to create a db of size specified by key '-n';
     (currently not implemented).

=cut

sub make_reference_db(@) {
    my %options = @_;
    my ($num, $maxpd_num, $maxpd_compute_num, $gene_symbol, 
	$refdb_basename, $maxpd_outfile, 
	$maxpd_input_tree, $dir, $file_seqs_basename);
    my $append_flag = $options{'-a'};
    $num = $options{'-n'};
    $maxpd_num = $options{'-m'};
    $maxpd_compute_num = $options{'-c'};
    $gene_symbol = $options{'-g'};
    $refdb_basename = $options{'-o'};
    $maxpd_outfile = $options{'-p'};
    $maxpd_input_tree = $options{'-t'};
    $file_seqs_basename = $options{'-f'};
    $dir = $options{'-d'};
    
    unless ($gene_symbol or $file_seqs_basename) {die "Require a gene or file of sequences to be specified.\n";}
    
    my ($seqs_file_pep, $seqs_file_dna, $outfile_pep, $outfile_dna);
    if ($append_flag) {
	die "make_reference_db: The append option is currently not implemented!\n";
    }
    if (! $dir) {       
	$dir = $amphora_ref_seq_path;
    }
    if (defined($gene_symbol) and !(-d $dir)) {
	die "Cannot find directory $dir containing the AMPHORA Reference Sequences.\n"; 
    }        
    if (!defined($file_seqs_basename)) {
	$file_seqs_basename = $gene_symbol;
    }
    $seqs_file_pep = $file_seqs_basename.$pep_fasta_ext;
    $seqs_file_dna = $file_seqs_basename.$dna_fasta_ext;
    if (!defined($file_seqs_basename)) {
	$seqs_file_pep = File::Spec->catfile($dir, $seqs_file_pep);    
	$seqs_file_dna = File::Spec->catfile($dir, $seqs_file_dna);
    }    
    $outfile_pep = $refdb_basename . $pep_fasta_ext;
    $outfile_dna = $refdb_basename . $dna_fasta_ext;
    
    unless ($num > 0) {die "Require a positive number of sequences to be specified.\n";}
    if ($maxpd_num) {
	unless (($maxpd_num > 0) and ($num >= $maxpd_num)) {
	    die "Require a positive number of sequences less than $num to be chosen via maxPD.\n";
	}
	if (! $maxpd_outfile) {
	    if ($maxpd_input_tree) {
		die "Please specify a filename for the output of maxPD calculations.\n";
	    } else {
		die "Please specify an input tree to feed maxPD_list or a file containing maxPD_list output.\n";
	    }
	}
    }
    
    if (!$maxpd_compute_num) {
	$maxpd_compute_num = $maxpd_default_num;
    }
    unless ($refdb_basename) {die "Require base name for reference database output files to be specified.\n";}
    
    if ($maxpd_input_tree) {
	maxPD_list('-t', $maxpd_input_tree, '-n', $maxpd_compute_num, '-o', $maxpd_outfile);
	if (!(-e $maxpd_outfile)) {
	    die "Problem writing maxPD output to $maxpd_outfile!\n";
	}
	print "MaxPD list calculated for input tree $maxpd_input_tree; results written to file $maxpd_outfile.\n";	
    } elsif ($maxpd_num) {
	print "Using existing maxPD_list output file $maxpd_outfile.\n";
    }
    
    # choose a bunch of sequences via maxPD and the rest randomly 
    # (none of these sequences have descriptions containing words in @taxa_to_avoid and all of them are from taxa 
    # listed in the maxpd output file);
    # my @cur_ref_ids;
    # my @cur_ref_headers;
    my @seq_headers;
    if (-e $seqs_file_dna) {
	@seq_headers = get_sequence_headers($seqs_file_dna);
	print "Getting sequence headers from file $seqs_file_dna.\n";
    } else {
	@seq_headers = get_sequence_headers ($seqs_file_pep);
	print "Getting sequence headers from file $seqs_file_pep.\n";
    } 
    my ($chosen_headers_ar, $rem_headers_ar, $spots_left);
    if ($maxpd_num) {
	($chosen_headers_ar, $rem_headers_ar) = 
	    process_maxpd_output(\@seq_headers, $maxpd_num, $maxpd_outfile);
	$spots_left = $num - scalar(@{$chosen_headers_ar});
    } else {
	$spots_left = $num;
	$rem_headers_ar = \@seq_headers;
    }
    my @rem_ids = get_gene_ids_from_headers($rem_headers_ar);
    my @rem_avail_hdrs;

    if ($spots_left > 0) {
	my $num_rem = scalar(@{$rem_headers_ar});
	if ($num_rem == $spots_left) {
	    push(@{$chosen_headers_ar}, @{$rem_headers_ar});
	} elsif ($num_rem < $spots_left) {
	    die "Reference database size parameter is too large:\n".
		"Number of sequences to sample randomly (without replacement) for reference database can be at most $num_rem.\n";
	} else {
	    my @random_sample = get_random_subset($spots_left, $rem_headers_ar);
	    print ''.scalar(@random_sample)." sequences chosen randomly (without replacement) from remaining sequences.\n";
	    push(@{$chosen_headers_ar}, @random_sample);
	}
    } 
    
    my $chosen_gene_ids_ar = get_gene_ids_from_headers($chosen_headers_ar);
    if (-e $seqs_file_dna) {
	write_seqs($chosen_gene_ids_ar, $seqs_file_dna, $outfile_dna);
	if (-e $outfile_dna) {
	    print "DNA reference database written to file $outfile_dna.\n";	    
	} else {
	    die "Could not write reference database $outfile_dna!\n";
	}	
    } 
    if (-e $seqs_file_pep) {
	write_seqs($chosen_gene_ids_ar, $seqs_file_pep, $outfile_pep);
	if (-e $outfile_pep) {
	    print "Peptide reference database written to file $outfile_pep.\n";
	} else {
	    die "Could not write reference database $outfile_pep!\n";
	}
    }
    return;
}

=head2 sample_seqs

 Usage: sample_seqs(@args);

 Function: Samples sequences randomly from a fasta file of sequences,
   but so that the correct number are taken from the Reference
   Database, if a file containing those headers (with or without
   sequences) is provided.  
 Note: Right now, this function avoids sampling @taxon_to_avoid, set
   in simPipeVars.pm.

 Example: 
   my $sample_fna = sample_seqs(100, 'sample', 'seqs', 0, 'RefDB');

 Returns: The name of the output file containing DNA sequences.

 Args: Array of 5 arguments; the last 2 are optional.
  arg 1: positive integer = size of the sample;
  arg 2: basename for output files;
  arg 3: base name of fasta files that contain unaligned sequences and
         unaligned dna sequences, which have extensions '.pep' and
         '.fna', respectively; (the dna sequence file is the only one
         that *must* exist);
  arg 4: (optional) integer = the number of sequences that should be
         sampled uniformly at random from the reference database (arg
         5); the remaining number (arg 1 minus arg 4) are sampled from
         the database files given by arg 3); must be used in
         conjuction with arg 5;
  arg 5: (optional) basename of "reference database" file (or file
         containing those headers), for use with arg 4.

=cut

sub sample_seqs (@) {
    my ($sample_size, $outfile_basename, $db_basename, 
	$num_from_ref_db, $ref_db_hdrs_file) = @_;    
    unless ($sample_size>0) { 
	die "Specify positive sample size.\n"; }

    my $db_file_pep = $db_basename . $pep_fasta_ext;
    my $db_file_dna = $db_basename . $dna_fasta_ext;
    my @all_headers = get_sequence_headers($db_file_dna);  
    my $all_ids_ar = get_gene_ids_from_headers(\@all_headers);
    
    my (@ref_headers, @rem_headers);
    if ($ref_db_hdrs_file) {
	@ref_headers=get_sequence_headers($ref_db_hdrs_file);
    }
    
    my @copy_all_ids = @{$all_ids_ar};
    foreach my $header (@all_headers) {
	my $id = shift(@copy_all_ids);
	if ($header !~ $id) {	    
	    die "In sample_seqs:  Problem coordinating headers and ids.\n";
	}
	if ( (! grep ($header =~ $_, @taxa_to_avoid)) and (! grep( $_ =~ $id, @ref_headers))) {
	    push(@rem_headers, $header);
	}
    }
    
    my ($outfile_pep, $outfile_dna);
    $outfile_pep = $outfile_basename . $pep_fasta_ext;
    $outfile_dna = $outfile_basename . $dna_fasta_ext;
    
    my (@chosen_ref_headers, @chosen_rem_headers, @chosen_headers);
    if (defined($num_from_ref_db) and ($num_from_ref_db >= 0)) {
	if ($num_from_ref_db > 0) {
	    print "Sampling $num_from_ref_db sequences from the ".scalar(@ref_headers). " sequences in the Reference Database.\n";
	    @chosen_ref_headers = get_random_subset($num_from_ref_db, \@ref_headers);
	}
	
	my $num_rem = $sample_size - $num_from_ref_db;
	if ($num_rem > 0) {
	    print "Sampling $num_rem sequences from the ".scalar(@rem_headers). " available sequences outside the Reference Database.\n";
	    @chosen_rem_headers = get_random_subset($num_rem, \@rem_headers);
	}
	@chosen_headers = (@chosen_ref_headers, @chosen_rem_headers);
    } else {
	@all_headers = (@ref_headers, @rem_headers);
	@chosen_headers = get_random_subset($sample_size, \@all_headers);
    }
    my $chosen_ids_ar = get_gene_ids_from_headers(\@chosen_headers);
    if (-e $db_file_pep) {
	write_seqs($chosen_ids_ar, $db_file_pep, $outfile_pep);
    }
    write_seqs($chosen_ids_ar, $db_file_dna, $outfile_dna);
    # unless ($no_dna_flag) {
    #    my @chosen_gene_ids = cut_to_gene_ids(\@chosen);
    #    write_seqs(\@chosen_gene_ids, $db_file_dna, $outfile_dna, 1);
    # }
		 
    return $db_file_dna;
}

# first arg: positive integer = SIZE of random subset to be returned;
# second arg: reference to array that contains elements of set to be sampled
# returns an array of length SIZE containing a random subset of the elements.
sub get_random_subset ($$) {
    my ($size, $array_ref) = @_;
    my $range = scalar(@$array_ref);
    if ($range == $size) {
	my @ret_array=@{$array_ref};
	return @ret_array;
    }
    unless ($range > $size) {die "In get_random_subset:  Size of array must be at least the desired size of random subset.\n";}
    my @set = (0..($range-1));
    my @chosen;
    foreach my $count (0..($size-1)) {
	my $random_num = int(rand($range-$count-1));	
	swap($count, ($random_num+$count+1),\@set);	
	my $index = $set[$count];
	# print "count: $count; at index $index: ". $$array_ref[$index]."\n";
	push(@chosen, $$array_ref[$index]);
    }
    return @chosen;
}

# first arg: positive integer = index
# second arg: positive integer = index
# third arg: reference to an array
sub swap ($$$) {
    my ($index1, $index2, $array_ref) = @_;
    unless ( ($index1 < scalar(@$array_ref)) and ($index2 < scalar(@$array_ref)) ) {
	die "Trying to reference element outside array boundary.\n";}
    # print "Swapping elements at indices $index1 and $index2.\n";
    my $temp = $$array_ref[$index1];
    $$array_ref[$index1]=$$array_ref[$index2];
    $$array_ref[$index2]=$temp;
    # print "at $index1: ". $$array_ref[$index1]."; at $index2: ".
	$$array_ref[$index2]."\n";
}

sub get_gene_ids_from_headers($) {
    my $headers_ar = shift(@_);
    my @gene_ids;
    my $pattern = qr{^>(\w+)};
    foreach my $header (@{$headers_ar}) {
	if ($header =~ /$pattern/) {
	    my $id = $1;
	    push(@gene_ids, $id);
	} else {
	    die "Badly formatted header!\n";
	}	
    }
    return \@gene_ids;
}

sub cut_hdr_to_id($) {
    my $header = shift();
    my $pattern = qr{^>(\w+)};
    if ($header =~ /$pattern/) {
	my $id = $1;
	return $id;
    } else {
	die "Badly formatted header!\n";
    }
		  }

# first arg: reference to array of sequence gene ids
# second arg: filename of fasta file in which to find sequences
# third arg: filename in which to write sequences
sub write_seqs($$$;$) {
    my ($seqs_array_ref, $infile, $outfile) = @_;
    my $header = qr{^>};
    open(INF, "$infile") or die "Cannot open $infile: $!\n";
    open(OUTF, ">$outfile") or die "Cannot open $outfile: $!\n";
    # print "Writing sequences from $infile to $outfile...\n";
    my $line = <INF>;
    my @hits;
    # print "Expecting to write ". scalar(@{$seqs_array_ref})." sequences.\n";
    while($line) {	  
	my @new_hits = grep($line =~ "^>$_", @$seqs_array_ref);
	if (scalar(@new_hits)) {
	    if (scalar(@new_hits)>1) {
		print "More than one gene id in list matches: ".join(' ', @new_hits)."\n";
		print "Writing from line: $line.\n";
	    }
	    push(@hits, @new_hits);	    
	    do {		
		# print "$line";
		print OUTF $line;
		$line = <INF>;
	    } while ($line and !($line =~ /$header/));
	} else {
	    $line = <INF>;
	}
    }
    # print ''.scalar(@hits)." sequences identified in file $infile and written to file $outfile.\n";
    my @nonhits;
    if (scalar(@hits) < scalar(@{$seqs_array_ref})) {
	@nonhits = @{$seqs_array_ref};
	foreach my $hit(@hits) {
	    @nonhits = grep($_ ne $hit, @nonhits);
	}
	if (!scalar(@nonhits)) {
	    print ''.scalar(@nonhits)." gene ids for which sequences were not identified: ".join(' ', @nonhits)."\n"; 
	}
    }
    close(INF);
    close(OUTF);
	       }

# first arg: filename for fasta file of sequences with headers
# returns an array of the sequence (gene) ids
sub get_sequence_headers($) {
    my $infile = shift(@_);
    open(INF, "$infile") or die "Cannot open $infile : $!\n";
    # my $header = '^>([\S]*)\s.*';
    # my $header = '^>(\w+)-'; 
    my $header = qr{^>};
    my @seq_headers;
    while( my $line = <INF> ) {
	if ($line =~ /$header/) {
	    push(@seq_headers, $line);
	}
    }
    close(INF);
    return @seq_headers;
}


1;
