#!/usr/bin/perl
# S.J. Riesenfeld
# MetaPASSAGE module: seqsPad
# Updated: Nov 2010
#
# Prints to $file_out every sequence in $file_in but with pads of 'N's
# of length $pad_size on each end of the sequence.

package seqsPad;
	       
use strict;
use warnings;

use simPipeVars qw(:simpipe :metasim);
use List::Util qw(min);

require Exporter;
our @ISA = qw(Exporter);
our @EXPORT_OK = qw(seqs_pad seqs_remove_pad);

sub seqs_pad ($;$$$) {
    
    my ($file_in, $pad_size, $file_out) = @_;
    if (!defined($file_in)) {
	die "Please provide a file name to seqs_pad!\n";
    }
    if (!defined ($pad_size)) {
	$pad_size = $seq_pad_size_default;
    }
    if (!defined($file_out)) {
	if ($file_in =~ /(.+)(\.[^\.]+)/) {
	    $file_out = $1.$pad_filename_part.$2;
	} else {
	    $file_out = $file_in.$pad_filename_part;
	}
    }

    # avoid dependence on BioPerl for now
    open(IN, "$file_in") or die "Cannot open $file_in for read: $!\n";
    open(OUT, ">$file_out") or die "Cannot open $file_out for write: $!\n";
    my $pad = 'N'x$pad_size;

    my ($hdr, $seq);
    my $first_hdr = 1;
    my $line_len;
    while (my $line = <IN>) {
	if ($line =~ /^>/) {
	    if (! $first_hdr) {
		print OUT $hdr;
		$seq = $pad.$seq.$pad;
		write_seq($seq, $line_len, \*OUT);
	    } else {
		$first_hdr=0;
	    }
	    $seq='';
	    $hdr = $line;
	} else {   
	    chomp($line);
	    if (!defined($line_len) or (length($line) > $line_len)){
		$line_len = length($line);
	    }
	    $seq .= $line;
	}
    }
    
    if (! $first_hdr ) {
	$seq = $pad.$seq.$pad;
	print OUT $hdr;
	write_seq($seq, $line_len, \*OUT);
    }
    # print "$dropped_ctr sequences dropped.\n";
    close(OUT);
    close(IN);
    return ($file_out, $pad_size);
	      }
   

## If 'N' is found in middle of sequence, sequence is dropped, as pads
## should only appear on ends
sub seqs_remove_pad ($;$) {
    my ($file_in, $file_out) = @_;
    my $bad_seqs=0;
    if (!defined($file_in)) {
	die "Please provide a file name to seqs_remove_pad!\n";
    }
    if (!defined($file_out)) {
	if ($file_in =~ /(.+)(\.[^\.]+)/) {
	    $file_out = $1.$nopad_filename_part.$2;
	} else {
	    $file_out = $file_in.$nopad_filename_part;
	}
    }

    # avoid dependence on BioPerl for now
    open(IN, "$file_in") or die "Cannot open $file_in for read: $!\n";
    open(OUT, ">$file_out") or die "Cannot open $file_out for write: $!\n";
    
    my $bad_seq_pattern = qr/[^N]+N+[^N]+/;
    my $ok_seq_pattern = qr/^N*([^N]*)N*$/;
    my ($hdr, $seq);
    my $first_hdr = 1;
    my $line_len;
    while (my $line = <IN>) {
	if ($line =~ /^>/) {
	    if (! $first_hdr) {
		if ($seq =~ /$bad_seq_pattern/) {
		    # warn "Bad seq is: $seq\n";
		    $seq = '';
		}
		$seq =~ s/$ok_seq_pattern/$1/;
		if ($seq eq '') {  # catches seqs that match $bad_seq_pattern and also seqs that have only 'N's
		    $bad_seqs++;
		} else {
		    print OUT $hdr;
		    write_seq($seq, $line_len, \*OUT);
		}
	    } else {
		$first_hdr=0;
	    }
	    $seq='';
	    $hdr = $line;
	} else {   
	    chomp($line);
	    if (!defined($line_len) or (length($line) > $line_len)){
		$line_len = length($line);
	    }
	    $seq .= $line;
	}
    }
    
    if (! $first_hdr) {
	if ($seq =~ /$bad_seq_pattern/) {
	    # warn "Bad seq is: $seq\n";
	    $seq = '';
	}
	$seq =~ s/$ok_seq_pattern/$1/;
	if ($seq eq '') {
	    $bad_seqs++;
	} else {
	    print OUT $hdr;	    
	    write_seq($seq, $line_len, \*OUT);
	}
    }
    # print "$bad_seqs sequences dropped during pad removal\n";
    close(OUT);
    close(IN);
    return ($file_out, $bad_seqs);
		     }

sub write_seq ($$$) {
    my ($seq, $line_len, $fh)= @_;
    if ($line_len < 0) {
	die "Line length must be nonnegative.\n";
    }
    while ($seq) {
	my $len = min(length($seq), $line_len);
	print $fh substr($seq, 0, $len). "\n";
	$seq = substr($seq, $len);
    }
	       }

1;
