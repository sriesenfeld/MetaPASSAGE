#!/usr/bin/perl

# S. Riesenfeld

use strict;
use warnings;

my ($file_in, $file_out) = @ARGV;
open(IN, "$file_in") or die "Cannot open $file_in: $!\n";
open(OUT, ">$file_out") or die "Cannot open $file_out: $!\n";
while (my $line = <IN>) {
    if ($line !~ /^>/) {
	$line =~ s/T/U/g;
	$line =~ s/t/u/g;
    }
    print OUT $line;
}
