#!/usr/bin/perl
 
use strict;
use warnings;
use Bio::AlignIO;

# Iker Irisarri @ University of Konstanz, Jan 2015
# Script for processing multi phylip files generated for multilocus bootstrap in MP-EST
# Given an input file, it prints all lines in the file after a given number of alignments have been counted
# E.g., some phyml runs would not finish completely and stop after 2000 trees have been estimated. This script allows extracting the alns after aln #2000 in order to estimate only the trees that are missing
# I tried using Bio::AlignIO but it does not process the alns in order so it was more complicated with bioperl

my $usage = "multiphylip_get_alns.pl multi_phylip_file.phy tree_number > stdout\n";
my $infile = $ARGV[0] or die $usage;
my $aln_num = $ARGV[1] or die $usage;

open (IN, "<", $infile);

my $count = 0;

while ( my $line = <IN> ) {
    
    chomp $line;
    
    # match begining of each aln, in this case all have 22 taxa
    if ( $line =~ /^22 / ) {
	
	$count++;
    
    }

    if ( $count > $aln_num ) {

	print "$line\n";

    }
}


print STDERR "done!\n";
