#!/usr/bin/perl

use strict;
use warnings;
use Bio::TreeIO;
use Data::Dumper;
use Statistics::Descriptive;

# Iker Irisarri, University of Konstanz. Feb 2016
# simple script to get summary stats from tree branch lengths
# originally for identifying genes with higher and lower among-lineage rate heterogeneity
# using single-gene trees (after outgroup trimming for MP-EST)
# If used in a for loop, it creates a table for all trees


my $usage = "get_BL_stats_from_tree.pl tree monophyly_file outgroup> stdout\n";
my $in_tree = $ARGV[0] or die $usage;
my $outgroup = "Callorhinc";

# read in trees
my $treeio = new Bio::TreeIO(-file   => "$in_tree",
						     -format => "newick"
						    );

my $total_length = 0;
my @branch_lengths = ();
my $stat = Statistics::Descriptive::Full->new();

# print out heading

while( my $tree = $treeio->next_tree ) {


	# root the tree
    my $root = $tree->find_node( -id => $outgroup );
    $tree->reroot( $root );

	# get total branch length
	$total_length = $tree->total_branch_length;		
	
	for my $node ( $tree->get_nodes ) {
	
		my $length = $node->branch_length;
		# scape undefined branch lengths (root)
		next if (!defined $length);
		push (@branch_lengths, $length);
	}

	# calculate stats from @branch_lengths
	$stat->add_data(@branch_lengths);
	my $mean = $stat->mean();
	my $sd = $stat->standard_deviation();
	
	# print out stuff
	print "File\tmean_BL\tsd_BL\ttotal_BL\n";
	print "$in_tree\t$mean\t$sd\t$total_length\n";
	
}
	
print STDERR "\nDone!\n\n";	

