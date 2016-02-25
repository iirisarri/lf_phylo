#!/usr/bin/perl

use strict;
use warnings;
use Bio::TreeIO;
use Data::Dumper;
use Statistics::Descriptive;

# Iker Irisarri, University of Konstanz. Feb 2016
# Modif of get_BL_stats_from_tree.pl to get means of three clades
# simple script to get summary stats from tree branch lengths
# originally for identifying genes with higher and lower among-lineage rate heterogeneity
# using single-gene trees (after outgroup trimming for MP-EST)
# If used in a for loop, it creates a table for all trees


my $usage = "get_BL_stats_by_clade_from_tree.pl tree monophyly_file outgroup> stdout\n";
my $in_tree = $ARGV[0] or die $usage;
my $outgroup = "Callorhinc";

# define clades
my @teleosts = qw (Takifugu_r Danio_reri Oreochromi);
my @sarc_fishes = qw (Lepidosire Neoceratod Protopteru Latimeria_);
my @tetrapods = qw (Anolis_car Canis_fami Dasypus_no Gallus_gal Homo_sapie Loxodonta_ Macropus_e Meleagris_ Monodelphi Mus_muscul Ornithorhy Rana_chine Taeniopygi Xenopus_tr);

# read in trees
my $treeio = new Bio::TreeIO(-file   => "$in_tree",
						     -format => "newick"
						    );

my $total_length = 0;
my @branch_lengths = ();
my @depths_teleosts = ();
my @depths_tetrapods = ();
my @depths_sarc_fishes = ();

my $stat_tot = Statistics::Descriptive::Full->new();
my $stat_tet = Statistics::Descriptive::Full->new();
my $stat_sar = Statistics::Descriptive::Full->new();
my $stat_tel = Statistics::Descriptive::Full->new();

# print out heading

while( my $tree = $treeio->next_tree ) {


	# root the tree
    my $root = $tree->find_node( -id => $outgroup );
    $tree->reroot( $root );

	# get total branch length
	$total_length = $tree->total_branch_length;		
	
	for my $node ( $tree->get_nodes ) {

		# get tetrapod branches
		foreach my $t ( @tetrapods ) {
	
			my $tet = $tree->find_node(-id => $t);
			# $tet is a Bio::Tree::Node object
			# depth takes calculates sum of BL till the root
			my $tet_depth = $tet->depth;
			print "$t depth: $tet_depth\n";
			# scape if a particular species is not present in the tree
			next if (!defined $tet_depth);
			push (@depths_tetrapods, $tet_depth);
		}
		# get sarcopterygian fish branches
		foreach my $u ( @sarc_fishes ) {
	
			my $sar = $tree->find_node(-id => $u);
			# $tet is a Bio::Tree::Node object
			my $sar_depth = $sar->depth;
			print "$u depth: $sar_depth\n";
			next if (!defined $sar_depth);
			push (@depths_sarc_fishes, $sar_depth);
		}
		# get teleost branches
		foreach my $v ( @teleosts ) {
	
			my $tel = $tree->find_node(-id => $v);
			# $tet is a Bio::Tree::Node object
			my $tel_depth = $tel->depth;
			print "$v depth: $tel_depth\n";
			next if (!defined $tel_depth);
			push (@depths_teleosts, $tel_depth);
		}
		# get all nodes together (also previous ones!)

		my $length = $node->branch_length;
		# scape undefined branch lengths (root)
		next if (!defined $length);
		push (@branch_lengths, $length);
	}

	# calculate stats from @branch_lengths
	$stat_tot->add_data(@branch_lengths);
	my $mean = $stat_tot->mean();
	my $sd = $stat_tot->standard_deviation();
	
	# calcualte stats from @depths_tetrapods
	$stat_tet->add_data(@depths_tetrapods);
	my $tet_mean = $stat_tet->mean();
	my $tet_sd = $stat_tet->standard_deviation();
	my $tet_min = $stat_tet->min();
	my $tet_max = $stat_tet->max();

	# calcualte stats from @depths_sarc_fishes
	$stat_tet->add_data(@depths_sarc_fishes);
	my $sar_mean = $stat_tet->mean();
	my $sar_sd = $stat_tet->standard_deviation();
	my $sar_min = $stat_tet->min();
	my $sar_max = $stat_tet->max();

	# calcualte stats from @depths_teleosts
	$stat_tet->add_data(@depths_teleosts);
	my $tel_mean = $stat_tet->mean();
	my $tel_sd = $stat_tet->standard_deviation();
	my $tel_min = $stat_tet->min();
	my $tel_max = $stat_tet->max();

	# calculate mean differences among groups
	my $diff1 = abs ( $tet_mean - $sar_mean);
	my $diff2 = abs ( $tel_mean - $sar_mean );
	my $diff3 = abs ( $tet_mean - $tel_mean );
	
	
	
	# print out stuff
	print 	"File\tmean_BL\tsd_BL\ttotal_BL\t",
			"Tetrapods_mean\tTetrapods_sd\tTetrapods_min\tTetrapods_max\t",
			"Sarc_fish_mean\tSarc_fish_sd\tSarc_fish_min\tSarc_fish_max\t",
			"Teleosts_mean\tTeleosts_sd\tTeleosts_min\tTeleosts_max\t",
			"Diff_TetvsSarcFish\tDiff_TelvsSarcFish\tDiff_TetvsTel\n";
	print 	"$in_tree\t$mean\t$sd\t$total_length\t",
			"$tet_mean\t$tet_sd\t$tet_min\t$tet_max\t",
			"$sar_mean\t$sar_sd\t$sar_min\t$sar_max\t",
			"$tel_mean\t$tel_sd\t$tel_min\t$tel_max\t",
			"$diff1\t$diff2\t$diff3\n";
	
}
	
print STDERR "\nDone!\n\n";	

