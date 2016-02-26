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
my %teleosts = ("Takifugu_r" => "1",
				"Danio_reri" => "1",
				"Oreochromi" => "1");

my %sarc_fishes = (	"Lepidosire" => "1",
					"Neoceratod" => "1",
					"Protopteru" => "1",
					"Latimeria_" => "1");

my %tetrapods = (	"Canis_fami" => "1",
					"Dasypus_no" => "1",
					"Gallus_gal" => "1",
					"Homo_sapie" => "1",
					"Loxodonta_" => "1",
					"Macropus_e" => "1",
					"Meleagris_" => "1",
					"Monodelphi" => "1",
					"Mus_muscul" => "1",
					"Ornithorhy" => "1",
					"Rana_chine" => "1",
					"Taeniopygi" => "1",
					"Xenopus_tr" => "1");
				
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

       	# get node leafs
       	if ( $node->is_Leaf ) {
			
		my $species = $node->id;
		my $species_node = $tree->find_node(-id => $species);
		my $depth = $species_node->depth;
						
		# scape undefined depths (often related to rooting)
		next if (!defined $depth);
		
		# get depths for tetrapods, teleosts and sarcopterygian fishes
		if ( exists $tetrapods{$species} ) {

			push (@depths_tetrapods, $depth);
		}
		if ( exists $sarc_fishes{$species} ) {

			push (@depths_sarc_fishes, $depth);
		}
		if ( exists $teleosts{$species} ) {

			push (@depths_teleosts, $depth);
		}
		
	# get branch length (for all leaves and internal branches)
	my $length = $node->branch_length;
	# scape if branch length is undefined (often related to rooting)
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

