#!/usr/bin/perl

use strict;
use warnings;
use Bio::TreeIO;
use Data::Dumper;

########################################################

##   Iker Irisarri, University of Konstanz. Jun 2015  ##

########################################################

# modif from parse_trees_test_alt_topologies.pl
# Script will test for monophyly for 4 clades and print a messages for each of the tests
# A total of 6 monophyly tests are performed

# Intended to be run in a for loop (bash) to test multiple trees (in individual files)
# possibility to run also a single file with multiple trees (check options bellow for tracking trees)
# Prints file name and checks if all the taxa in $clade are really present in the tree and remove missing taxa accordingly
# This allows running the script for a bunch of trees that might have a slightly different taxon sampling

# is_monophyletic in TreeIO module is broken for some reason (J. Stajich G. Jordan) so I implemented my own topology test
# see parse_tree_check_monophyly.pl for more details

#########################################################################

my $usage = "parse_trees_test_alt_topologies.pl tree > stdout\n";
my $in_tree = $ARGV[0] or die $usage;

my $tree_num = 0;
# boolean variables to control if monophyly tests succeed ("=0": monophyly is FALSE; "=1" monophyly is TRUE)
my $monophyly_test = 0;
my $monophyly_sarcs = 0;
my $monophyly_tets = 0;
my $monophyly_lfs = 0;
my $monophyly_tet_lf = 0;
my $monophyly_tet_lat = 0;
my $monophyly_lfs_lat = 0;

my $sarcopterygii = 0;
my $tetrapoda = 0;
my $dipnoi = 0,
my $sarcs_and_tets = 0;
my $sarcs_and_lfs = 0;
my $tets_and_lfs = 0;
my $T1 = 0;
my $T2 = 0;
my $T3 = 0;
my $unresolved = 0;

# hard-coded arrays for dataset 251_amemiya_lf
#my $outgroup = "Chondricht";
#my @sarcs = qw (Latimeria_ Lepidosire Neoceratod Protopteru Anolis_car Canis_fami Dasypus_no Gallus_gal Homo_sapie Loxodonta_ Macropus_e Meleagris_ Monodelphi Mus_muscul Ornithorhy Rana_chine Taeniopygi Xenopus_tr);
#my @tets = qw (Anolis_car Canis_fami Dasypus_no Gallus_gal Homo_sapie Loxodonta_ Macropus_e Meleagris_ Monodelphi Mus_muscul Ornithorhy Rana_chine Taeniopygi Xenopus_tr);
#my @lfs = qw (Lepidosire Neoceratod Protopteru);
#my @tet_lf = qw (Lepidosire Neoceratod Protopteru Anolis_car Canis_fami Dasypus_no Gallus_gal Homo_sapie Loxodonta_ Macropus_e Meleagris_ Monodelphi Mus_muscul Ornithorhy Rana_chine Taeniopygi Xenopus_tr);
#my @tet_lat = qw (Latimeria_ Anolis_car Canis_fami Dasypus_no Gallus_gal Homo_sapie Loxodonta_ Macropus_e Meleagris_ Monodelphi Mus_muscul Ornithorhy Rana_chine Taeniopygi Xenopus_tr);
#my @lfs_lat = qw (Latimeria_ Lepidosire Neoceratod Protopteru);

# hard-coded arrays for dataset 1821_smx_3lf
my $outgroup = "Callorhinchus_milii";
#my $outgroup = "outgroup";
my @sarcs = qw (Latimeria_chalumnae Lepidosiren_paradoxa Neoceratodus_forsteri Protopterus_annectens Gallus_gallus Anolis_carolinensis Macropus_eugenii Taeniopygia_guttata Homo_sapiens Xenopus_Silurana_tropicalis Monodelphis_domestica Ornithorhynchus_anatinus Loxodonta_africana Pelodiscus_sinensis Dasypus_novemcinctus Canis_lupus_familiaris Meleagris_gallopavo Mus_musculus);
my @tets = qw (Gallus_gallus Anolis_carolinensis Macropus_eugenii Taeniopygia_guttata Homo_sapiens Xenopus_Silurana_tropicalis Monodelphis_domestica Ornithorhynchus_anatinus Loxodonta_africana Pelodiscus_sinensis Dasypus_novemcinctus Canis_lupus_familiaris Meleagris_gallopavo Mus_musculus);
my @lfs = qw (Lepidosiren_paradoxa Neoceratodus_forsteri Protopterus_annectens);
my @tet_lf = qw (Lepidosiren_paradoxa Neoceratodus_forsteri Protopterus_annectens Gallus_gallus Anolis_carolinensis Macropus_eugenii Taeniopygia_guttata Homo_sapiens Xenopus_Silurana_tropicalis Monodelphis_domestica Ornithorhynchus_anatinus Loxodonta_africana Pelodiscus_sinensis Dasypus_novemcinctus Canis_lupus_familiaris Meleagris_gallopavo Mus_musculus);
my @tet_lat = qw (Latimeria_chalumnae Gallus_gallus Anolis_carolinensis Macropus_eugenii Taeniopygia_guttata Homo_sapiens Xenopus_Silurana_tropicalis Monodelphis_domestica Ornithorhynchus_anatinus Loxodonta_africana Pelodiscus_sinensis Dasypus_novemcinctus Canis_lupus_familiaris Meleagris_gallopavo Mus_musculus);
my @lfs_lat = qw (Latimeria_chalumnae Lepidosiren_paradoxa Neoceratodus_forsteri Protopterus_annectens);

# read in trees
my $treeio = new Bio::TreeIO(-file   => "$in_tree",
						     -format => "newick");

# print out trees that pass the four monophyly tests
# IMPORTANT NOTE: all trees are printed together in the outfile (no newlines are added)
# 				  and there is an issue with how the outgroup is written (need to check) 
#my $out = new Bio::TreeIO(-file   => ">outfile.tre",
#                          -format => "newick");

# print file name to stdout
print "# monophyly\tsarcs\ttets\tlfs\ttet_lf\ttet_lat\tlf_lat\n";
print "File: ", $in_tree, "\t";

while( my $tree = $treeio->next_tree ) {

	# initialize test results for each tree
	$monophyly_test = 0;
	$monophyly_sarcs = 0;
	$monophyly_tets = 0;
	$monophyly_lfs = 0;
	$monophyly_tet_lf = 0;
	$monophyly_tet_lat = 0;
	$monophyly_lfs_lat = 0;

    # track tree number if multiple trees are stored in a single file
    $tree_num++;
    if ( $tree_num > 1) {
	print "Tree # $tree_num\t";
    }
    
    # root the tree with outgroup
	my $root = $tree->find_node( -id => $outgroup );
	$tree->reroot( $root );

	# 1st monophyly test
	# send $tree and @sarcs to subroutine
	$monophyly_sarcs = monophyly_test($tree, \@sarcs);

	print "$monophyly_sarcs\t";
	
	$monophyly_tets = monophyly_test($tree, \@tets);
	
	print "$monophyly_tets\t";
	
	$monophyly_lfs = monophyly_test($tree, \@lfs);
	
	print "$monophyly_lfs\t";

	$monophyly_tet_lf = monophyly_test($tree, \@tet_lf);
	
	print "$monophyly_tet_lf\t";

	$monophyly_tet_lat = monophyly_test($tree, \@tet_lat);
	
	print "$monophyly_tet_lat\t";

	$monophyly_lfs_lat = monophyly_test($tree, \@lfs_lat);
	
	print "$monophyly_lfs_lat\n";
		
	if ( $monophyly_sarcs == 1 && $monophyly_tets == 1 && $monophyly_lfs == 1 && $monophyly_tet_lf ) {
		
		$T1++;
	}
	if ( $monophyly_sarcs == 1 && $monophyly_tets == 1 && $monophyly_lfs == 1 && $monophyly_lfs_lat ) {
		
		$T2++;
	}
	if ( $monophyly_sarcs == 1 && $monophyly_tets == 1 && $monophyly_lfs == 1 && $monophyly_tet_lat ) {
		
		$T3++;
	}
	if ( $monophyly_sarcs == 1 && $monophyly_tets == 1 && $monophyly_lfs == 0 ) {
		
		$sarcs_and_tets++;
	}
	if ( $monophyly_sarcs == 1 && $monophyly_tets == 0 && $monophyly_lfs == 1 ) {
		
		$sarcs_and_lfs++;
		#print "$tree_num\n";
	}
	if ( $monophyly_sarcs == 0 && $monophyly_tets == 1 && $monophyly_lfs == 1 ) {
		
		$tets_and_lfs++;
	}
	if ( $monophyly_sarcs == 1 && $monophyly_tets == 0 && $monophyly_lfs == 0 ) {
		
		$sarcopterygii++;
	}
	if ( $monophyly_sarcs == 0 && $monophyly_tets == 1 && $monophyly_lfs == 0 ) {
		
		$tetrapoda++;
	}
	if ( $monophyly_sarcs == 0 && $monophyly_tets == 0 && $monophyly_lfs == 1 ) {
		
		$dipnoi++;
	}
	if ( $monophyly_sarcs == 0 && $monophyly_tets == 0 && $monophyly_lfs == 0 ) {
		
		$unresolved++;
	}
	#else {
		
	#	print "error: non-tracked topology. Check $tree_num\n";
	#}
}

# Print summary

print "\nSummary of results:\n";

print "\nFully resolved topologies:\n";
print "$T1\t(coelacanth,(lungfishes,tetrapods))\n";
print "$T2\t((coelacanth,lungfishes),tetrapods)\n";
print "$T3\t(lungfishes,(coelacanth,tetrapods))\n";

print "\nRemaining topologies:\n";
print "$sarcs_and_tets\tsarcopterygians and tetrapods\n";
print "$sarcs_and_lfs\tsarcopterygians and lungfishes\n";
print "$tets_and_lfs\ttetrapods and lungfishes\n";
print "$sarcopterygii\tonly sarcopterygians\n";
print "$tetrapoda\tonly tetrapods\n";
print "$dipnoi\tonly lungfishes\n";
print "$unresolved\tnon-monophyletic sarcs, tets and lungfsihes\n\n";

## SUBROUTINES ##


sub monophyly_test {

	# passed to the subroutine: tree and array with taxa for monophyly test
	my $current_tree = shift;
	my $clade_ref = shift;
	# de-reference array (arrays cannot be passed to subroutines)
	my @clade = @{ $clade_ref };

	# initialized
	my $tax_num = 0;
	my @clade2;				# contains taxa
	my @test_taxa;			# contains Bio::Tree::Node objects
	my @lca_descend;		# contains taxa
	my @lca_desc;			# contains Bio::Tree::Node objects
	my @sorted_clade;		# contains taxa
	my @sorted_lca_descend;	# contains taxa

	# initialize matrices & counts for every tree (if multiple trees in one file)
	#@test_taxa = ();
	#$tax_num = 0;

	# store taxa into a hash where keys are taxa names
	my %leaves = %{ &get_all_leaves($current_tree->get_root_node) };

	# store taxa from @clade that is present in the current tree (remove taxa not present in the tree)
	foreach my $taxa (@clade) {	
		if ( exists $leaves{$taxa} ) {
	    	push (@clade2, $taxa);
		}
	}

=pod
	# (re)root tree 
	foreach my $key (keys %leaves) {
		if (!exists $leaves{$outgroup} ) {
			print "\tError: cannot find $outgroup in the tree!\n";
			# exit foreach loop after error
			exit;
		}
		if ( $key eq $outgroup ) {
			my $root = $current_tree->find_node( -id => $key );
			$current_tree->reroot( $root );
		}
	}    

=cut

	# 1st, check that we have at least 2 taxa in @clade2
	$tax_num = scalar @clade2;
	if ( $tax_num < 2 ) {
		print "\tAt least two taxa are required for testing monophyly and only $tax_num is present\n\n";
	}
	else {

		# 2nd, need to find the nodes for taxa in the monophyly test & outgroup!
		# find node for taxa 
		foreach my $i ( @clade2 ) {
			my $node = $current_tree->find_node(-id => $i);
			# $node is a Bio::Tree::Node object
			push (@test_taxa, $node);
		}


	### Custom monophyly test ###

		# (1) get the mrca of species for which we want to test monophyly
		# function returns a NodeI object
		my $lca = $current_tree->get_lca(-nodes => \@test_taxa);

		# (2) get all descendents from that mrca
		my @lca_desc_obj = $lca->get_all_Descendents();
	
		# (3) store all descendent leaves into an array
		foreach my $desc ( @lca_desc_obj ) {
			if ( $desc->is_Leaf ) {
			# get leaf id (actual taxa name) and store it
			my $leaf = $desc->id();
			push ( @lca_descend, $leaf);
			}
		}

		# (4) sort both arrays
		@sorted_clade = sort (@clade);
		@sorted_lca_descend = sort (@lca_descend);
	
		# (5) compare sorted arrays
		if (@sorted_clade == @sorted_lca_descend){
		
			# if monophyletic, $monophyly_test is TRUE
			$monophyly_test = 1;
			#print "\tMonophyly of the following taxa:\n\t";
			#print join (" ", @sorted_clade), "\n";
		}
		else {
			
			# otherwise, $monophyly_test is FALSE
			$monophyly_test = 0;
			#print "\tTaxa are not monophyletic :-(\n\n";
		}

	### End of monophyly test ###

	}
	
	# rapid check for monophyly
	#if ( $monophyly_test == 1 ) {
	#	print "monophyly\n";
	#}
	#else { 
	#	print "non-monophyly\n";
	#}

	# subroutine returns true if taxa in array are monophyletic
	return $monophyly_test;
#	print "\tresult: $monophyly_test\n";
}


sub get_all_leaves {
    #Esta rutina devuelve todas las "hojas" descendientes de un nodo
    my $initial_node = $_[0];
    my %nodes;
    if( $initial_node->is_Leaf ) {
	$nodes{ $initial_node->id } = 1;
	return \%nodes;
    }
    foreach my $node ( $initial_node->get_all_Descendents() ) {
	if( $node->is_Leaf ) { 
	    # for example use below
	    $nodes{$node->id} = 1;
	}
    }
    return \%nodes;
}
