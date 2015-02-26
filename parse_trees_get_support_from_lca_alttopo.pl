#!/usr/bin/perl

use strict;
use warnings;
use Bio::TreeIO;
use Data::Dumper;
use Statistics::Descriptive;

########################################################

##   Iker Irisarri, University of Konstanz. Feb 2015  ##

########################################################

# modif from parse_trees_test_alt_topologies.pl
# Script gets support values for the relevant node in 3 possible alternative topologies (rather than printing true/false message)
# 3 tests are predefined, while the 4th one must be chosen (comment out non-relevant options) to test for 3 possible topologies

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
my $monophyly_test_1 = 0;
my $monophyly_test_2 = 0;
my $monophyly_test_3 = 0;
my $monophyly_test_4 = 0;

# hard-coded arrays for dataset 251_amemiya_lf
#my $outgroup = "Chondricht";
#my @sarcs = qw (Latimeria_ Lepidosire Neoceratod Protopteru Anolis_car Canis_fami Dasypus_no Gallus_gal Homo_sapie Loxodonta_ Macropus_e Meleagris_ Monodelphi Mus_muscul Ornithorhy Rana_chine Taeniopygi Xenopus_tr);
#my @tets = qw (Anolis_car Canis_fami Dasypus_no Gallus_gal Homo_sapie Loxodonta_ Macropus_e Meleagris_ Monodelphi Mus_muscul Ornithorhy Rana_chine Taeniopygi Xenopus_tr);
#my @lfs = qw (Lepidosire Neoceratod Protopteru);
#my @tet_lf = qw (Lepidosire Neoceratod Protopteru Anolis_car Canis_fami Dasypus_no Gallus_gal Homo_sapie Loxodonta_ Macropus_e Meleagris_ Monodelphi Mus_muscul Ornithorhy Rana_chine Taeniopygi Xenopus_tr);
#my @tet_lat = qw (Latimeria_ Anolis_car Canis_fami Dasypus_no Gallus_gal Homo_sapie Loxodonta_ Macropus_e Meleagris_ Monodelphi Mus_muscul Ornithorhy Rana_chine Taeniopygi Xenopus_tr);
#my @lfs_lat = qw (Latimeria_ Lepidosire Neoceratod Protopteru);

# hard-coded arrays for dataset 1821_strict_3lf
my $outgroup = "outgroup";
my @sarcs = qw (Latimeria_chalumnae Lepidosiren_paradoxa Neoceratodus_forsteri Protopterus_annectens Anolis_carolinensis Canis_lupus_familiaris Dasypus_novemcinctus Gallus_gallus Homo_sapiens Loxodonta_africana Macropus_eugenii Meleagris_gallopavo Monodelphis_domestica Mus_musculus Ornithorhynchus_anatinus Pelodiscus_sinensis Taeniopygia_guttata Xenopus_Silurana_tropicalis);
my @tets = qw (Anolis_carolinensis Canis_lupus_familiaris Dasypus_novemcinctus Gallus_gallus Homo_sapiens Loxodonta_africana Macropus_eugenii Meleagris_gallopavo Monodelphis_domestica Mus_musculus Ornithorhynchus_anatinus Pelodiscus_sinensis Taeniopygia_guttata Xenopus_Silurana_tropicalis);
my @lfs = qw (Lepidosiren_paradoxa Neoceratodus_forsteri Protopterus_annectens);
my @tet_lf = qw (Lepidosiren_paradoxa Neoceratodus_forsteri Protopterus_annectens Anolis_carolinensis Canis_lupus_familiaris Dasypus_novemcinctus Gallus_gallus Homo_sapiens Loxodonta_africana Macropus_eugenii Meleagris_gallopavo Monodelphis_domestica Mus_musculus Ornithorhynchus_anatinus Pelodiscus_sinensis Taeniopygia_guttata Xenopus_Silurana_tropicalis);
my @tet_lat = qw (Latimeria_chalumnae Anolis_carolinensis Canis_lupus_familiaris Dasypus_novemcinctus Gallus_gallus Homo_sapiens Loxodonta_africana Macropus_eugenii Meleagris_gallopavo Monodelphis_domestica Mus_musculus Ornithorhynchus_anatinus Pelodiscus_sinensis Taeniopygia_guttata Xenopus_Silurana_tropicalis);
my @lfs_lat = qw (Latimeria_chalumnae Lepidosiren_paradoxa Neoceratodus_forsteri Protopterus_annectens);

# select 4th monophyly test that will discriminate between the 3 possible hypotheses
my @custom_test = @tet_lf;
#my @custom_test = @tet_lat;
#my @custom_test = @lfs_lat;

=pod

# get taxa from $clade file
open (IN, "<", $clade) or die "Cannot open file $clade: $!\n";
open (OUT, ">", $outfile) or die "Cannot create output file: $!\n";

while ( my $line = <IN> ) {
    chomp $line;
    next if ( $line =~ /^#.+/ );
    push (@clade, $line);
}

close(IN);

=cut


# read in trees
my $treeio = new Bio::TreeIO(-file   => "$in_tree",
						     -format => "newick");

my $stat = Statistics::Descriptive::Full->new();

# print out trees that pass the four monophyly tests
# IMPORTANT NOTE: all trees are printed together in the outfile (no newlines are added)
# 				  and there is an issue with how the outgroup is written (need to check) 
#my $out = new Bio::TreeIO(-file   => ">outfile.tre",
#                          -format => "newick");



# print file name to stdout
print "File: ", $in_tree, "\n";

while( my $tree = $treeio->next_tree ) {

	# initialize test results for each tree
	$monophyly_test = 0;
	$monophyly_test_1 = 0;
	$monophyly_test_2 = 0;
	$monophyly_test_3 = 0;
	$monophyly_test_4 = 0;

    # track tree number if multiple trees are stored in a single file
    $tree_num++;
    if ( $tree_num > 1) {
	print "Tree number $tree_num\n";
    }
    
    # root the tree with outgroup
	my $root = $tree->find_node( -id => $outgroup );
	$tree->reroot( $root );


	# 1st monophyly test
	# send $tree and @sarcs to subroutine
	my $monophyly_test_1 = monophyly_test($tree, \@sarcs);
	
	# If 1st test successful, do 2nd monophyly test
	# send $tree and taxa @array ref to subroutine
	if ( $monophyly_test_1 == 1 ) {
		
		#print "\t1st monophyly test succeeded (sarcs)\n";
		my $monophyly_test_2 = monophyly_test($tree, \@tets);
		
		# If 2nd test successful, do 3rd monophyly test
		# send $tree and taxa @array ref to subroutine
		if ( $monophyly_test_2 == 1 ) {
	
			#print "\t2nd monophyly test succeeded (tets)\n";
			my $monophyly_test_3 = monophyly_test($tree, \@lfs);
			
			# If 3rd test successful, do 4th monophyly test
			# send $tree and taxa @array ref to subroutine
			if ( $monophyly_test_3 == 1 ) {
	
				#print "\t3rd monophyly test succeeded (lfs)\n";
				my $monophyly_test_4 = monophyly_test($tree, \@custom_test);
			
				if ( $monophyly_test_4 == 1 ) {
	
					#print "\t4th monophyly test succeeded (custom)\n";
					print "\tall 4 monophyly tests succeeded!!\n";
					#$out->write_tree($tree);
			
					## GET SUPPORT VALUE FROM THE LCA OF THE TAXA IN @custom_test
					# get Bio::Tree::Node object for the taxa
					my @custom_test_obj = ();
					foreach my $i (@custom_test) {
						my $node = $tree->find_node(-id => $i);
						# $node is a Bio::Tree::Node object
						push (@custom_test_obj, $node);
					}

					# bootstrap value for the MRCA for the taxa in @custom_test
					# function returns a NodeI object
					my $lca = $tree->get_lca(-nodes => \@custom_test_obj);
					#print $lca->id, "\n";
					my $lca_support = $lca->id;
					
					print "\tsupport for relevant node: $lca_support\n";
					
				}
				else {
					print "\t4th monophyly test failed (custom)\n";
				}
			}
			else {
				print "\t3rd monophyly test failed (lungfishes)\n";
			}
		}	
		else {
			print "\t2nd monophyly test failed (tetrapods)\n";
		}
	}
	else {
		print "\t1st monophyly test failed (sarcopterygii)\n";
	}
}
	
	

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
		foreach my $i (@clade2) {
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
