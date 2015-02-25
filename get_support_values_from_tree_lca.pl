#!/usr/bin/perl

use strict;
use warnings;

use Bio::TreeIO;
use Data::Dumper;
use Statistics::Descriptive;

#########################################################################

# Iker Irisarri, University of Konstanz. Feb 2015
# get_suppot_values_from_tree.pl
# requires (1) one tree; (2) file with taxa for which we want to get the support of their MRCA
# (3) outgroup to roo the tree
# Prints mean support values and the support for the MRCA of the reference taxa
# Rooting is required to make sure the MRCA node is the correct one

# It gives some warnings when nodes do not have support values. To avoid, filter the output:
# get_support_values_from_tree.pl 251_trees_with_shlike_support.tre clade | grep LCA | sort | cut -f2 | uniq -c

#########################################################################

my $usage = "get_support_values_from_tree_lca.pl tree lca_file outgroup> stdout\n";
my $in_tree = $ARGV[0] or die $usage;
my $clade = $ARGV[1] or die $usage;
my $outgroup = $ARGV[2] or die $usage;

my @clade;
	
open (IN, "<", $clade) or die "Cannot open file $clade: $!\n";

while ( my $line = <IN> ) {
   	chomp $line;
   	next if ( $line =~ /^#.+/ );
   	push (@clade, $line);
}

close(IN);



# read in trees
my $treeio = new Bio::TreeIO(-file   => "$in_tree",
							-format => "newick");

my $stat = Statistics::Descriptive::Full->new();


while( my $tree = $treeio->next_tree ) {
	
	    my $root = $tree->find_node( -id => $outgroup );
	    $tree->reroot( $root );
	
	
	# initialize matrices & counts for every tree (if multiple trees in one file)
    # reasign taxa for monophyly test to @clade2
	
	my @bootstraps;
	my @clade2;
	my @clade2_obj;

    # store taxa into a hash where keys are taxa names
    my %leaves = %{ &get_all_leaves($tree->get_root_node) };

    # remove taxa from @clade2 if not present in this particular tree
   foreach my $taxa (@clade) {
		if ( exists $leaves{$taxa} ) {
		# push into a new array
	    	push (@clade2, $taxa);
		}
	}

	#print Dumper \@clade2;

## GET SUPPORT VALUE FROM THE LCA OF THE TAXA IN @clade2
	# get Bio::Tree::Node object for the taxa
	foreach my $i (@clade2) {
	    my $node = $tree->find_node(-id => $i);
	    # $node is a Bio::Tree::Node object
	    push (@clade2_obj, $node);
	}

	# bootstrap value for the MRCA for the taxa in @clade2
	# function returns a NodeI object
	my $lca = $tree->get_lca(-nodes => \@clade2_obj);
	#print $lca->id, "\n";
	my $lca_support = $lca->id;

## GET SUPPOR VALUE FOR ALL NODES TO CALCULATE MEAN		
	# coming back to the whole tree...
	# move the bootstrap values over from the internal node Ids
	$tree->move_id_to_bootstrap;

	for my $node ( $tree->get_nodes ) {
	
		my $bootstrap = $node->bootstrap;

		# filter out undefined bootstrap values (nodes that do not have bootstrap values)
		# This dumps a warning but it is not problematic
		# Use of uninitialized value $bootstrap in string ne at /Applications/Bioinformatics/scripts_github/lf_phylo/get_support_values_from_tree.pl line 111, <GEN0> line 1.
		my $length = length $bootstrap;

		# exclude nodes without bootstrap values && discard outgroup name
		# (It is treated as support value, probably because $tree->move_id_to_bootstrap
		if ($bootstrap ne "" && $bootstrap ne $outgroup) {
			push (@bootstraps, $bootstrap);
		}
	}

	# calculate mean for all nodes
	$stat->add_data(@bootstraps);
	my $mean = $stat->mean();

	# PRINT INFO
	print "mean support values:\t$mean\n";
	print "support value for LCA:\t$lca_support\n";
	#print join ('--', @bootstraps);
		
}
	

## SUBROUTINES ##

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

__END__

example of tree used as input

((Takifugu_r,Oreochromi)54,((Latimeria_,Chondricht)31,((Lepidosire,Protopteru)99,(Xenopus_tr,(Anolis_car,((Taeniopygi,(Meleagris_,Gallus_gal)89)99,((Dasypus_no,((Loxodonta_,Mus_muscul)82,(Canis_fami,Homo_sapie)83)71)100,(Ornithorhy,(Macropus_e,Monodelphi)84)79)98)52)93)93)91)98,Danio_reri);

