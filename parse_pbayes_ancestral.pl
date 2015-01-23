#!/usr/bin/perl

use strict;
use warnings;
use Bio::TreeIO;
use Data::Dumper;

# parse_pbayes_ancestral.pl
# Iker Irisarri, University of Konstanz, Jan 2014
# It parses output of ancestral (phylobayes package) that generates a complete history of substitution along the phylogeny, and outputs replacements that ocurr in >=90% of the points/replicates for a number of nodes of interest
# It seems phylobayes already removes the initial cycles (20%), so all points that are written to *.sub are post-burnin

# ancestral should be run with a fixed tree topology, which it is fixed to tet+lf, but could be modified
# IMPORTANT!! We can only know the TYPE of substitutions, NOT the DIRECTION. This is because of the way ancestral outputs trees (first and last elements of the string are not the ancestral and derived states, respectively. It depends on how the tree is written by ancestral.

my $usage = "parse_pbayes_ancestral.pl infile.sub > stdout\n";
my $in_tree = $ARGV[0] or die $usage;


open (IN, "<", $in_tree) or die "Cannot open file $in_tree: $!\n";

my $point = 0;   # generations
my $rep = 0;     # replicates
my $aln_pos = 0; # position in alignment (after invariants are removed)
my $tree;

my %data;
my %tet;
my %lfs;
my %lfs_tet;
my %sarco;

# parse-in file from pbayes ancestral
while ( my $line = <IN> ) {
    chomp $line;
    # get points
    if ( $line =~ /^point/ ) {
	$point++;
    }
    # get replicates
    elsif ( $line =~ /^rep/ ) {
	$rep++;
    }
    # get actual info for each aln position
    elsif ( $line =~ /^\w+\t.+/ ) {

	($aln_pos, $tree) = split ("\t", $line);


## NEED TO CHOOSE ONE TOPOLOGY AND COMMENT OUT THE OTHER TOPOLOGY ##

	
	# match synapomorphic changes for tets in tree lf+tet
	$tree =~ /\)([\w\d\.:-]+),\(\(Protopteru_/;

	# example of information from node:
	# M:0.0078287:M:0.0249226:L:0.130534:M

	# submit to subroutine and de-reference array
	my $replacements_tet = get_replacements($1); 

	my @replacements_tet = @{ $replacements_tet}; 

	# save results in new matrix %tet
	$tet{$aln_pos}{$point}{$rep} = [@replacements_tet];


	# match synapomorphic changes for lfs (1) and lf+tet (2) in tree lf+tet
	$tree =~/Neoceratod_[\w\d\.:-]+\)([\w\d\.:-]+)\)([\w\d\.:-]+),Latimeria_/;

	my $replacements_lfs = get_replacements($1);
	my $replacements_lfs_tet = get_replacements($2);

        my @replacements_lfs = @{ $replacements_lfs };
	my @replacements_lfs_tet = @{ $replacements_lfs_tet };

        $lfs{$aln_pos}{$point}{$rep} = [@replacements_lfs];
        $lfs_tet{$aln_pos}{$point}{$rep} = [@replacements_lfs_tet];


	# match synapomorphic changes for sarco in tree lf+tet
	$tree =~/Latimeria__[\w\d\.:-]+\)([\w\d\.:-]+),\(\(Oreochromi/;

	my $replacements_sarco = get_replacements($1);

        my @replacements_sarco = @{ $replacements_sarco };

	$sarco{$aln_pos}{$point}{$rep} = [@replacements_lfs];


    }
}

# send hashes to subroutine to analyze and print

print "\ntetrapods\n";
count_replacements_and_print(\%tet);

print "\nlungfishes\n";
count_replacements_and_print(\%lfs);

print "\nlungfish+tetrapods\n";
count_replacements_and_print(\%lfs_tet);

print "\nsarcopterygians\n";
count_replacements_and_print(\%sarco);



print "\ndone!\n\n";

## SUBROUTINES ##

sub get_replacements {

    my $node = shift;
    my @node_infos = split (":", $node);
    my @replacements = ();

    # if there is a replacmeent in that node the list should have >3 elements
    # and element 0 (ancestral) and 5 (derived) should be different

    if ( scalar @node_infos > 3 ) {
        if ( $node_infos[0] ne $node_infos[4] ) {

            # store amino acid changes (pair elements in list) in new array
            for ( my $i=0; $i<scalar(@node_infos)+1; $i+=2 ) {

                push (@replacements, $node_infos[$i]);
            }

        }
    }

    # if there is no replacement, store 'none'
    else {

        @replacements = 'none';

    }

    # returns reference to array
    return \@replacements;

}



sub count_replacements_and_print {

    my $hash_ref = shift;
    my $point_count;
    my $change_count;
    my %synapom;

    # de-reference hash that contains info on replacements for that node
    my %clade_results = %{ $hash_ref };


    foreach my $aln ( keys %clade_results ) {

	# reinitialize counts and @t_synapom for each new aln position
	$point_count = 0;
	$change_count = 0;
	%synapom = ();

	foreach my $point ( keys %{ $clade_results{$aln} } ) {
	    
	    #here it would be possible to remove initial cycles (burnin) by setting something like this and define burnin as an argument    
	    #if ( $point > $burnin ) {
	    
	    foreach my $rep (keys %{ ${ $clade_results{$aln} }{$point} } ) {

		# count empty arrays, withouth changes ('none')
		$point_count++;

		my @value_array = @{ ${ ${ $clade_results{$aln} }{$point} }{$rep} };

		# process non-empty arrays (that contain changes)
		if ( $value_array[0] ne 'none' ) {

		    # take first and last states (i.e. ignore substitution history)
		    # VERY IMPORTANT!! The way ancestral is programmed, the first element is no necessary the ancestral state and the last element the derived state. It depends on how the tree is written. So it shows us the replacement type not direction!!
		    #e.g. (Protopteru_A:0.0403415:A,Lepidosire_A:0.0337664:A)A:0.113037:A,Neoceratod_A:0.0607286:A)#A:0.0668141:A:0.0655169:V#)
		    # for node ## V is ancestral (present in all other species) and A is derived (in lungfishes)

		    my $first = shift @value_array;
		    my $last = pop @value_array;

		    my $synapomorphy = $first . '<->' . $last;

		    if ( !exists $synapom{$synapomorphy} ) {
			$change_count = 1;
		    }
		    else {

		    $change_count++;
		    }

		    # non-emtpy arrays into a hash
		    $synapom{$synapomorphy} = $change_count;

		}
	    }
	}

	foreach my $syn ( keys %synapom ) {

	    my $percent = ($synapom{$syn}/$point_count)*100;
	
	    # print out replacement if present in at least 90% of the points/replicates
	    if ( $percent > 90 ) {

		print "\n\taln pos: $aln\t";
		print $syn, "\tin $percent % of points/replicates";

	    }
	}

    }

    print "\n";

}


__END__

# example of file from sub

point 1
rep 1
1	(((((((((((Mus_muscul_V:0.1453:V,Homo_sapie_V:0.0413422:V)V:0.000799158:V,Canis_fami_V:0.0488305:V)V:0.0217549:V,(Dasypus_no_V:0.0533447:V,Loxodonta__V:0.0368459:V)V:0.00154624:V)V:0.0280085:V,(Monodelphi_V:0.0286486:V,Macropus_e_V:0.0226434:V)V:0.0472848:V)V:0.0230245:V,Ornithorhy_V:0.0344183:V)V:0.0642463:V,(((Meleagris__V:0.0131901:V,Gallus_gal_V:0.00870232:V)V:0.0269581:V,Taeniopygi_V:0.0462125:V)V:0.0366126:V,Anolis_car_V:0.181893:V)V:0.0117471:V)V:0.0495589:V,(Xenopus_tr_V:0.0730028:V,Rana_chine_V:0.188929:V)V:0.0493053:V:0.0318468:A:0.084293:V)V:0.0736126:V,((Protopteru_A:0.0403415:A,Lepidosire_A:0.0337664:A)A:0.113037:A,Neoceratod_A:0.0607286:A)A:0.0668141:A:0.0655169:V)V:0.00686853:V,Latimeria__V:0.235074:V)V:0.0034714:V,((Oreochromi_V:0.0303549:V,Takifugu_r_V:0.141675:V)V:0.10266:V,Danio_reri_I:0.0583341:I:0.0331744:V)V:0.449703:V)V:0:V,((Leucoraja__V:0.199734:V,Scyliorhin_V:0.158381:V)V:0.180696:V,Callorhinc_V:0.211086:V)V:0.230053:V)V;
2	(((((((((((Mus_muscul_M:0.1453:M,Homo_sapie_M:0.0413422:M)M:0.000799158:M,Canis_fami_M:0.0488305:M)M:0.0217549:M,(Dasypus_no_M:0.0533447:M,Loxodonta__M:0.0368459:M)M:0.00154624:M)M:0.0280085:M,(Monodelphi_M:0.0286486:M,Macropus_e_M:0.0226434:M)M:0.0472848:M)M:0.0230245:M,Ornithorhy_M:0.0344183:M)M:0.0642463:M,(((Meleagris__M:0.0131901:M,Gallus_gal_M:0.00870232:M)M:0.0269581:M,Taeniopygi_M:0.0462125:M)M:0.0366126:M,Anolis_car_M:0.181893:M)M:0.0117471:M)M:0.0495589:M,(Xenopus_tr_M:0.0730028:M,Rana_chine_M:0.188929:M)M:0.165445:M)M:0.0736126:M,((Protopteru_M:0.0403415:M,Lepidosire_M:0.0337664:M)M:0.113037:M,Neoceratod_M:0.0607286:M)M:0.132331:M)M:0.00686853:M,Latimeria__M:0.235074:M)M:0.0034714:M,((Oreochromi_M:0.0303549:M,Takifugu_r_M:0.141675:M)M:0.10266:M,Danio_reri_M:0.0915085:M)M:0.449703:M)M:0:M,((Leucoraja__V:0.199734:V,Scyliorhin_V:0.158381:V)V:0.180696:V,Callorhinc_V:0.211086:V)V:0.192052:V:0.0380006:M)M;

etc.
