#!/usr/bin/perl

use strict;
use warnings;
#use diagnostics;

use Data::Dumper;
use Bio::SeqIO;


# Iker Irisarri, Mar 2016. University of Konstanz.
# site_pattern_support_for_clades.pl
# Identifies alignment positions supporting a particular clade (ingroup), divided in 3 classes:
# BINARY: both ingroup and outgroup are homogeneous (have the same amino acid, different between clades)
# ASYMMETRICAL: ingroup is homogeneous and its amino acid is not allowed in outgroup
# NOISY: ingroup is homogeneous and its amino acid is allowed in outgroup up to 25%

# Based on the method presented in: Wägele J-W, Rödding F. 1998. A priori estimation of 
# phylogenetic information conserved in aligned sequences. Mol Phylogenet Evol 9:358-365.

# Requires prior trimming of data set to remove invariant positions and any position without gaps. E.g.:
# sed 's/X/-/g' infile.fas (BE CAREFUL WITH TAXON NAMES CONTAINING X!)
# java -jar /Applications/Phylogeny/BMGE-1.12/BMGE.jar -t AA -g 0 -h 1E-5:1 -i 251_amemiya_lf.fas -of 251_amemiya_lf_nogap_noinv.fas



my $usage = "site_pattern_support_for_clades.pl in.fasta > STDOUT\n";
my $in_fasta = $ARGV[0] or die $usage;

my $proportion_noise = "0.25";

# 251 data set
# Sarcopterygii
#my @ingroup = qw( Latimeria_ Lepidosire Protopteru Neoceratod Xenopus_tr Rana_chine Anolis_car Loxodonta_ Monodelphi Dasypus_no Macropus_e Homo_sapie Gallus_gal Mus_muscul Taeniopygi Meleagris_ Ornithorhy Canis_fami);
#my @outgroup = qw( Callorhinc Scyliorhin Leucoraja_ Danio_reri Takifugu_r Oreochromi);
# Tetrapoda
#my @ingroup = qw( Xenopus_tr Rana_chine Anolis_car Loxodonta_ Monodelphi Dasypus_no Macropus_e Homo_sapie Gallus_gal Mus_muscul Taeniopygi Meleagris_ Ornithorhy Canis_fami);
#my @outgroup = qw( Callorhinc Scyliorhin Leucoraja_ Danio_reri Takifugu_r Oreochromi Latimeria_ Lepidosire Protopteru Neoceratod);
# T1
#my @ingroup = qw( Lepidosire Protopteru Neoceratod Xenopus_tr Rana_chine Anolis_car Loxodonta_ Monodelphi Dasypus_no Macropus_e Homo_sapie Gallus_gal Mus_muscul Taeniopygi Meleagris_ Ornithorhy Canis_fami);
#my @outgroup = qw( Latimeria_ Callorhinc Scyliorhin Leucoraja_ Danio_reri Takifugu_r Oreochromi);
# T2
#my @ingroup = qw( Latimeria_ Xenopus_tr Rana_chine Anolis_car Loxodonta_ Monodelphi Dasypus_no Macropus_e Homo_sapie Gallus_gal Mus_muscul Taeniopygi Meleagris_ Ornithorhy Canis_fami);
#my @outgroup = qw( Callorhinc Scyliorhin Leucoraja_ Danio_reri Takifugu_r Oreochromi Lepidosire Protopteru Neoceratod);
# T3
#my @ingroup = qw( Lepidosire Protopteru Neoceratod Latimeria_);
#my @outgroup = qw( Danio_reri Takifugu_r Scyliorhin Callorhinc Leucoraja_ Xenopus_tr Rana_chine Anolis_car Loxodonta_ Monodelphi Dasypus_no Macropus_e Homo_sapie Gallus_gal Mus_muscul Taeniopygi Meleagris_ Ornithorhy Canis_fami Oreochromi);

# 1821 data set
# Sarcopterygii
#my @ingroup = qw( Latimeria_chalumnae Neoceratodus_forsteri Lepidosiren_paradoxa Protopterus_annectens Gallus_gallus Anolis_carolinensis Macropus_eugenii Taeniopygia_guttata Homo_sapiens Xenopus_Silurana_tropicalis Monodelphis_domestica Ornithorhynchus_anatinus Loxodonta_africana Pelodiscus_sinensis Dasypus_novemcinctus Canis_lupus_familiaris Meleagris_gallopavo Mus_musculus);
#my @outgroup = qw( Callorhinchus_milii Lepisosteus_oculatus Danio_rerio Takifugu_rubripes Oreochromis_niloticus);
# Tetrapoda
#my @ingroup = qw( Gallus_gallus Anolis_carolinensis Macropus_eugenii Taeniopygia_guttata Homo_sapiens Xenopus_Silurana_tropicalis Monodelphis_domestica Ornithorhynchus_anatinus Loxodonta_africana Pelodiscus_sinensis Dasypus_novemcinctus Canis_lupus_familiaris Meleagris_gallopavo Mus_musculus);
#my @outgroup = qw( Latimeria_chalumnae Neoceratodus_forsteri Lepidosiren_paradoxa Protopterus_annectens Callorhinchus_milii Lepisosteus_oculatus Danio_rerio Takifugu_rubripes Oreochromis_niloticus);
# T1
#my @ingroup = qw( Neoceratodus_forsteri Lepidosiren_paradoxa Protopterus_annectens Gallus_gallus Anolis_carolinensis Macropus_eugenii Taeniopygia_guttata Homo_sapiens Xenopus_Silurana_tropicalis Monodelphis_domestica Ornithorhynchus_anatinus Loxodonta_africana Pelodiscus_sinensis Dasypus_novemcinctus Canis_lupus_familiaris Meleagris_gallopavo Mus_musculus);
#my @outgroup = qw( Latimeria_chalumnae Callorhinchus_milii Lepisosteus_oculatus Danio_rerio Takifugu_rubripes Oreochromis_niloticus);
# T2
#my @ingroup = qw( Latimeria_chalumnae Gallus_gallus Anolis_carolinensis Macropus_eugenii Taeniopygia_guttata Homo_sapiens Xenopus_Silurana_tropicalis Monodelphis_domestica Ornithorhynchus_anatinus Loxodonta_africana Pelodiscus_sinensis Dasypus_novemcinctus Canis_lupus_familiaris Meleagris_gallopavo Mus_musculus);
#my @outgroup = qw( Neoceratodus_forsteri Lepidosiren_paradoxa Protopterus_annectens Callorhinchus_milii Lepisosteus_oculatus Danio_rerio Takifugu_rubripes Oreochromis_niloticus);
# T3
#my @ingroup = qw( Neoceratodus_forsteri Lepidosiren_paradoxa Protopterus_annectens Latimeria_chalumnae);
#my @outgroup = qw( Latimeria_chalumnae Callorhinchus_milii Lepisosteus_oculatus Danio_rerio Takifugu_rubripes Oreochromis_niloticus Gallus_gallus Anolis_carolinensis Macropus_eugenii Taeniopygia_guttata Homo_sapiens Xenopus_Silurana_tropicalis Monodelphis_domestica Ornithorhynchus_anatinus Loxodonta_africana Pelodiscus_sinensis Dasypus_novemcinctus Canis_lupus_familiaris Meleagris_gallopavo Mus_musculus);

my %fasta_hash;
my $header;
my @seq;
my $seq_length = 0;

# read fasta file with SeqIO
my $seqio_obj = Bio::SeqIO->new('-file' => "<$in_fasta",
								'-format' => "fasta");

while ( my $seq_obj = $seqio_obj->next_seq ) {

    $header = $seq_obj->primary_id;
    my $sequence = $seq_obj->seq;
    
    # sanity check if any sequence has a different length
    if ( $seq_obj->length != $seq_length && $seq_length != 0 ) {
    	
    	print STDERR "ERR: sequences might not be aligned!\n";
    }   
    $seq_length = $seq_obj->length;
    @seq = split (//, $sequence);
	$fasta_hash{$header} = [@seq];
}

my %results_binary = ();
my %results_asymmetrical = ();
my %results_noisy = ();

# iterate through each alignment position (marked by $i)
SITE: for ( my $i=0; $i<$seq_length; $i++ ) {
	
	# initialize variables for each new site
	my $ig_highest_freq_aa = "X";
	my $ig_highest_freq_freq = "0";
	my $og_highest_freq_aa = "X";
	my $og_highest_freq_freq = "0";
	my %proportion_ingroup = ();
	my %proportion_outgroup = ();
	
	# get each site pattern for ingroup and outgroups
	INGROUP: foreach my $ig ( @ingroup ) {

		# get site for tetrapods	
		my $site_i = ${ $fasta_hash{$ig} }[$i];

		# get proportion of amino acids into %proportion_ingroup
		if ( !exists $proportion_ingroup{$site_i} ) {
		
			$proportion_ingroup{$site_i} = 1;
		}
		else {
		
			$proportion_ingroup{$site_i}++;
		}
	}
	# get most frequent ingroup amino acid
	foreach my $kig ( sort { $proportion_ingroup{$b} <=> $proportion_ingroup{$a} } keys %proportion_ingroup ) {
		
		$ig_highest_freq_aa = $kig;
		$ig_highest_freq_freq = $proportion_ingroup{$kig};
		last;
	}
	OUTGROUP: foreach my $og ( @outgroup ) {

		# get site for tetrapods	
		my $site_o = ${ $fasta_hash{$og} }[$i];

		# get proportion of amino acids into %proportion_ingroup
		if ( !exists $proportion_outgroup{$site_o} ) {
		
			$proportion_outgroup{$site_o} = 1;
		}
		else {
		
			$proportion_outgroup{$site_o}++;
		}
	}
	# get most frequent outgroup amino acid
	foreach my $kog ( sort { $proportion_outgroup{$b} <=> $proportion_outgroup{$a} } keys %proportion_outgroup ) {
	
		$og_highest_freq_aa = $kog;
		$og_highest_freq_freq = $proportion_outgroup{$kog};
		last;
	}
	
	# skip cases without replacements and when the most frequent amino acid is < ( 1 - noise )
	next if ( $ig_highest_freq_aa eq $og_highest_freq_aa );
	next if ( 	$ig_highest_freq_freq < ( scalar @ingroup * (1 - $proportion_noise ) ) ||
				$og_highest_freq_freq < ( scalar @outgroup * (1 - $proportion_noise ) ) );

	#print "INGROUP: $ig_highest_freq_aa ($ig_highest_freq_freq)\n";
	#print "OUTGROUP: $og_highest_freq_aa ($og_highest_freq_freq)\n";

	# LOOK FOR BIPARTITION CLASSES

	my $replacement =  "$og_highest_freq_aa" . "=>" . "$ig_highest_freq_aa";

	# BINARY BIPARTITIONS
	if ( $ig_highest_freq_freq == scalar @ingroup &&
		 $og_highest_freq_freq == scalar @outgroup ) {
	
		#print "BINARY: site $i: $og_highest_freq_aa ($og_highest_freq_freq) => $ig_highest_freq_aa ($ig_highest_freq_freq)\n";

		# summarize results
		if ( !exists $results_binary{$replacement} ) {
		
			$results_binary{$replacement} = 1;
		}
		else {
			$results_binary{$replacement}++;
		}
		# go to next site to avoid storing same site in more "relaxed" categories
		next SITE;
	}

	# ASYMMETRICAL BIPARTITIONS: ingroup state is homogeneous and not present in outgroup
	if ( 	$ig_highest_freq_freq == scalar @ingroup &&
			!exists $proportion_outgroup{$ig_highest_freq_aa} )  {
	
		#print "ASYMMETRICAL: site $i: $og_highest_freq_aa ($og_highest_freq_freq) => $ig_highest_freq_aa ($ig_highest_freq_freq)\n";

		# summarize results
		if ( !exists $results_asymmetrical{$replacement} ) {
		
			$results_asymmetrical{$replacement} = 1;
		}
		else {
			$results_asymmetrical{$replacement}++;
		}
		# go to next site to avoid storing same site in more "relaxed" categories
		next SITE;
	}
	
	# NOISY BIPARTITIONS: ingroup state is homogeneous and allowed in outgroup up to 25%
	if ( 	$ig_highest_freq_freq == scalar @ingroup &&
			$proportion_outgroup{$ig_highest_freq_aa} < ( scalar @outgroup * $proportion_noise ) )  {
	
		#print "NOISY: site $i: $og_highest_freq_aa ($og_highest_freq_freq) => $ig_highest_freq_aa ($ig_highest_freq_freq)\n";

		# summarize results
		if ( !exists $results_noisy{$replacement} ) {
		
			$results_noisy{$replacement} = 1;
		}
		else {
			$results_noisy{$replacement}++;
		}
		# go to next site to avoid storing same site in more "relaxed" categories
		next SITE;
	}
	# all other positions are ignored.
	# In previous steps, sites were removed if:
	# (i) have equal most freq amino acid in in- and out-group, (ii) most frequent amino 
	# acids are >75% frequent, (iii) belong to either binary, asymmetrical or noisy bipartitions
}
	
	

# PRINT OUT RESULTS
print "\nBINARY BIPARTITIONS\n";

my $results_binary_count = 0;
foreach my $k ( sort { $results_binary{$b} <=> $results_binary{$a} } keys %results_binary ) {

	print "$k\t$results_binary{$k}\n";
	$results_binary_count += $results_binary{$k};
}

print "\nASYMMETRICAL BIPARTITIONS\n";

my $results_asymmetrical_count = 0;

foreach my $l ( sort { $results_asymmetrical{$b} <=> $results_asymmetrical{$a} } keys %results_asymmetrical ) {

	print "$l\t$results_asymmetrical{$l}\n";
	$results_asymmetrical_count += $results_asymmetrical{$l} ;
}

print "\nNOISY BIPARTITIONS\n";

my $results_noisy_count = 0;

foreach my $m ( sort { $results_noisy{$b} <=> $results_noisy{$a} } keys %results_noisy ) {

	print "$m\t$results_noisy{$m}\n";
	$results_noisy_count += $results_noisy{$m};
}

print "\ntotal BINAR: $results_binary_count\n";
print "total ASYMM: $results_asymmetrical_count\n";
print "total NOISY: $results_noisy_count\n";

print STDERR "\ndone!\n";