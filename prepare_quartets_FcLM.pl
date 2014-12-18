#!/usr/bin/perl

use strict;
use warnings;
use Bio::SeqIO;

# Iker Irisarri, University of Kosntanz, December 2014
# 
# Prepares sets of sequences for Four Cluster Likelihood Mapping (FcLM)
# Given 4 predefined groups, the script generates all possible combinations of taxa (one from each group)
# and extracts sequences from input fasta, creating as many outfiles as necessary (in fasta)
# 4 nested foreach loops produce all possible combinations of elements from each array (checked with smaller arrays)

my $usage = "prepare_quartets_FcLM.pl fasta_file\n";
my $fasta = $ARGV[0] or die $usage;



# define the four groups
my @outgroup = qw ( Danio_reri Takifugu_r Scyliorhin Callorhinc Leucoraja_ Oreochromi );
my @coelacanth= qw ( Latimeria_ );
my @lungfishes = qw ( Lepidosire Protopteru Neoceratod );
my @tetrapods = qw ( Monodelphi Dasypus_no Macropus_e Homo_sapie Gallus_gal Mus_muscul Taeniopygi Meleagris_ Ornithorhy Canis_fami Loxodonta_ Anolis_car Xenopus_tr Rana_chine );

my %quartets;
my $quartet_num = 0;

# try out all the combinations and store them as array into %quartets, using successive numbers as keys
foreach my $o ( @outgroup ) {
    foreach my $c ( @coelacanth ) {
	foreach my $l ( @lungfishes ) {
	    foreach my $t ( @tetrapods ) {
		$quartet_num++;
		$quartets{$quartet_num} = [$o, $c, $l, $t];
	    }
	}
    }
}


# to print out all possible quartets
#foreach my $key ( sort keys %quartets ) {
#    print join (' ', @{ $quartets{$key} }), "\n";
#}

my %fasta;

# read fasta file with SeqIO && store sequences in a hash
my $seqio_obj = Bio::SeqIO->new('-file' => "<$fasta",
				'-format' => "fasta");

while (my $seq_obj = $seqio_obj->next_seq){
    my $name = $seq_obj->primary_id;
    my $seq = $seq_obj->seq;
    $fasta{$name} = $seq;
}

# extract sequences for each quartet && write to new output file (add numbers to filename with count)
my $count = 0;

foreach my $q ( keys %quartets ) {
    $count++;
    my $outfile = "$fasta" . ".$count";
    open (OUT, ">", $outfile);
    foreach my $e ( @{ $quartets{$q} } ) {
	if ( exists ( $fasta{$e} ) ) {
	    print OUT ">$e\n";
	    print OUT "$fasta{$e}\n";
	}
    }
    close(OUT);
}
	
