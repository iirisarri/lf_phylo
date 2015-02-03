#!/usr/bin/perl
 
use strict;
use Bio::SeqIO;
use Data::Dumper;

# Iker Irisarri. University of Konstanz, Februrary 2015
# prints out name of fasta if at least one species from the three groups is present
# ideally to be run in a for loop for multiple files

my $usage = "filter_fasta_by_taxon_sampling.pl in_fasta > stdout";
my $fasta = $ARGV[0] or die $usage; 

# file in with seqio

#READ_IN: 
my $seqio_obj = Bio::SeqIO->new('-file' => "<$fasta", 
				'-format' => "fasta", 
				'-alphabet' => "protein"
				);

# store sequences in a hash

my %hash;

while (my $inseq = $seqio_obj->next_seq) {

    my $name = $inseq->primary_id;
    my $seq = $inseq->seq;
    # hash where $seqname are keys and take the value sof $seq
    $hash{$name}=[$seq];
}

# check which keys are present

if (exists $hash{'Latimeria_chalumnae'}) {

    if (exists $hash{'Lepidosiren_paradoxa'} || exists $hash{'Neoceratodus_forsteri'} || exists $hash{'Protopterus_annectens'} ) {

	if (exists $hash{'Danio_rerio'} || exists $hash{'Lepisosteus_oculatus'} || exists $hash{'Oreochromis_niloticus'} || exists $hash{'Takifugu_rubripes'} || exists $hash{'Petromyzon_marinus'} ||exists $hash{'Callorhinchus_milii'}) {

	    print "$fasta\n";

	}
    }
}
