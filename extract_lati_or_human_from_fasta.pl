#!/usr/bin/perl -w

use strict;

use Bio::SeqIO;
use Data::Dumper;

# Iker Irisarri, Univeristy of Konstanz, January 2014
# modification of extract_from_fasta_by_name.pl to:
# (1) add gene names to fasta headers
# (2) extract lati sequences, and if not present, xenopus, human or mouse sequences

# example of usage:
# files=*.fa
# for f in $files; do perl /computation/iker/software/scripts/phylogm/extract_lati_or_human_from_fasta.pl $f ; done > out

my $usage = "extractSeqFromFasta fasta_file > stdout\n";
my $fasta = $ARGV[0] or die $usage;

# read fasta file with SeqIO
my $seqio_obj = Bio::SeqIO->new('-file' => "<$fasta",
                		        '-format' => "fasta");

# get gene name
$fasta =~ /(\w+)\.fa\w*/;
my $name = $1;

my %info;           		        

while (my $seq_obj = $seqio_obj->next_seq){

    my $seqname = $seq_obj->primary_id;
    my $seq = $seq_obj->seq;

    $info{$seqname} = $seq;

}


if ( exists $info{"Latimeria_chalumnae"} ) {

    print ">Latimeria_chalumnae", "_", "$name\n";
    print "$info{'Latimeria_chalumnae'}\n";

   }
    
elsif ( exists $info{"Xenopus_Silurana_tropicalis"} ) {

    print ">Xenopus_Silurana_tropicalis", "_", "$name\n";
    print "$info{'Xenopus_Silurana_tropicalis'}\n";

   }   

elsif ( exists $info{"Homo_sapiens"} ) {

    print ">Homo_sapiens", "_", "$name\n";
    print "$info{'Homo_sapiens'}\n";

}

elsif ( exists $info{"Mus_musculus"} ) {

    print ">Mus_musculus", "_", "$name\n";
    print "$info{'Mus_musculus'}\n";

   }

else {

    print "ERROR: sequence not found for file $fasta\n";

}



