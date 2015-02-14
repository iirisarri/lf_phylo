#!/usr/bin/perl

# local bioperl library for hpc2
use lib "/home/iirisar/iirisar/lib/perl5/";
use Bio::Perl;                                                                             


use strict;
use Bio::SeqIO;
use Data::Dumper;

# Iker Irisarri, University of Konstanz. Feb 2015
# modif of reduce_outgroups_from_fasta.pl

my $usage = "reduce_outgroups_from_fasta.pl infile.fa > outfile.fa";
my $fasta = $ARGV[0] or die $usage; 


# read file in with seqio

#READ_IN: 
my $seqio_obj = Bio::SeqIO->new('-file' => "<$fasta", 
				'-format' => "fasta", 
				'-alphabet' => "protein"
				);

# store sequences in a hash

my %hash;

while (my $inseq = $seqio_obj->next_seq) {

	#$inseq->primary_id=~/(.{10}).*/g;
	#my $name = $1;
    my $name = $inseq->primary_id;
    my $seq = $inseq->seq;
	# hash where $seqname are keys and take the value sof $seq
	$hash{$name}=[$seq];
}

#print Dumper \%hash;


if (exists $hash{Callorhinchus_milii}) {

    # re-assing callorhinchus to a new key "outgroup" and remove callorhinchus
    $hash{outgroup}=$hash{Callorhinchus_milii};

    delete $hash{Callorhinchus_milii};

    foreach my $name (sort keys %hash) {
	print ">", $name, "\n";
	print $hash{$name}[0], "\n";
    }

}

elsif (exists $hash{Lepisosteus_oculatus}) {

    $hash{outgroup}=$hash{Lepisosteus_oculatus};

    delete $hash{Lepisosteus_oculatus};

    foreach my $name (sort keys %hash) {
        print ">", $name, "\n";
        print $hash{$name}[0], "\n";
    }

}

elsif (exists $hash{Danio_rerio}) {

    $hash{outgroup}=$hash{Danio_rerio};

    delete $hash{Danio_rerio};

    foreach my $name (sort keys %hash) {
	print ">", $name, "\n";
	print $hash{$name}[0], "\n";
    }

}
else {

    print STDERR "OUTGROUP MISSING!!\n";

}



