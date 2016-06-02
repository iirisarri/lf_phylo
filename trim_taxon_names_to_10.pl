#!/usr/bin/perl
 
use strict;
use Bio::SeqIO;
use Data::Dumper;

# modified from reduce_outgroups_from_fasta.pl
# trims taxon names to a maximum of 10 characters (between 1 and 10 chars)

my $usage = "trim_taxon_names_to_10.pl infile.fa > outfile.fa";
my $fasta = $ARGV[0] or die $usage;

# read file in with seqio

#READ_IN: 
my $seqio_obj = Bio::SeqIO->new('-file' => "<$fasta", 
				'-format' => "fasta", 
				'-alphabet' => "protein"
				);

while (my $inseq = $seqio_obj->next_seq) {
#	my $name = $inseq->primary_id;
	# reduce number of characters to 10
#	$inseq->primary_id=~/(.{10}).*/g; # matches exactly 10 chars (names, if shorter, are lost)
	$inseq->primary_id=~/(.{1,10}).*/g;
	print ">", $1, "\n";
	print $inseq->seq, "\n";
}


