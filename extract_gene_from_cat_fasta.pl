#!/usr/bin/perl
 
use strict;
use Bio::SeqIO;

# extract protopterus sequences from fasta files by matching seq name
# use single gene alignments from Amemiya 2013 for blasting against lf transcriptomes

my $usage = "extract_proann_from_fasta.pl infile.fa gene_name> outfile.fa";
my $fasta = $ARGV[0] or die $usage; 
my $gene = $ARGV[1] or die $usage;

# read file in with seqio
my $seqio_obj = Bio::SeqIO->new('-file' => "<$fasta", 
				'-format' => "fasta", 
				'-alphabet' => "protein"
				);

while (my $seqio_obj = $seqio_obj->next_seq) {
#    print $seqio_obj->primary_id, "\n";
    if ($seqio_obj->primary_id=~/$gene/g) {
	print ">", $seqio_obj->primary_id, "\n";
	print $seqio_obj->seq, "\n";
    }
}


