#!/usr/bin/perl
 
use strict;
use Bio::SeqIO;
use Data::Dumper;

# WORKS PERFECTLY!!
# select sequences from fasta files ("GXXX_all.fa") that include Amemiya 2013's alns + lungfish matches
# select only one outgroup from the three available

my $usage = "reduce_outgroups_from_fasta.pl infile.fa > outfile.fa";
my $fasta = $ARGV[0] or die $usage; 
my $summary = "reduce_outgroups_from_fasta.summary";

# read file in with seqio

#READ_IN: 
my $seqio_obj = Bio::SeqIO->new('-file' => "<$fasta", 
				'-format' => "fasta", 
				'-alphabet' => "protein"
				);

# store sequences in a hash

my %hash;

while (my $inseq = $seqio_obj->next_seq) {
#	my $name = $inseq->primary_id;
	# reduce number of characters to 10
	# to remove gene name from Lepidosiren and Neoceratodus that I added for the BLAST step	
	$inseq->primary_id=~/(.{10}).*/g;
	my $name = $1;
	my $seq = $inseq->seq;
	# hash where $seqname are keys and take the value sof $seq
	$hash{$name}=[$seq];
}

#print Dumper \%hash;

open (OUT, ">>", $summary);

# check which keys are present (from outgroup) and delete the rest
# order of outgroups preferred: Callorhinc > Scyliorhin > Leucoraja_

# if Callorhinc is present, remove other outgroups and print
if (exists $hash{Callorhinc}) {
	# delete removes a key/value pair from a hash
	delete $hash{Scyliorhin};
	delete $hash{Leucoraja_};
	print OUT $fasta, "\tCallorhinc\n";
	foreach my $name (sort keys %hash) {
		print ">", $name, "\n";
		print $hash{$name}[0], "\n";
	}
# if Callorhinc is not present but Scyliorhin yes, remove other outgroups and print
}	elsif (exists $hash{Scyliorhin}) {
		delete $hash{Leucoraja_};
		print OUT $fasta, "\tScyliorhin\n";
		foreach my $name (sort keys %hash) {
			print ">", $name, "\n";
			print $hash{$name}[0], "\n";
		}
# if only Leucoraja_ is present, then print
	}	else {
			print OUT $fasta, "\tLeucoraja_\n";
			foreach my $name (sort keys %hash) {
			print ">", $name, "\n";
			}
		}



close(OUT);

__END__

# print the rest


__END__

# THIS WORKS FINE but I decided to do it otherwise
# It would print out all sequences that are not either Callorhinc, Scyliorhin or Leucoraja_

while (my $seqio_obj = $seqio_obj->next_seq) {
	if ($seqio_obj->primary_id!~/Callorhinc/g && $seqio_obj->primary_id!~/Scyliorhin/g && $seqio_obj->primary_id!~/Leucoraja_/g) {
# small code just to remove gene names from headers of Lepidosiren and Neoceratodus (that I added for blasting)
		if ($seqio_obj->primary_id=~/(Lepidosire).+/g)	{
			print ">", $1, "\n";
			print $seqio_obj->seq, "\n";
		} elsif ($seqio_obj->primary_id=~/(Neoceratod).+/g) {
			print ">", $1, "\n";
			 $seqio_obj->seq, "\n";
			print $seqio_obj->seq, "\n";
		} else {
			print ">", $seqio_obj->primary_id, "\n";
			print $seqio_obj->seq, "\n";
		}
	}
}