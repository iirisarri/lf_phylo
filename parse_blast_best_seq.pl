#!/usr/bin/perl -w
use strict;
use Bio::SeqIO;
use Bio::SearchIO;
use Data::Dumper;

# blast parser takes best hit and prints the sequence
# check origin of best seq and create statistics?

my $usage = "parse_blast_best_seq.pl blast_report > output.fa";
my $blast_report = $ARGV[0] or die $usage;

# read blast report
my $report = new Bio::SearchIO(-format => "blast", 
                           -file   => "<$blast_report",
                          );

my %hash;
my $bits0 = 0;
my $hsp_num = 0;

while( my $result = $report->next_result ) {
    # get result info
    $result->query_name =~ /Protopteru_(G\d{5})/g;
    my $query=$1;

    while( my $hit = $result->next_hit ) {
	# get hit info
	my $hit_name = $hit->name;
	my $bits1 = $hit->bits; # bit score of best hsp for that hit                                      
	$hsp_num = 0;

	# process only hsps if bit score of new hit is higher
	# this should not happen because hits are ordered by bit score, but just in case :-)

	if ($bits1 > $bits0 ) {
	    
	    # $bits0 = bit score of previous hit
	    # $bits1 = bit score of current hit
	    # reassing new bits for next comparison

	    $bits0 = $bits1;

	    while( my $hsp = $hit->next_hsp) {

		$hsp_num++;
		my $seq = $hsp->query_string;
		# remove dashes from sequence
		$seq =~ s/-//g;
		my @query_coord = $hsp->range('query');
		my @hit_coord = $hsp->range('hit');

		# create new hash of hashes with the above info

		$hash{$query}{$hit_name}{'bits'} = $bits1;
		$hash{$query}{$hit_name}{'hsp#'.$hsp_num}{'query_coord'} = [ @query_coord ];
		$hash{$query}{$hit_name}{'hsp#'.$hsp_num}{'hit_coord'} = [ @hit_coord ];
		$hash{$query}{$hit_name}{'hsp#'.$hsp_num}{'seq'} = $seq;
	    }
	}
    }
}


#print Dumper \%hash;


foreach my $q ( keys %hash ) {
    my %inner_hash1 = % { $hash{$q}  };
    foreach my $h ( keys %inner_hash1 ) {
	my %inner_hash2 = % { $inner_hash1{$h} };
	if 
	foreach my $p ( keys %inner_hash2 ) {
	    my %inner_hash3 = % { $inner_hash2{$p} };
	    print "$inner_hash3{'query_coord'}[0]\n";
	}
    }
}
#	foreach my $p ( sort { $hash{$q}{$h}{$p}{'hit_coord'}[0]{$a} 
#			<=> $hash{$q}{$h}{$p}{'hit_coord'}[0]{$b} }
#			keys $hash{$q}{$h} ) {


#	    $hash{$q}{$h}{$p}{'hit_coord'}[0];


#    print $keys,  "\n";
# for each query print best hit fasta format
#    print ">", $hash{$keys}[0], "\n";
#    print $hash{$keys}[2],"\n";
#}
__END__
# print sequences, not in order (because $query is alphanumeric)
# latter to be accessed by $query number
foreach my $keys(sort keys %hash ) {
# normal header
#    print ">", $hash{$keys}[1], "_", $hash{$keys}[0], "\n";
# modified header to include species name and gene name
    print ">", "Lepidosire", "_", $hash{$keys}[0], "\n";
#    print ">", "Neoceratod", "_", $hash{$keys}[0], "\n";
    print $hash{$keys}[3],"\n";
}
