use strict ; 
use warnings ; 

my %prop ; 
my %total ;

while (<STDIN>) { 

	chomp ; 

	my @split = split ( /\t/, $_ ) ; 

	if ( $split[2] < 10 ) { 
		$total{10} ++ ; 
		$prop{10} += $split[3] ; 
	}
	elsif ( $split[2] < 100 ) { 
		$total{100} ++ ; 
		$prop{100} += $split[3] ; 
	}
	elsif ( $split[2] < 1000 ) { 
		$total{1000} ++ ; 
		$prop{1000} += $split[3] ; 
	}
}

foreach my $bin ( sort {$a<=>$b} keys %prop ) { 
	print $bin, "\t", $prop{$bin}/$total{$bin}, "\t", $total{$bin}, "\n" ;
}
