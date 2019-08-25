use strict ; 
use warnings ; 

my $min_prop = 0.05 ; 
my %data ; 

foreach my $file ( `ls Bathy_Lao*` ) { 
	chomp $file ; 
	open IN, "<$file" ;
	<IN> ; 
	while (<IN>) { 
		chomp ; 
		my @split = split ( /\t/, $_ ) ; 

		if ( $split[0] ne $split[2] ) { 
			next ;
		}

		if ( $split[3] < $split[1] ) { 
			my $swap = $split[1] ;
			$split[1] = $split[3] ; 
			$split[3] = $swap ; 
		}

		$data{$split[0]."_".$split[1]."_".$split[2]."_".$split[3]}{"00"} += $split[6] ; 
		$data{$split[0]."_".$split[1]."_".$split[2]."_".$split[3]}{"01"} += $split[7] ;
		$data{$split[0]."_".$split[1]."_".$split[2]."_".$split[3]}{"10"} += $split[8] ;
		$data{$split[0]."_".$split[1]."_".$split[2]."_".$split[3]}{"11"} += $split[9] ;
		$data{$split[0]."_".$split[1]."_".$split[2]."_".$split[3]}{"TOTAL"} += $split[6] + $split[7] + $split[8] + $split[9] ;
		$data{$split[0]."_".$split[1]."_".$split[2]."_".$split[3]}{"DIST"} = $split[3] - $split[1] ; 
	}
	close IN ; 
}

foreach my $s ( sort keys %data ) { 
	my $fd = 0 ; 
	if ( $data{$s}{"00"}/$data{$s}{"TOTAL"} > $min_prop && $data{$s}{"01"}/$data{$s}{"TOTAL"} > $min_prop && $data{$s}{"10"}/$data{$s}{"TOTAL"} > $min_prop && $data{$s}{"11"}/$data{$s}{"TOTAL"} > $min_prop ) { 
		$fd ++ ; 
	}
	print $s, "\t", $data{$s}{"TOTAL"}, "\t", $data{$s}{"DIST"}, "\t", $fd, "\n" ;
}
