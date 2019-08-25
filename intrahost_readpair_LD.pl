use strict ;
use warnings ;
use File::Basename ;
use Sort::Naturally ;

#Usage: perl test_reads_for_haplotypes.pl positions2test.txt sample.bam
##positions2test.txt == format for output of intra-host_variants_from_mpileup.pl (_all_AF_count.txt)

### NOTES:
######### Optical duplicate reads are filtered out
######### No need to look for reverse complement matches - they are already complemented to the reference and the flag is marked 16 for reverse complement
######### The remaining duplicate reports for read matches at a site are real - some reads overlap (e.g. PJRI47 ILLUMINA-D00365:314:HA5YEADXX:2:1116:2462:4203) - These are accounted for at line 119-142: if 2 reads in pair report the same site, they are counted as separate reads

my$positions = $ARGV[0] ;
my$bam = $ARGV[1] ;

open POSITIONS, "<$positions" or die "cannot open $positions\n" ;

my%positions ;
# Get positions to test
# Only using minor count because only considering biallelic, so only have 2 haplotype states (associated variants (2 haplotypes) vs. not)
while (<POSITIONS>) {
	chomp ;
	my@split = split(/\t/, $_) ;
	my$scaffold = $split[0] ;
	my$pos = $split[1] ;
	my$ref_base = $split[2] ;
	my$ref_count = $split[3] ;
	my$alt_base = $split[4] ;
	my$alt_count = $split[5] ;
	my$allele_freq = $split[6] ;

	$positions{$scaffold}{$pos}{"REF"} = $ref_base ;
	$positions{$scaffold}{$pos}{"REF_COUNT"} = $ref_count ;
	$positions{$scaffold}{$pos}{"ALT"} = $alt_base ;
	$positions{$scaffold}{$pos}{"ALT_COUNT"} = $alt_count ;
}

close POSITIONS ;


##Now, read through bam file, recording all

open BAM, "samtools view $bam |" ;
my%read_pairs ;
my$dups = 0 ;
my$processed_reads = 0 ;

# Get reads that map to those positions
while (<BAM>) {
	next if (/^(\@)/) ;
	s/\n// ;
	s/\r// ;
	my@sam = split(/\t/, $_) ;
	my$read = $sam[0] ;
	my$flag = $sam[1] ;
	my$scaffold = $sam[2] ;

	##Skip optical/PCR duplicates
	my@dup_flags = (1024, 1040, 1089, 1097, 1105, 1107, 1113, 1123, 1137, 1153, 1161, 1171, 1177, 1185, 1187, 1201) ;
	if (grep {$_ == $flag} @dup_flags) {$dups ++ ; next;}

	##Process reads
	##Collect data from read pairs
	my$results ;
	##if cigar string contains an I or D, indicating an insertion or deletion, respectively
	if ($sam[5] =~ m/[ID]/) {$results = read_with_indel($_, \%positions) ;}
	else {$results = read_no_indel($_, \%positions) ;}
	my%results = %{$results} ;

	##Add read calls to pair
	foreach my$read (keys %results) {
		foreach my$scaffold (keys %{$results{$read}}) {
			foreach my$site (keys %{$results{$read}{$scaffold}}) {
				$read_pairs{$read}{$scaffold}{$site} = $results{$read}{$scaffold}{$site} ;
			}
		}
	}
	$processed_reads ++ ;
}

close BAM ;


##Print read stats
print "Total reads processed: ", $processed_reads, "\n" ;
print "Optical/PCR duplicate reads removed: ", 	$dups, "\n" ;

foreach my$read (keys %read_pairs) {
	foreach my$scaff (keys %{$read_pairs{$read}}) {
		foreach my$pos (keys %{$read_pairs{$read}{$scaff}}) {
		}
	}
}

##Gather genotype calls
my%genotypes ;
foreach my$read (keys %read_pairs) {
	##make list of sites (scaff+position) on each read for pairwise iteration below
	my@sites = () ;
	foreach my$scaffold (keys %{$read_pairs{$read}}) {
		foreach my$site (keys %{$read_pairs{$read}{$scaffold}}) {
			push @sites, "${scaffold},${site}" ;
		}
	}

	##record genotype calls for each pair of sites on each read
	foreach my$position1 ( 0..$#sites-1 ) {
		foreach my $position2 ( $position1+1..$#sites ) {
			my@coord1 = split(",", $sites[$position1]) ;
			my@coord2 = split(",", $sites[$position2]) ;

			my$type = $read_pairs{$read}{$coord1[0]}{$coord1[1]} . $read_pairs{$read}{$coord2[0]}{$coord2[1]} ;
			my$type_rev = $read_pairs{$read}{$coord2[0]}{$coord2[1]} . $read_pairs{$read}{$coord1[0]}{$coord1[1]} ;

			if (exists $genotypes{"${sites[$position1]}\t${sites[$position2]}"}) {
				if (exists $genotypes{"${sites[$position1]}\t${sites[$position2]}"}{$type}) {
					$genotypes{"${sites[$position1]}\t${sites[$position2]}"}{$type} ++ ;
				}
				else {
					$genotypes{"${sites[$position1]}\t${sites[$position2]}"}{$type} = 1 ;
				}
			}
			elsif (exists $genotypes{"${sites[$position2]}\t${sites[$position1]}"}) {
				if (exists $genotypes{"${sites[$position2]}\t${sites[$position1]}"}{$type}) {
					$genotypes{"${sites[$position2]}\t${sites[$position1]}"}{$type_rev} ++ ;
				}
				else {
					$genotypes{"${sites[$position2]}\t${sites[$position1]}"}{$type_rev} = 1 ;
				}
			}
			else {
				$genotypes{"${sites[$position1]}\t${sites[$position2]}"}{$type} = 1 ;
			}
		}
	}
}

##Calculate LD statistics across reads (i.e., samples) at each set of pairwise sites
my$site_pairs = scalar (keys %genotypes) ;

##print to output file
my$out = basename($bam) ;
$out =~ s/.bam/_LD.out/ ;
open OUT, ">$out" or die "cannot open $out\n" ;
print OUT "#scaffold1\tsite1\tscaffold2\tsite2\ttotal\tx00\t00\t01\t10\t11\tfreq0\tfreq1\td\tr\n" ;


foreach my$pair (keys %genotypes) {
	my@paired_sites = split(/\t/, $pair) ;
	my@site1 = split(/,/, $paired_sites[0]) ;
	my@site2 = split(/,/, $paired_sites[1]) ;

	##add up total genotype counts
	my$total = 0 ;

	foreach my$type (keys %{$genotypes{$pair}}) {$total += $genotypes{$pair}{$type} ;}

	if ($total < 1) {print "zero genotypes recorded for paired sites: ", $pair, "\n" ;}

	##add in zeros for missing $genotypes
	if (! exists $genotypes{$pair}{"00"}) {$genotypes{$pair}{"00"} = 0 ;}
	if (! exists $genotypes{$pair}{"01"}) {$genotypes{$pair}{"01"} = 0 ;}
	if (! exists $genotypes{$pair}{"10"}) {$genotypes{$pair}{"10"} = 0 ;}
	if (! exists $genotypes{$pair}{"11"}) {$genotypes{$pair}{"11"} = 0 ;}

	my$x00 = $genotypes{$pair}{"00"} / $total ;
	my$freq0 = ($genotypes{$pair}{"00"} + $genotypes{$pair}{"01"} ) / $total ;
	my$freq1 = ($genotypes{$pair}{"00"} + $genotypes{$pair}{"10"} ) / $total ;

	if ( $freq0 == 0 || $freq0 == 1 || $freq1 == 0 || $freq1 == 1 ) {next ;}

	my$d = $x00 - $freq0 * $freq1 ;
	my$sq = sqrt($freq0*(1-$freq0)*$freq1*(1-$freq1)) ;
	my$r = "NA" ;
	if ( $sq > 0 ) {$r = $d/$sq ;}


	print OUT $site1[0], "\t", $site1[1], "\t", $site2[0], "\t", $site2[1], "\t" ;
	print OUT $total, "\t", $x00, "\t", $genotypes{$pair}{"00"}, "\t", $genotypes{$pair}{"01"}, "\t", $genotypes{$pair}{"10"}, "\t", $genotypes{$pair}{"11"},  "\t", $freq0, "\t", $freq1, "\t", $d, "\t", $r, "\n" ;
}

close OUT ;



sub read_no_indel {
	my@sam = split(/\t/, $_[0]) ;
	my%positions = %{$_[1]} ;

	my$read = $sam[0] ;
	my$scaffold = $sam[2] ;
	my$read_start = $sam[3] ;
	my$sequence = $sam[9] ;

	##Get read sequence and coordinates
	my@sequence = split("", $sequence) ;
	my$read_end = $read_start + length($sequence)-1 ;
	my@range = $read_start .. $read_end ;

	##iterate though all pairwise read sites, looking for reads with 2+ variant sites
	##record into hash
	my%call ;
	foreach my $position1 ( 0..$#range-1 ) {
		foreach my $position2 ( $position1+1..$#range ) {
			##Check if position 1 and 2 are recorded as variant sites
			if ($positions{$scaffold}{$range[$position1]} && $positions{$scaffold}{$range[$position2]}) {
				##record calls on read
				##the REF : cases
				if ($sequence[$position1] eq $positions{$scaffold}{$range[$position1]}{"REF"}) {
					$call{$read}{$scaffold}{$range[$position1]} = 0 ;
				}
				if ($sequence[$position2] eq $positions{$scaffold}{$range[$position2]}{"REF"}) {
					$call{$read}{$scaffold}{$range[$position1]} = 0 ;
				}

				##the ALT : cases
				if ($sequence[$position1] eq $positions{$scaffold}{$range[$position1]}{"ALT"}) {
					$call{$read}{$scaffold}{$range[$position1]} = 1 ;
				}
				if ($sequence[$position2] eq $positions{$scaffold}{$range[$position2]}{"ALT"}) {
					$call{$read}{$scaffold}{$range[$position1]} = 1 ;
				}
			}
		}
	}

	return \%call ;
}

sub read_with_indel {
	my@sam = split(/\t/, $_[0]) ;
	my%positions = %{$_[1]} ;

	my$read = $sam[0] ;
	my$flag = $sam[1] ;
	my$scaffold = $sam[2] ;
	my$read_start = $sam[3] ;
	my$cigar = $sam[5] ;
	my$sequence = $sam[9] ;

	## Parse indels (add/subtract them from postions to recreate reference sequence length)
	my%indels ;
	#my@insertions = () ;
	#my@deletions = () ;
	my@cigars = () ;
	##Expand cigars into list matching @sequence site by site
	my@expanded_cigars = () ;

	if ($cigar =~ m/[ID]/) {
		while ($cigar =~ m/([0-9]+[MIDS=X])/g) {push @cigars, $1 ;}

		my$read_counter = 0 ;
		my$new_sequence ;
		foreach my$cigar (@cigars) {
			my$bp = $cigar ;
			$bp =~ s/([MIDS=X])//g ;
			my$type = $1 ;
			my$current_ref_position = $read_start + $read_counter ;

			if ($type eq "M" || $type eq "=" || $type eq "S" || $type eq "X") {
				$new_sequence .= substr($sequence, $read_counter, $bp) ;
				$read_counter += $bp ;
				foreach my$i (1..$bp) {push @expanded_cigars, $type ;}
			}

			elsif ($type eq "I") {
				##record indel
				##keep new bases (read insertion sequence "doesn't consume the reference" according to the SAM format). Reference base at start of insertion will be added in later.
				$indels{$scaffold}{$current_ref_position}{"INSERTION"} = substr($sequence, $read_counter, $bp) ;
				#skip adding indel sequence to new sequence to generate reference length
				$read_counter += $bp ;
			}

			else {
				#add Ds to new sequence for deleted region
				foreach my$i (1..$bp) {
					$new_sequence .= "D" ;
					push @expanded_cigars, $type ;
				}
				##we don't have the deletion sequence without reading in the reference fasta, so this only records length of deletion (not including ref base still at site - that's added in later)
				$indels{$scaffold}{$current_ref_position}{"DELETION"} =  $bp ;
				#don't increment counter because deletions aren't in read sequence ;
			}
		}
		##update read $sequence to match reference sequence length (no insertions with deleted segments marked by Ds)
		$sequence = $new_sequence ;
	}

	##########
	##Search read sequence for variant sites
	my@sequence = split("", $sequence) ;
	my$read_end = $read_start + length($sequence)-1 ;
	my@range = $read_start .. $read_end ;

	#print "compare sequence length: ", scalar@sequence, " to range length: ", scalar@range, " to cigar list length: ", scalar@cigars, " to expanded cigar list length: ", scalar@expanded_cigars, "\n" ;
	#print "cigar string: ", join(",", @cigars), "\n\n" ;

	##iterate though all pairwise read sites, looking for reads with 2+ variant sites
	##record into hash
	my%call ;
	foreach my $position1 ( 0..$#range-1 ) {
		foreach my $position2 ( $position1+1..$#range ) {

			##These are indels,  so what is the variant type? Insertion, deletion, or snp?
			##Approach: check the cigar string to know how to deal with variation
			##Matches are checked for snps, insertions are checked for identical indel sequence, and deletions are checked for start site and length relative to the reference

			##Check if position 1 and 2 are recorded as variant sites
			if ($positions{$scaffold}{$range[$position1]} && $positions{$scaffold}{$range[$position2]}) {
				##record position1 and position2 genotypes for later entry
				my$geno1 ;
				my$geno2 ;

				##check for insertions first, as they aren't in the new read sequence or expanded_cigar list
				if (exists $indels{$scaffold}{$range[$position1]}{"INSERTION"}) {
					my$full_insertion = $positions{$scaffold}{$range[$position1]}{"REF"} . $indels{$scaffold}{$range[$position1]}{"INSERTION"} ;
					if ($positions{$scaffold}{$range[$position1]}{"ALT"} eq $full_insertion) {
						$geno1 = "1" ;
					}
					else {$geno1 = "0" ;}
				}
				elsif (exists $indels{$scaffold}{$range[$position1]}{"DELETION"}) {
					##if length of deleted reference sequence in vcf equals length of deleted sequence on read, the deletion records match
					if (length($positions{$scaffold}{$range[$position1]}{"REF"}) == ($indels{$scaffold}{$range[$position1]}{"DELETION"}+1)) {
						$geno1 = "1" ;
					}
					else {$geno1 = "0" ;}
				}
				else {
					##check that snp identity matches
					if ($positions{$scaffold}{$range[$position1]}{"REF"} eq $sequence[$position1]) {
						$geno1 = "0" ;
					}
					elsif ($positions{$scaffold}{$range[$position1]}{"ALT"} eq $sequence[$position1]) {
						$geno1 = "1" ;
					}
					else {
						#print "NOT MATCHED: read variant of type ", $expanded_cigars[$position1], " with sequence ", $sequence[$position1], " (REF:", $positions{$scaffold}{$range[$position1]}{"REF"}, " ALT:", $positions{$scaffold}{$range[$position1]}{"ALT"}, ")", "\n" ;
					}
				}

				if (exists $indels{$scaffold}{$range[$position2]}{"INSERTION"}) {
					my$full_insertion = $positions{$scaffold}{$range[$position2]}{"REF"} . $indels{$scaffold}{$range[$position2]}{"INSERTION"} ;
					if ($positions{$scaffold}{$range[$position2]}{"ALT"} eq $full_insertion) {
						$geno2 = "1" ;
					}
					else {$geno2 = "0" ;}
				}
				elsif (exists $indels{$scaffold}{$range[$position2]}{"DELETION"}) {
					##if length of deleted reference sequence in vcf equals length of deleted sequence on read, the deletion records match
					if (length($positions{$scaffold}{$range[$position2]}{"REF"}) == ($indels{$scaffold}{$range[$position2]}{"DELETION"}+1)) {
						$geno2 = "1" ;
					}
					else {$geno2 = "0" ;}
				}
				else {
					##check that snp identity matches
					if ($positions{$scaffold}{$range[$position2]}{"REF"} eq $sequence[$position2]) {
						$geno2 = "0" ;
					}
					elsif ($positions{$scaffold}{$range[$position2]}{"ALT"} eq $sequence[$position2]) {
						$geno2 = "1" ;
					}
					else {
						#print "NOT MATCHED: read variant of type ", $expanded_cigars[$position2], " with sequence ", $sequence[$position2], " (REF:", $positions{$scaffold}{$range[$position2]}{"REF"}, " ALT:", $positions{$scaffold}{$range[$position2]}{"ALT"}, ")", "\n" ;
					}
				}

				##Record genotypes to hash if genotypes at both sites were recorded (i.e., biallelic)
				if ($geno1 && $geno2) {
					$call{$read}{$scaffold}{$range[$position1]} = $geno1 ;
					$call{$read}{$scaffold}{$range[$position2]} = $geno2 ;
				}
			}
		}
	}

	return \%call ;
}
