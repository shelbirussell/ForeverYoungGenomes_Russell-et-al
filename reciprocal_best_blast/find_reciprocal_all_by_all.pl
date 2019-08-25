use strict ;
use warnings ;
use File::Basename ;
use Sort::Naturally ;

##Usage:
##perl find_reciprocal_all_by_all.pl Svsym_top4 /home/shelbirussell/shelbilrussell@gmail.com/lab/projects/symbiont_coalescent/empirical_samples/new_genomes/final_genomes/Solemya_velum/Solemya_velum_Sym_top4.gff filtered_rbb_hits/ > rbb_all_seqs.output

my$base_seq = $ARGV[0] ;
my$base_gene_set = $ARGV[1] ;
my$blast_dir = $ARGV[2] ;

print "\nbase sequence: ", $base_seq, "\n" ;
print "base gene set gff file: ", $base_gene_set, "\n" ;
print "blast directory: ", $blast_dir, "\n" ;

my$base_genes = read_gff($base_gene_set) ;
my%base_genes = %{$base_genes} ;

opendir(GENOMES, $blast_dir) ;
my@all_files = readdir(GENOMES) ;

##only keep blast hit file names
my@files = () ;
foreach my$i (@all_files) {if ($i =~ m/_hits.txt/) {push @files, $i ;}}
print "blast file count: ", scalar @files, "\n\n" ;

##read in all blast results
my%blasts ;
foreach my$blast (@files) {
  my$reciprocal ;
  if ($blast =~ m/^(.+_vs_.+)_nt_filtered_hits.txt/) {$reciprocal = $1 ;}
  $blasts{$reciprocal} = parse_blast("${blast_dir}/${blast}") ;
}

## Find reciprocal best matches for each reciprocal blast
## record in base_matches if $base_seq is one of pair, otherwise, save to $other_matches
my%base_matches ;
my%other_matches ;
my@spp = () ; # make list of species, in sorted order for output

foreach my$reciprocal (keys %blasts) {
  $reciprocal =~ m/(.+)_vs_(.+)/ ;
  my$sp1 = $1 ;
  my$sp2 = $2 ;
  my$opposite = "${sp2}_vs_${sp1}" ;
	#print $sp1, "\t", $sp2, "\n" ;
  #print $reciprocal, "\t", $opposite, "\n\n" ;

  ##count total genes and not rbb matched genes
  my$total = 0 ;
  my$not_matches = 0 ;

  my%sp1vssp2 = %{$blasts{$reciprocal}} ;   ## sp1vsp2
  my%sp2vssp1 = %{$blasts{$opposite}} ;   ## sp2vsp1

  foreach my$query (keys %{$blasts{$reciprocal}}) {
    $total ++ ;

    my$subject1 = $sp1vssp2{$query} ; ## subject 1 = query for sp 2
    my$subject2 = $sp2vssp1{$subject1} ; ## subject 2 = query for sp 1

    ## Skip genes that didn't have a BLAST hit that passed the filters
    if (! $subject2 || ! $subject1) {next ;}

    #print $subject1, "\t", $subject2, "\n" ;
    #print "QUERY: ", $query, "\n", $sp1, "\t", $subject2, "\n", $sp2, "\t", $subject1, "\n\n" ;

    if ($query eq $subject2) {
      ##record $base_seq gene name as primary level of hash

      ##if $base_seq was one of the two genomes compared in this rbb, record the species and gene name
      if ($sp1 eq $base_seq) {
        $base_matches{$subject1}{$sp2} = $subject2 ;
        if (! grep(/$sp2/, @spp)) {push @spp, $sp2} ;
      }
      elsif ($sp2 eq $base_seq) {
        $base_matches{$subject2}{$sp1} = $subject1 ;
        if (! grep(/$sp1/, @spp)) {push @spp, $sp1} ;
      }

      ##otherwise, record entry separately
      else {
        if ($other_matches{$sp2}{$subject2}) {$other_matches{$sp2}{$subject2}{$sp1} = $subject1 ;}
        else {$other_matches{$sp1}{$subject1}{$sp2} = $subject2 ;}
      }
    }

    else {
      #print "NOT MATCHED QUERY: ", $query, "\n", $sp1, "\t", $subject2, "\n", $sp2, "\t", $subject1, "\n\n" ;
      $not_matches ++ ;
    }
  }
  print "FILE: ", $reciprocal, "\n", "genes not matched: ", $not_matches, " of a total ", $total, " genes\n\n" ;
}

##########
##Check that rbb hits are consistent across all taxa
##I.e., A = B & A = C, so does B = C?
##If not, delete gene from record
foreach my$gene (keys %base_matches) {
  ##make list of homologs matching $base_seq $gene
  my@homologs = () ;
  foreach my$taxon2 (keys %{$base_matches{$gene}}) {push @homologs, "${taxon2}\t${base_matches{$gene}{$taxon2}}" ;}
  ##check that all homolog pairs had each other as their rbb
  foreach my$i (0..$#homologs-1) {
    foreach my$j (1..$#homologs) {
      my@is = split(/\t/, $homologs[$i]) ;
      my@js = split(/\t/, $homologs[$j]) ;
      ##Rbb is consistent if either of these two if statements are true
      if ($other_matches{$is[0]}{$is[1]}{$js[0]} = $js[1]) {
        #print "TRUE!\n" ;
      }
      elsif ($other_matches{$js[0]}{$js[1]}{$is[0]} = $is[1]) {
        #print "ALSO TRUE\n" ;
      }
      ##Rbb has another gene as the hit if not
      else {
        ##some reporting
        if ($other_matches{$is[0]}{$is[1]}{$js[0]}) {
          print "Real ", $is[0], "\t", $is[1], " hit to ", $js[0], ":\n" ;
          print $other_matches{$is[0]}{$is[1]}{$js[0]}, "\n" ;
        }
        if ($other_matches{$js[0]}{$js[1]}{$is[0]}) {
          print "Real ", $js[0], "\t", $js[1], " hit to ", $is[0], ":\n" ;
          print $other_matches{$js[0]}{$js[1]}{$is[0]}, "\n" ;
        }
        ##delete homolog from hash to be printed if homologs mismatch across taxa
        delete $base_matches{$gene} ;
      }
    }
  }
}

####################################################
##Print all genes, regardless of hits across spp
my$output = "multireciprocal.out" ;
open OUT, ">$output" or die "can't open $output" ;
print OUT "#${base_seq}\t", join("\t", @spp), "\n" ;

foreach my$gene (sort keys %base_matches) {
  print OUT $gene, "\t" ;

  foreach my$sp (@spp) {
    if ($base_matches{$gene}{$sp}) {print OUT $base_matches{$gene}{$sp}, "\t", ;}
    else {print OUT ".\t" ;}
  }

  print OUT "\n" ;
}

close OUT ;

###########################################################
##Print only those genes that have orthologs in all spp
##Save data for genes with fewer hits than scalar@spp
my$output2 = "all_out_of_all_multireciprocal.out" ;
open OUT, ">$output2" or die "can't open $output2" ;
print OUT "#${base_seq}_gene\t", join("\t", @spp), "\n" ;

##record number of genes with X number of hits (1->all spp) + their records
my%hit_count ;
foreach my$i (0..$#spp) {
  my$count = $i + 1 ;
  $hit_count{$count}{"COUNT"} = 0 ;
  @{$hit_count{$count}{"LINES2PRINT"}} = () ;
}

foreach my$gene (sort keys %base_matches) {
  my@new_spp = () ;
  my@line = () ;
  push @line, $gene ;

  foreach my$sp (@spp) {
    if (exists $base_matches{$gene}{$sp}) {
      push @new_spp, $sp ;
      push @line, $base_matches{$gene}{$sp} ;
    }
  }

  ##only print line for homolog if all taxa were represented (@spp+reference)
  if (scalar@line == scalar@spp+1) {print OUT join("\t", @line), "\n" ;}
  else {
    ##add count to bin of scalar@line hits
    $hit_count{scalar@line}{"COUNT"} ++ ;
  }
}

close OUT ;

print "number of hits\tcount\n" ;
foreach my$bin (nsort keys %hit_count) {
  print $bin, "\t", $hit_count{$bin}{"COUNT"}, "\n" ;
}
print "\n" ;


###########################################################


sub read_gff {
  open GFF, "<$_[0]" ;

  my %data ;

  while (<GFF>) {

    if ( $_ =~ m/^#/ ) {next ;}

    chomp $_ ;
    my @split = split ( /\t/, $_ ) ;

    my $id = "" ;
    if ( $split[8] =~ m/ID=(.+);/ ) {$id = $1 ;}

    my $gene = "" ;
    if ( $split[8] =~ m/Name=(.+);/ ) {$gene = $1 ;}

    $data{$split[0]}{$id}{"START"} = $split[3] ;
    $data{$split[0]}{$id}{"STOP"} = $split[4] ;
    $data{$split[0]}{$id}{"STRAND"} = $split[6] ;
    $data{$split[0]}{$id}{"GENE"} = $gene ;

  }
  close GFF ;
  return \%data ;
}

sub parse_blast {
  my%hits ;
  open REPORT, "<$_[0]" ;

  while (<REPORT>) {
    if ($_ =~ m/^#/) {next ;}

    my@info = split(/\t/, $_) ;
    my$query = $info[0] ;
    my$subject ;

		if ($info[1] =~ m/#/) {
			my@split = split("#", $info[1]) ;
			$subject = $split[0] ;
		}

		else {$subject = $info[1] ;}

    my$product = $info[2] ;
    $product =~ s/ /_/g ;
    my$pident = $info[4] ;
    my$evalue = $info[12] ;

    $hits{$query} = $subject ;
  }

  close REPORT;

  return \%hits ;
}
