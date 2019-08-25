use strict ;
use warnings ;
use Sort::Naturally ;

##Purpose: Read in Mauve xmfa alignment blocks,
##remove blocks with few samples,
##and output phylip file for tree inference

my$xmfa = $ARGV[0] ;
open XMFA, "<$xmfa" or die "cannot open $xmfa\n" ;

my%alignment ;
my$block = 0 ;
my$current_seq = "" ;

while (<XMFA>) {
  chomp ;

  ##Set up names
  if ($_ =~ m/^#Sequence(\d+)File/) {
    my$number = $1 ;
    my@split = split(/\t/, $_) ;
    my$sample = $split[1] ;
    $sample =~ s/.+\/(.+)$/$1/g ;

    ##eliminate some characters to make names shorter
    $sample =~ s/^NC_// ;
    $sample =~ s/.gb$// ;
    $sample =~ s/_mitochondrian$// ;
    $sample =~ s/\(reversed\)// ;
    $sample =~ s/_mito// ;
    $sample =~ s/\s+// ;
    $sample =~ s/Solemya_/S/ ;
    $sample =~ s/Calyptogena_/C/ ;
    $sample =~ s/Bathymodiolus_/B/ ;
    ##after all that, now just lop off the end of the header to =9+space before alignment starts
    $sample =~ s/([A-Za-z0-9_]{9}).+/$1/ ;

    ##check that all sample names are <= 10 characters long
    if (length($sample) > 10) {print "SHORTEN: ", $sample, "\n" ;}

    print $number, "\t", $sample, "\n" ;
    $alignment{$number}{"NAME"} = $sample ;
  }

  ##Skip other commented lines
  if ($_ =~ m/^#/) {next;}

  ##Start new seq with new header
  elsif ($_ =~ m/^>/) {
    my$header = $_ ;
    $header =~ s/// ;
    $current_seq = $header ;

    ##get sample number
    my@info = split(/\s+/, $header) ;
    my$align = $info[1] ;
    my@sample = split(/:/, $align) ;

    ##initialize sequence for sample and block
    $alignment{$sample[0]}{"ALIGNMENTS"}{$block} = "" ;
  }

  ##End alignment block
  elsif ($_ =~ m/^=/) {$block ++ ;}

  ##Read lines of alignment into record for $current_seq and $block
  else {
    ##get sample number
    my@info = split(/\s+/, $current_seq) ;
    my$align = $info[1] ;
    my@sample = split(/:/, $align) ;
    $alignment{$sample[0]}{"ALIGNMENTS"}{$block} .= $_ ;
  }
}

close XMFA ;

##Iterate through alignments and delete blocks with sequence content for less than 100% of sequences
my%blocks ;
foreach my$sample (keys %alignment) {
  foreach my$block (keys %{$alignment{$sample}{"ALIGNMENTS"}}) {
    if ($alignment{$sample}{"ALIGNMENTS"}{$block} =~ m/[ATCG]+/g) {
      if (exists $blocks{$block}) {$blocks{$block} ++ ;}
      else {$blocks{$block} = 1 ;}
    }
  }
}

my$total_samples = scalar keys %alignment ;
foreach my$block (keys %blocks) {
  print "Number of samples in block ${block}: ", $blocks{$block}, "\n" ;
  ##check for blocks with few samples aligned
  if ($blocks{$block} < $total_samples) {
    ##and delete them from %alignment
    foreach my$sample (keys %alignment){
      delete $alignment{$sample}{"ALIGNMENTS"}{$block} ;
    }
  }
}


##Concatenate alignments for each sample
my%concat_alignments ;
my$full_length_alignment ;

foreach my$sample (keys %alignment) {
  $concat_alignments{$sample} = "" ;
  #print "SAMPLE: ", $sample, "\t", $alignment{$sample}{"NAME"}, "\n" ;
  foreach my$block (nsort keys %{$alignment{$sample}{"ALIGNMENTS"}}) {
    $concat_alignments{$sample} .= $alignment{$sample}{"ALIGNMENTS"}{$block} ;
  }
  if ($full_length_alignment) {
    if (length($concat_alignments{$sample}) == $full_length_alignment) {next ;}
    else {
      print "sample ", $sample, " alignment length: ", length($concat_alignments{$sample}), " does not equal global alignment length: ", $full_length_alignment, "!\n" ;
    }
  }
  else {
    $full_length_alignment = length($concat_alignments{$sample}) ;
    print "FULL ALIGNMENT LENGTH: ", length($concat_alignments{$sample}), "\n" ;
  }
}


##Print phylip format
my$phylip = $xmfa ;
$phylip =~ s/.xmfa/.phylip/ ;
open OUT, ">$phylip" or die "cannot open $phylip\n" ;
print OUT "\t", $total_samples, "\t", $full_length_alignment, "\n" ;

foreach my$sample (nsort keys %concat_alignments) {
  print OUT $alignment{$sample}{"NAME"}, " ", $concat_alignments{$sample}, "\n" ;
}

close OUT ;
