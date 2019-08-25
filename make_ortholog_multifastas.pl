use strict ;
use warnings ;
use File::Basename ;
use Sort::Naturally ;

## Usage: perl make_ortholog_multifastas.pl multireciprocal-4of4.out nt_dir/

my$ortholog_file = $ARGV[0] ;
my$nt_dir = $ARGV[1] ;

opendir(SEQS, $nt_dir) ;
my@all_files = readdir(SEQS) ;

##only keep fasta file names
my@files = () ;
foreach my$i (@all_files) {if ($i =~ m/.fasta/) {push @files, $i ;}}

open IN, "<$ortholog_file" ;
my%orthologs ;
my$ortholog_count = 0 ;

while (<IN>) {
  chomp ;
  if ($_ =~ m/^#/) {
    my$header = $_ ;
    $header =~ s/^#// ;
    my@split = split(/\t/, $header) ;
    ##name spp in hash by their column number so subsequent rows get entered appropriately
    foreach my$column (0..$#split) {$orthologs{$column}{"NAME"} = $split[$column] ;}
  }

  else {
    ##Add gene info to the correct species entry by using column order
    my@split = split(/\t/, $_) ;
    foreach my$spp (0..$#split) {
      ##skip blank entries
      if ($split[$spp] eq ".") {next ;}
      else {$orthologs{$spp}{"GENES"}{$ortholog_count} = $split[$spp] ;}
    }
    $ortholog_count ++ ;
  }
}

close IN ;


my%output ;

foreach my$spp (nsort keys %orthologs) {
#  #print $orthologs{$spp}{"NAME"}, "\n" ;
  foreach my$file (@files) {
    if ($file =~ m/$orthologs{$spp}{"NAME"}/) {
      ##File name check
      #print "FILE: ", $file, " matches header: ", $orthologs{$spp}{"NAME"}, "\n\n" ;

      ##Read fasta nt seq file in for spp and record seq for each ortholog to hash
      my$seq = read_fasta($nt_dir,$file) ;
      my%seq = %{$seq} ;

      foreach my$gene (nsort keys %{$orthologs{$spp}{"GENES"}}) {
        foreach my$header (keys %seq) {
          ##eq not =~ m// because enough of headers are long enough to match a substantial proportion of mismatched gene names
          if ($header eq $orthologs{$spp}{"GENES"}{$gene}) {
            #print $header, "\t", $orthologs{$spp}{"GENES"}{$gene}, "\n" ;
            ##record sequences in hash by ortholog number
            ##no need to retain fasta header, as the hash info will be used to make a new one
            $output{$gene}{$orthologs{$spp}{"NAME"}}{"SEQ"} = $seq{$header} ;
            $output{$gene}{$orthologs{$spp}{"NAME"}}{"NAME"} = $orthologs{$spp}{"GENES"}{$gene} ;
          }
          else {
            #print $orthologs{$spp}{"NAME"}, "\t", $gene, "\t", $header, " does not match ", $orthologs{$spp}{"GENES"}{$gene}, "\n\n" ;
          }
        }
      }
    }
  }
}

foreach my$gene (keys %output) {
  open OUT, ">ortholog_${gene}_nt.fasta" or die "cannot open ortholog_${gene}_nt.fasta\n" ;
  print "Gene: ", $gene, " represented by ", scalar(keys%{$output{$gene}}), " species\n" ;

  foreach my$spp (nsort keys %{$output{$gene}}) {
    print OUT ">", ${spp}, "_", $output{$gene}{$spp}{"NAME"}, "\n" ;
    print OUT $output{$gene}{$spp}{"SEQ"}, "\n" ;
  }

  close OUT ;
}



sub read_fasta {
  open FASTA, "<$_[0]/$_[1]" ;
  my%seqs ;
  my$header ;
  my$seq ;

  while (<FASTA>) {
    chomp ;
    if ( $_ =~ m/^#/ ) {next ;}
    if ( $_ =~ m/>/ ) {
      if ($seq) {$seqs{$header} = $seq ;}
      my@full_header = split(/\s+/, $_) ;
      $header = $full_header[0] ;
      $header =~ s/^>// ;

      $seq = "" ;
    }

    else {
      $_ =~ s/\s+//g ;
      $seq .= $_ ;
    }
  }

  close FASTA ;

  if ($seq) {$seqs{$header} = $seq ;}

  return \%seqs ;
}
