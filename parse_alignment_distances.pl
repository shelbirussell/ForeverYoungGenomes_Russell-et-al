use strict ;
use warnings ;
use Sort::Naturally ;

##Usage:perl /home/shelbirussell/shelbilrussell@gmail.com/lab/projects/symbiont_coalescent/empirical_samples/cs_population_analyses/relative_divergence/perl_scripts/parse_alignment_distances.pl /home/shelbirussell/shelbilrussell@gmail.com/lab/projects/symbiont_coalescent/empirical_samples/cs_population_analyses/relative_divergence/taxa_traits.txt ../4.5-ortholog_distances/

my$taxa = $ARGV[0] ;
my$matrix_dir = $ARGV[1] ;

print "TAXON REFERENCE FILE: ", $taxa, "\n\n" ;
print "MATRIX FILE DIRECTORY: ", $matrix_dir, "\n\n" ;

open TAXA, "<$taxa" or die "cannot open $taxa\n" ;
my%record ;

while (<TAXA>) {
  chomp ;
  if ($_ =~ m/^#/) {next;}
  my@split = split(/\s+/, $_) ;
  $record{$split[0]}{"LIFESTYLE"} = $split[1] ;
  $record{$split[0]}{"SYM_CLADE"} = $split[2] ;
}

close TAXA ;

opendir(MATRICIES, $matrix_dir) ;
my@matrix_files = readdir(MATRICIES) ;

foreach my$i (0..$#matrix_files) {
  if ($matrix_files[$i] =~ m/\.csv/) {
    print $matrix_files[$i], "\n" ;

    ##read in matrix
    my%matrix ;
    open MATRIX, "<$matrix_dir/$matrix_files[$i]" or die "cannot open $matrix_dir/$matrix_files[$i]\n" ;
    while (<MATRIX>) {
      chomp ;
      my@split = split(',', $_) ;

      ##top row and first column (i.e., first entry of every line) are taxon ids
      ##top row starts with an empty tab
      if ($_ =~ m/^,/) {
        ##enter each name in order by column number to reference on subsequent rows
        foreach my$column (1..$#split) {
          #print "SAMPLE NAME ", $column, " : ", $split[$column], "\n" ;
          my$sample = $split[$column] ;
          $sample = s/ \(reversed\)// ;
          $matrix{$column}{"NAME"} = $split[$column] ;
        }
      }
      ##subsequent lines: percent identities get entered into hash
      else {
        foreach my$column (1..$#split) {
          $matrix{$column}{"PERCIDS"}{$split[0]} = $split[$column] ;
          #print "PAIRWISE ID BTWN ", $matrix{$column}{"NAME"}, " and ", $split[0], " = ", $split[$column], "\n" ;
        }
      }
    }

    ##record data in %matrix to %record
    ##$column1 corresponds with $taxon1, and is the taxon the record is recorded under
    foreach my$column (keys %matrix) {
      ##Find taxon name in %record that matches matrix column taxon name in %matrix
      foreach my$taxon (keys %record) {
        if ($matrix{$column}{"NAME"} =~ m/$taxon/) {
          ##record % identities for each gene for each taxon pair by the file number
          ##this way, an array doesn't have to be used
          ##record formatting is: taxon1 -> IDENTITIES -> taxon2 -> alignment# -> %id
          foreach my$taxon2 (keys %{$matrix{$column}{"PERCIDS"}}) {
            ##get taxon name that matches gene name (first half of gene name)
            foreach my$tax (keys %record) {
              if ($taxon2 =~ m/$tax/) {
                $record{$taxon}{"IDENTITIES"}{$tax}{$i} = $matrix{$column}{"PERCIDS"}{$taxon2} ;
              }
            }
          }
        }
      }
    }
  }
}

##Make hash to organize records by sym clade
my%groups ;
foreach my$taxon (keys %record) {
  @{$groups{$record{$taxon}{"SYM_CLADE"}}{"SAME_GROUP"}} = () ;
  @{$groups{$record{$taxon}{"SYM_CLADE"}}{"ALLSYMS"}} = () ;
  @{$groups{$record{$taxon}{"SYM_CLADE"}}{"FREELIVING"}} = () ;
}

foreach my$taxon (keys %record) {
  foreach my$category (keys %{$record{$taxon}}) {
    #print $taxon, "\t", $category, "\n" ;
  }
  foreach my$taxon2 (keys %{$record{$taxon}{"IDENTITIES"}}) {
        #print $taxon2, "\n" ;

  }
}

#Record values
##iterate through $taxon1 and $taxon2 pairs, and record %ID according to clade and lifestyle
foreach my$taxon1 (nsort keys %record) {
  foreach my$taxon2 (nsort keys %{$record{$taxon1}{"IDENTITIES"}}) {
    ##Skip self
    if ($taxon1 eq $taxon2) {next;}
    #print $taxon1, "\t", $taxon2, "\n" ;
    my$taxon1_group = $record{$taxon1}{"SYM_CLADE"} ;
    my$taxon2_group = $record{$taxon2}{"SYM_CLADE"} ;

    foreach my$gene (nsort keys %{$record{$taxon1}{"IDENTITIES"}{$taxon2}}) {
      ##record %IDs for pairs belonging to same clade
      ##don't add same group syms to LIFESTYLE="symbiont" by making elsif statements
      ##(restart if on second, if you want to count same group symbionts counted twice)
      if ($taxon1_group eq $taxon2_group) {
        push @{$groups{$taxon1_group}{"SAME_GROUP"}}, 100-$record{$taxon1}{"IDENTITIES"}{$taxon2}{$gene} ;
      }
      ##record IDs to lists for taxon1 clade compared to free-living and symbiont
      elsif ($record{$taxon2}{"LIFESTYLE"} eq "free-living") {
        push @{$groups{$taxon1_group}{"FREELIVING"}}, 100-$record{$taxon1}{"IDENTITIES"}{$taxon2}{$gene} ;
      }
      elsif ($record{$taxon2}{"LIFESTYLE"} eq "symbiont") {
        push @{$groups{$taxon1_group}{"ALLSYMS"}}, 100-$record{$taxon1}{"IDENTITIES"}{$taxon2}{$gene} ;
      }
      else {
        print "NO match for: ", $taxon2, "\t", $record{$taxon2}{"LIFESTYLE"}, "\n" ;
      }
    }
  }
}


my$output = "parsed_alignment_matrices-divergence.txt" ;
open OUT, ">$output" or die "cannot open $output\n" ;

##Print a tab delimited row for each group
##Enter "NA"s so all rows are equal length
my$max_length = 0 ;
foreach my$group (keys %groups) {
  foreach my$division (keys %{$groups{$group}}) {
    if ($max_length < scalar(@{$groups{$group}{$division}})) {
      $max_length = scalar(@{$groups{$group}{$division}}) ;
    }
  }
}

foreach my$group (keys %groups) {
  foreach my$division (keys %{$groups{$group}}) {
    if (scalar(@{$groups{$group}{$division}}) < $max_length) {
      my$NAs_to_add = $max_length - scalar(@{$groups{$group}{$division}}) ;
      foreach my$i (1..$NAs_to_add) {
        push @{$groups{$group}{$division}}, "NA" ;
      }

    }
    print OUT $group, "_", $division, "\t", join("\t", @{$groups{$group}{$division}}), ")\n" ;
  }
}


close OUT ;
