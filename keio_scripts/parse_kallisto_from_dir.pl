#!/usr/bin/env perl
use strict;
use warnings;

use Data::Dumper;

my $input_dir = $ARGV[0];

opendir( my $DH, $input_dir);
my @dirs = readdir $DH;
closedir $DH;

my $tree;
my @samples;
foreach my $dir ( @dirs ) {
#    next if $dir =~ /^\./;
    next unless $dir =~ /.kallisto/;
    next if -f $dir ;
    my $sample = $dir;
    print STDERR $input_dir."/".$dir."/abundance.tsv\n";
    push @samples, $sample;

    open my $fh, "<", $input_dir."/".$dir."/abundance.tsv" || die "cannot open $dir/abundance.tsv : $!";

    while( my $line = <$fh> ) {
      next if $. == 1;
      chomp $line;
      my ( $gene_long, $length, $eff_length, $est_counts, $tpm ) = split /\t/, $line;
      $tree->{$gene_long}->{$sample}->{expression} = $tpm;
    }
    close $fh;
}


my @samples2 = sort  @samples;
print join("\t", ("Isoform",@samples2))."\n";
foreach my $isoform ( sort keys %{$tree} ) {
  my @out;
  push @out, $isoform;
  foreach my $timepoint (@samples2) {
    push @out, $_ foreach $tree->{$isoform}->{$timepoint}->{expression};
  }
  print join("\t", @out)."\n";
}
