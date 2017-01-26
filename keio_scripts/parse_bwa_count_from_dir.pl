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
    next if $dir =~ /^\./;
    next unless $dir =~ /count$/;
#    my $sample = (split /\.fastq/, $dir)[0];
    my $sample;
    if( $dir =~ /count$/ ) {
	    $sample = (split /\./, $dir)[0];
    }
    print STDERR $dir."\n";

    push @samples, $sample;

    open my $fh, "<", $dir || die "cannot open $dir : $!";

    while( my $line = <$fh> ) {
      next if $. == 1;
      chomp $line;
      my ( $gene_long, $counts) = split /\t/, $line;

      $tree->{$gene_long}->{$sample}->{expression} = $counts;
    }
    close $fh;
}


my @samples2 = sort @samples;
print join("\t", ("Isoform",@samples2))."\n";
foreach my $isoform ( sort keys %{$tree} ) {
  my @out;
  push @out, $isoform;
  foreach my $timepoint (@samples2) {
    push @out, defined $_ ? $_ : 0  foreach $tree->{$isoform}->{$timepoint}->{expression};
  }
  print join("\t", @out)."\n";
}
