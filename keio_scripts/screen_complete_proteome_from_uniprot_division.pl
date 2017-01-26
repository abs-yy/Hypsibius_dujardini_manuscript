#!/usr/bin/env perl
use strict;
use warnings;
use G;

# perl screen_complete_proteome_from_uniprot_division.pl EBML_format.dat
## EMBL_format.dat ex : uniprot_sprot_archaea.dat

my $input = shift;
my %out = &get_fasta($input);

sub get_fasta{
  my $input = $_[0];
  my $tree = readFile($input, -format=>"swiss" );
  my ($dat, $div) = (split /\_/, $input)[1,2];

  $div =~ s/.dat//;

  foreach my $entry ( sort keys %{$tree} ) {
    if( defined $tree->{$entry}->{KW} && $tree->{$entry}->{KW} =~ /Complete\sproteome/ ) {
      my %fasta;
      my $seq = $tree->{$entry}->{"  "};
      $seq =~ s/\s+//g;
      say $tree->{$entry}->{LOCUS}->{id}."|".$dat."|".$div if $tree->{$entry}->{OC} =~ /Metazoa/;
      $fasta{$tree->{$entry}->{LOCUS}->{id}."|".$dat."|".$div} = $seq;
      say to_fasta(%fasta);
    }
  }
}
