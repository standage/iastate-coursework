#!/usr/bin/env perl
use strict;

my $aas = <>;
$aas = <>;
$aas = <>;
$aas =~ s/^\s*(.+)\s*$/$1/;
my @aa = split(/\s+/, $aas);
@aa = @aa[0..19];

printf("%s\n", join(",", @aa));
while(<>)
{
  chomp();
  $_ =~ s/^\s+(.+)\s*$/$1/;
  if( m/^\d+\s/)
  {
    my @values = split(/\s+/);
    @values = @values[2..21];
    printf("%s\n", join(",", @values));
  }
}