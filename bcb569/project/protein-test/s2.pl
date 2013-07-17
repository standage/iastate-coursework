#!/usr/bin/env perl
use strict;

while(<>)
{
  my @values = split(/,/);
  my $total = 0;
  foreach my $p(@values[1..20])
  {
    $total += $p if($p > 0);
  }
  printf("%s\t%d\n", $values[0], $total);
}