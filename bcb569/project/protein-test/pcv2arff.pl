#!/usr/bin/env perl
use strict;

my $class_label = shift() or die("Error: must provide class label\n");
my @pcv = ();
while(<>)
{
  chomp();
  my($aa, $pcvi) = split(/\t/);
  push(@pcv, $pcvi);
}
printf("%s,%s\n", join(",", @pcv), $class_label);