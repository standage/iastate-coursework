#!/usr/bin/perl
use strict;
use Math::Trig;
use Math::Random qw(random_uniform);


my $sphere = sample_sphere(@ARGV);
foreach my $point(@$sphere)
{
  printf("x=%.3f, y=%.3f, z=%.3f\n", @$point);
}

sub sample_sphere
{
  my($x,$y,$z,$r,$n) = @_;
  my @zrand = random_uniform($n, $z - $r, $z + $r);
  my @phirand = random_uniform($n, 0, 2*pi);
  my @points;
  for my $i(0..$n-1)
  {
    my @point;
    
    my $theta = asin($zrand[$i] / $r);                 # theta represents latitude
    $point[0] = $r * cos($theta) * cos($phirand[$i]); # compute x coordinate of sampled point
    $point[1] = $r * cos($theta) * sin($phirand[$i]); # compute y coordinate of sampled point
    $point[2] = $zrand[$i];                           # compute z coordinate of sampled point
    
    $points[$i] = \@point;
  }
  return \@points;
}
