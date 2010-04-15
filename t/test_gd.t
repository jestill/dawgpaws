#!/usr/bin/perl -w

# James C. Estill
# 01/22/2010
# Simple test to see if GD is working

use GD::Simple;

$img = GD::Simple->new(640, 480);
$img->fgcolor('black');
$img->bgcolor('yellow');
$img->rectangle(10, 10, 50, 50);
$img->ellipse(50, 50);

print $img->png;

exit;
