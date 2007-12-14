#!/usr/bin/perl -w
#
# J. Estill
# 10/13/2007
# Working test script to look at stacked plot output
use strict;
use warnings;

use Bio::Graphics;
use Bio::SeqFeature::Generic;

#open (IMG, ">/home/jestill/test_img.png") ||
#    die "Can not open image for output\n";

my $segment  = Bio::Graphics::Feature->new(-start=>1,-end=>700);

my $snp1     = Bio::SeqFeature::Generic ->new (-start     => 500,-end=>590,
					       -display_name =>'fred',
					       -tag=> { series => [
								   [10,20,30],
								   [30,30,0],
								   [5,45,10],
								   [5,45,10],
								   [5,45,10],
								   [50,0,50],
								   ],
							},
					       -source=>'A test',
					       );

my $snp2     = Bio::SeqFeature::Generic->new(-start     => 300,
					     -end       => 301,
					     -display_name  => 'rs12345',
					     -tag=> {
						 series => [
							    [30,20,10 ],
							    [80,10,10 ],
							    ],
						 },
					     -source=>'Another test',
					     );

my $snp3     = Bio::SeqFeature::Generic->new(-start     => 300,
					     -end       => 301,
					     -display_name  => 'rs12345',
					     -tag=> {
						 series => [
							    [ 20 ],
							    [ 10 ],
							    ],
						 },
					     -source=>'Another test',
					     );


my $panel = Bio::Graphics::Panel->new(-segment=>$segment,-width=>800);

$panel->add_track($segment,-glyph=>'arrow',-double=>1,-tick=>2);
$panel->add_track([$snp1,$snp2],
		  -height    => 50,
		  -glyph     => 'stackedplot',
		  -fixed_gap => 12,
		  -series_colors    => [qw(red blue lavender)],
		  -column_labels => [qw(a b c d e f g)],
		  -min_score => 0,
		  -max_score => 100,
		  -column_width => 8,
		  -column_font  => 'gdMediumBoldFont',
		  -scale_font   => 'gdTinyFont',
		  -label    => 1,
		  -description=>1,
		  );

$panel->add_track([$snp3],
		  -height    => 50,
		  -glyph     => 'stackedplot',
		  -fixed_gap => 12,
		  -series_colors    => [qw(gray)],
		  -column_labels => [qw(a b c d e f g)],
		  -min_score => 0,
		  -max_score => 100,
		  -column_width => 8,
		  -column_font  => 'gdMediumBoldFont',
		  -scale_font   => 'gdTinyFont',
		  -label    => 1,
		  -description=>1,
		  );

#-----------------------------+
# DRAW OUTPUT TO DISPLAY      |
#-----------------------------+
my $png_data = $panel->png;
open (DISPLAY,"| display -") 
    || die "ERROR: Can not open display for output\n";
binmode DISPLAY;
print DISPLAY $png_data;
close DISPLAY;

#print $panel->png;
#print IMG $panel->png;
#close (IMG);
