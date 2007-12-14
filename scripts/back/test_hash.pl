#!/usr/bin/perl -w
# Test of hash references

my $color = tn_color("Ames");
print "$color\n";

exit;

sub tn_color {
    my $family=shift;
    my $col_name;
    
    my %fam_color = (
		     'Ames'     => 'aliceblue',
		     'Angela'   => 'blue',
		     'Bagy'     => 'bisque',
		     'Barbara'  => 'cadetblue',
		     'Claudia'  => 'chartreuse',
		     'Daniela'  => 'chocolate',
		     'Derami'   => 'cornflowerblue',
		     'Erika'    => 'cornsilk',
		     'Fatima'   => 'red',
		     'Gujog'    => 'crimson',
		     'Hawi'     => 'cyan',
		     'Ifis'     => 'darkblue',
		     'Inga'     => 'darkcyan',
		     'Jaku'     => 'darkgoldenrod',
		     'Jeli'     => 'darkgreen',
		     'Latidu'   => 'darkkhaki',
		     'Laura'    => 'darkmagenta',
		     'Liuling'  => 'darkolive',
		     'Madil'    => 'darkorange',
		     'Morgane'  => 'darkorchid',
		     'Nusif'    => 'darkred',
		     'Retrofit' => 'darksalmon',
		     'Romani'   => 'orange',
		     'Sabrina'  => 'darkslateblue',
		     'Sutear'   => 'deeppink',
		     'TAR1'     => 'dodgerblue',
		     'unnamed4' => 'firebrick',
		     'unnamed5' => 'forestgreen',
		     'Veju'     => 'fuchsia',
		     'WHAM'     => 'gold',
		     'Wilma'    => 'goldenrod',
		     'WIS'      => 'green'
		     );
    
		$col_name = $fam_color{'Claudia'};
#    print $fam_color{'WIS'};

#		if ($family =~ "Barbara" ) {
#			$col_name = "red";
#		}
#	 	elsif ($family =~ "Jeli" ) {
#			$col_name = "blue";
#		}
#		elsif ($family =~ "Fatima" ) {
#			$col_name = "yellow";
#		}
#		elsif ($family =~ "Nusif") {
#			$col_name = "green";
#                }
#		elsif ($family =~ "Laura") {
#			$col_name = "violet";
#                }
#		else {
#			$col_name = "gray";
#		}

		# Return the color name

		return $col_name;

	
	    }
