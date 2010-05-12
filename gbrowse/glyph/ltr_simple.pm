package Bio::Graphics::Glyph::ltr_simple;

#-----------------------------------------------------------+
# Simple LTR Retrotransposon Object
#-----------------------------------------------------------+
#
#  AUTHOR: James C. Estill
# STARTED: 03/13/2008
# UPDATED: 01/28/2010
# DESCRIPTION:
#  Glyph for a simple LTR Retrotransposon. This
#  is an object for use with GBrowse;
#
#-----------------------------------------------------------+

#-----------------------------+
# INCLUDES                    |
#-----------------------------+
use strict;
#use GD::Polygon;
# Since box is for simple features, I will want to use
# segmets since ltr retros have multiple segments
#use base 'Bio::Graphics::Glyph::segments';
use base 'Bio::Graphics::Glyph::segments';
#use base 'Bio::Graphics::Glyph::box';
#use base 'Bio::Graphics::Glyph::generic';

#-----------------------------+
# CONSTANTS                   |
#-----------------------------+
# The default TSD_Length is zero
use constant TSD_LEN => 0;

# CONTEXT FONT
# FONT CHOICES LIMITED TO
#   FONT NAME       WIDTH(pxs)    HEIGHT(pxs)
#   gdTinyFont           5            8
#   gdSmallFont          6           12
#   gdMediumBoldFont     7           12
#   gdLargeFont          8           16
#   gdGiantFont          9           15
# This assumes that we are using the GD image class
# SVG Would probably not work with the following font choice
# It would be more stable to get font_class first as is used in
# Bio::Graphics::Panel.pm. In Panel.pm this is done as
#  $gd->string($self->{key_font},$left,KEYPADTOP+$top,"KEY:",$text_color);
# However this works for me for now

# Sequence context font
my $con_font = GD::gdSmallFont;
# Sequence feature font
my $feat_font = GD::gdMediumBoldFont;

sub roundup {
    my $n = shift;
    return(($n == int($n)) ? $n : int($n + 1))
}


#-----------------------------+
# PADDING                     |
#-----------------------------+
# TSDs are used to PAD since these are set to be outside
# of the LTR Retrotransposon feature
# FROM http://www.perlmonks.org/?node_id=66948
sub tsd_len {
    my $self = shift;
    my $requested_tsd_len = $self->option('tsd_len');
    if (defined $requested_tsd_len) {
	return $requested_tsd_len;
    }
    else {
	return TSD_LEN;
    }
}

sub pad_left {

    # Left needs to be padded by the drawn length of the
    # target site duplication as well as the length of 
    # the sequence context string
    my $self = shift;
    my $pad_left = 1;
    my $req_tsd_len = $self->option('tsd_len');
    my $req_seq_context = $self->seq_con_ltr5;
    
    if ($self->tsd_len) {
	# Since this can give fractional values, I will always
	# round this up using the internal roundup function
	$pad_left = $pad_left  + ($self->tsd_len*$self->scale);
	$pad_left = roundup($pad_left);
    }
    
    if ($self->seq_con_ltr5) {
	my $char_width = $con_font->width;
	my $seq_width = length($self->seq_con_ltr5) * $char_width;
	$pad_left = $pad_left + ($seq_width);
    }

    return $pad_left;
    
}

sub pad_right {
    my $self = shift;
    my $pad_right = 1;
    my $req_tsd_len = $self->option('tsd_len');
    my $req_seq_context = $self->seq_con_ltr3;
    
    if ($self->tsd_len) {
	# Since this can give fractional values, I will always
	# round this up using the internal roundup function
	$pad_right = $pad_right  + ($self->tsd_len*$self->scale);
	$pad_right = roundup($pad_right);
    }

    if ($self->seq_con_ltr3) {
	my $char_width = $con_font->width;
	my $seq_width = length($self->seq_con_ltr3) * $char_width;
	$pad_right = $pad_right + ($seq_width);
    }

    return $pad_right;

}

#-----------------------------+
# BASE LTR RETRO SPAN         |
#-----------------------------+
# This is the bare minimum that will be drawn
sub span_fg_color {
    my $self = shift;
    my $req_span_fg_color = $self->option('span_fg_color');
    if (defined $req_span_fg_color) {
	return $self->translate_color($req_span_fg_color);
    }
    else {
	return $self->fgcolor;
    }
}

sub span_bg_color {
    my $self = shift;
    my $req_span_bg_color = $self->option('span_bg_color');
    if (defined $req_span_bg_color) {
	return $self->translate_color($req_span_bg_color);
	
    }
    else {
	return $self->bgcolor;
    }
}

#-----------------------------+
# LABEL FEATURES              |
#-----------------------------+
sub label_feat {
    my $self = shift;
    my $req_label_feat = $self->option('label_feat');
    if (defined $req_label_feat) {
	return $req_label_feat;
    }
    else {
	return 0;
    }
}

#-----------------------------+
# SEQUENCE CONTEXT            |
#-----------------------------+
# The purpose of printing the sequence context is to show the putative TSDs
# on either side of the LTR. This could also be used to draw a larger context
# sequence such as that returned by LTR_Struc
sub seq_con_ltr5 {
    my $self = shift;
    my $req_seq_con_ltr5 = $self->option('seq_con_ltr5');
    if (defined $req_seq_con_ltr5) {
	return $req_seq_con_ltr5;
    }
    else {
	return 0;
    }
}

sub seq_con_ltr3 {
    my $self = shift;
    my $req_seq_con_ltr3 = $self->option('seq_con_ltr3');
    if (defined $req_seq_con_ltr3) {
	return $req_seq_con_ltr3;
    }
    else {
	return 0;
    }
}

sub seq_con_color {
    my $self = shift;
    my $req_seq_con_color = $self->option('seq_con_color');
    if (defined $req_seq_con_color) {
	return $self->translate_color($req_seq_con_color);
    }
    else {
	return $self->fgcolor;
    }
}

#-----------------------------+
# TARGET SITE DUPLICATION     |
#-----------------------------+
sub tsd_fg_color {
    my $self = shift;
    my $req_tsd_fg_color = $self->option('tsd_fg_color');
    if (defined $req_tsd_fg_color) {
	return $self->translate_color($req_tsd_fg_color);
    }
    else {
	return $self->fgcolor;
    }
}

sub tsd_bg_color {
    my $self = shift;
    my $req_tsd_bg_color = $self->option('tsd_bg_color');
    if (defined $req_tsd_bg_color) {
	return $self->translate_color($req_tsd_bg_color);
	
    }
    else {
	return $self->bgcolor;
    }
}


#-----------------------------+
# GENERAL TG/CA COLORS        |
#-----------------------------+
# Commented out 01/26/2010
#sub tg_bg_color {
#    my $self = shift;
#    my $req_tg_bg_color = $self->option('tg_bg_color');
#    if (defined $req_tg_bg_color) {
#	return $self->translate_color($req_tg_bg_color);
#    }
#    else {
#	return $self->bgcolor;
#    }
#}
#
#sub ca_bg_color {
#    my $self = shift;
#    my $req_ca_bg_color = $self->option('ca_bg_color');
#    if (defined $req_ca_bg_color) {
#	return $self->translate_color($req_ca_bg_color);
#    }
#    else {
#	return $self->bgcolor;
#    }
#}

#-----------------------------+
# 5' LONG TERMINAL REPEAT     |
#-----------------------------+
sub ltr_fg_color {
    my $self = shift;
    my $req_ltr_fg_color = $self->option('ltr_fg_color');
    if (defined $req_ltr_fg_color) {
	return $self->translate_color($req_ltr_fg_color);
    }
    else {
	return $self->fgcolor;
    }
}

sub ltr_bg_color {
    my $self = shift;
    my $req_ltr_bg_color = $self->option('ltr_bg_color');
    if (defined $req_ltr_bg_color) {
	return $self->translate_color($req_ltr_bg_color);
    }
    else {
	return $self->bgcolor;
    }
}


sub has_ltr5_tg {
    my $self = shift;
    my $req_has_ltr5_tg = $self->option('has_ltr5_tg');
    if (defined $req_has_ltr5_tg) {
	return $req_has_ltr5_tg;
    }
    else {
	return 0;
    }
}

sub has_ltr5_ca {
    my $self = shift;
    my $req_has_ltr5_ca = $self->option('has_ltr5_ca');
    if (defined $req_has_ltr5_ca) {
	return $req_has_ltr5_ca;
    }
    else {
	return 0;
    }
}

#-----------------------------+
# LTR TG AND CA               |
#-----------------------------+
#sub has_ltr3_tg {
#    my $self = shift;
#    my $req_has_ltr3_tg = $self->option('has_ltr3_tg');
#    if (defined $req_has_ltr3_tg) {
#	return $req_has_ltr3_tg;
#    }
#    else {
#	return 0;
#    }
#}
#
#sub has_ltr3_ca {
#    my $self = shift;
#    my $req_has_ltr3_ca = $self->option('has_ltr3_ca');
#    if (defined $req_has_ltr3_ca) {
#	return $req_has_ltr3_ca;
#    }
#    else {
#	return 0;
#    }
#}

#-----------------------------+
# LTR ARROW                   |
#-----------------------------+
sub ltr_arrow_fg_color {
    my $self = shift;
    my $req_ltr_arrow_fg_color = $self->option('ltr_arrow_fg_color');
    if (defined $req_ltr_arrow_fg_color) {
	return $self->translate_color($req_ltr_arrow_fg_color);
    }
    else {
	return $self->fgcolor;
    }
}

sub ltr_arrow_bg_color {
    my $self = shift;
    my $req_ltr_arrow_bg_color = $self->option('ltr_arrow_bg_color');
    if (defined $req_ltr_arrow_bg_color) {
	return $self->translate_color($req_ltr_arrow_bg_color);
    }
    else {
	return $self->bgcolor;
    }
}

# Boolean to draw the ltr with and arrow
sub ltr_arrow {
    my $self = shift;
    my $req_ltr_arrow = $self->option('ltr_arrow');
    if (defined $req_ltr_arrow) {
	return $req_ltr_arrow;
    }
    else {
	return 1;
    }
}

#-----------------------------+
# PRIMER BINDING SITE         |
#-----------------------------+
sub pbs_fg_color {
    my $self = shift;
    my $req_pbs_fg_color = $self->option('pbs_fg_color');
    if (defined $req_pbs_fg_color) {
	return $self->translate_color($req_pbs_fg_color);
    }
    else {
	return $self->fgcolor;
    }
}

sub pbs_bg_color {
    my $self = shift;
    my $req_pbs_bg_color = $self->option('pbs_bg_color');
    if (defined $req_pbs_bg_color) {
	return $self->translate_color($req_pbs_bg_color);
    }
    else {
	return $self->bgcolor;
    }
}

#-----------------------------+
# POLYPURINE TRACT            |
#-----------------------------+
sub rr_fg_color {
    my $self = shift;
    my $req_rr_fg_color = $self->option('rr_fg_color');
    if (defined $req_rr_fg_color) {
	return $self->translate_color($req_rr_fg_color);
    }
    else {
	return $self->fgcolor;
    }
}

sub rr_bg_color {
    my $self = shift;
    my $req_rr_bg_color = $self->option('rr_bg_color');
    if (defined $req_rr_bg_color) {
	return $self->translate_color($req_rr_bg_color);
    }
    else {
	return $self->bgcolor;
    }
}

#-----------------------------+
# MATURE PROTEIN COLORS       |
#-----------------------------+
sub protein_fg_color {
    my $self = shift;
    my $req_protein_fg_color = $self->option('protein_fg_color');
    if (defined $req_protein_fg_color) {
	return $self->translate_color($req_protein_fg_color);
    }
    else {
	return $self->fgcolor;
    }
}

sub protein_bg_color {
    my $self = shift;
    my $req_protein_bg_color = $self->option('protein_bg_color');
    if (defined $req_protein_bg_color) {
	return $self->translate_color($req_protein_bg_color);
    }
    else {
	return $self->bgcolor;
    }
}

#-----------------------------+
# ENVELOPE                    |
#-----------------------------+
#sub env_fg_color {
#    my $self = shift;
#    my $req_env_fg_color = $self->option('env_fg_color');
#    if (defined $req_env_fg_color) {
#	return $self->translate_color($req_env_fg_color);
#    }
#    else {
#	return $self->fgcolor;
#    }
#}
#
#sub env_bg_color {
#    my $self = shift;
#    my $req_env_bg_color = $self->option('env_bg_color');
#    if (defined $req_env_bg_color) {
#	return $self->translate_color($req_env_bg_color);
#    }
#    else {
#	return $self->bgcolor;
#    }
#}
#

#-----------------------------+
# GAG                         |
#-----------------------------+
sub gag_fg_color {
    my $self = shift;
    my $req_gag_fg_color = $self->option('gag_fg_color');
    if (defined $req_gag_fg_color) {
	return $self->translate_color($req_gag_fg_color);
    }
    else {
	return $self->fgcolor;
    }
}

sub gag_bg_color {
    my $self = shift;
    my $req_gag_bg_color = $self->option('gag_bg_color');
    if (defined $req_gag_bg_color) {
	return $self->translate_color($req_gag_bg_color);
    }
    else {
	return $self->bgcolor;
    }
}

#-----------------------------+
# Zinc Finger CCHC            |
#-----------------------------+
sub zf_fg_color {
    my $self = shift;
    my $req_zf_fg_color = $self->option('zf_fg_color');
    if (defined $req_zf_fg_color) {
	return $self->translate_color($req_zf_fg_color);
    }
    else {
	return $self->fgcolor;
    }
}

sub zf_bg_color {
    my $self = shift;
    my $req_zf_bg_color = $self->option('zf_bg_color');
    if (defined $req_zf_bg_color) {
	return $self->translate_color($req_zf_bg_color);
    }
    else {
	return $self->bgcolor;
    }
}

#-----------------------------+
# POL                         |
#-----------------------------+
sub pol_fg_color {
    my $self = shift;
    my $req_pol_fg_color = $self->option('pol_fg_color');
    if (defined $req_pol_fg_color) {
	return $self->translate_color($req_pol_fg_color);
    }
    else {
	return $self->fgcolor;
    }
}

sub pol_bg_color {
    my $self = shift;
    my $req_pol_bg_color = $self->option('pol_bg_color');
    if (defined $req_pol_bg_color) {
	return $self->translate_color($req_pol_bg_color);
    }
    else {
	return $self->bgcolor;
    }
}

#-----------------------------+
# PROTEASE                    |
#-----------------------------+
#sub pro_fg_color {
#    my $self = shift;
#    my $req_pro_fg_color = $self->option('pro_fg_color');
#    if (defined $req_pro_fg_color) {
#	return $self->translate_color($req_pro_fg_color);
#    }
#    else {
#	return $self->fgcolor;
#    }
#}
#
#sub pro_bg_color {
#    my $self = shift;
#    my $req_pro_bg_color = $self->option('pro_bg_color');
#    if (defined $req_pro_bg_color) {
#	return $self->translate_color($req_pro_bg_color);
#    }
#    else {
#	return $self->bgcolor;
#    }
#}
#
#-----------------------------+
# INTEGRASE                   |
#-----------------------------+
#sub int_fg_color {
#    my $self = shift;
#    my $req_int_fg_color = $self->option('int_fg_color');
#    if (defined $req_int_fg_color) {
#	return $self->translate_color($req_int_fg_color);
#    }
#    else {
#	return $self->fgcolor;
#    }
#}
#
#sub int_bg_color {
#    my $self = shift;
#    my $req_int_bg_color = $self->option('int_bg_color');
#    if (defined $req_int_bg_color) {
#	return $self->translate_color($req_int_bg_color);
#    }
#    else {
#	return $self->bgcolor;
#    }
#}

#-----------------------------+
# REVERSE TRANSCRIPTASE       |
#-----------------------------+
# sub rt_fg_color {
#     my $self = shift;
#     my $req_rt_fg_color = $self->option('rt_fg_color');
#     if (defined $req_rt_fg_color) {
# 	return $self->translate_color($req_rt_fg_color);
#     }
#     else {
# 	return $self->fgcolor;
#     }
# }

# sub rt_bg_color {
#     my $self = shift;
#     my $req_rt_bg_color = $self->option('rt_bg_color');
#     if (defined $req_rt_bg_color) {
# 	return $self->translate_color($req_rt_bg_color);
#     }
#     else {
# 	return $self->bgcolor;
#     }
# }

# #-----------------------------+
# # CHRromatin Organiation      |
# # MOdifier Domain             |
# #-----------------------------+
# sub chromo_fg_color {
#     my $self = shift;
#     my $req_chromo_fg_color = $self->option('chromo_fg_color');
#     if (defined $req_chromo_fg_color) {
# 	return $self->translate_color($req_chromo_fg_color);
#     }
#     else {
# 	return $self->fgcolor;
#     }
# }

# sub chromo_bg_color {
#     my $self = shift;
#     my $req_chromo_bg_color = $self->option('chromo_bg_color');
#     if (defined $req_chromo_bg_color) {
# 	return $self->translate_color($req_chromo_bg_color);
#     }
#     else {
# 	return $self->bgcolor;
#     }
# }

#-----------------------------+
# RNASEH                      |
#-----------------------------+
sub rh_fg_color {
    my $self = shift;
    my $req_rh_fg_color = $self->option('rh_fg_color');
    if (defined $req_rh_fg_color) {
 	return $self->translate_color($req_rh_fg_color);
    }
    else {
 	return $self->fgcolor;
    }
}

sub rh_bg_color {
    my $self = shift;
    my $req_rh_bg_color = $self->option('rh_bg_color');
    if (defined $req_rh_bg_color) {
 	return $self->translate_color($req_rh_bg_color);
    }
    else {
 	return $self->bgcolor;
    }
}

#-----------------------------+
# DRAW THE LTR RETRO          |
#-----------------------------+
sub draw_component {

    my $self = shift;
    my ($gd,$dx,$dy) = @_;

    # test self level
    #-----------------------------+
    # PRINT TO THE DRAWING AREA   |
    #-----------------------------+
    #my $level = $self->level;

    # Get the bounding box of the glyph
    my ($left,$top,$right,$bottom) = $self->bounds($dx,$dy);
    
    # Get the height
    my $height = $bottom - $top + 1;
    
    # Get the colors
    my $back_color = $self->bgcolor;
    my $fore_color = $self->fgcolor;

    #-----------------------------+
    # BASE LTR RETRO SPAN         |
    #-----------------------------+
    # This used the defauult background and foreground
    # color, it may be useful to make this a varaible
    # outside of the base.
    # Using with gbrowse the followings is always drawn
    # and the basic span is not being drawn
 #   if ( $level = 0 ) {
#	print STDOUT
  #  }


    #-----------------------------+
    # LTR Retrotransposon Span    |
    #-----------------------------+
    # The following is not working
    if ( $self->feature->type =~ /LTR_retrotransposon/i ) {

#	my $poly = GD::Polygon->new;
#	my $span_top = $top+($height * 0.25);
#	my $span_bottom = $bottom - ($height * 0.25);
#	$poly->addPt($left,$top);
#	$poly->addPt($right,$top);
#	$poly->addPt($right,$bottom);
#	$poly->addPt($left,$bottom);
#	$poly->addPt($left,$top);
#	$gd->filledPolygon($poly,$self->ltr_bg_color);
#	$gd->polygon($poly,$self->ltr_fg_color);

	$self->filled_box($gd,
			  $left,
			  $top+($height * 0.05 ),
			  $right,
			  $bottom-($height * 0.05 ),
			  $self->ltr_bg_color,
			  $self->ltr_fg_color);

    }


    #-----------------------------+
    # LTR                         |
    #-----------------------------+
    if ($self->feature->type =~ /five_prime_LTR||three_prime_LTR||solo_LTR/i) {

	$self->filled_box($gd,
			  $left,
			  $top+($height * 0.05 ),
			  $right,
			  $bottom-($height * 0.05 ),
			  $self->ltr_bg_color,
			  $self->ltr_fg_color);
	
	# Top strand Arrow
	my $ltr5_tri = GD::Polygon->new;
	$ltr5_tri->addPt($left, $top+($height * 0.05 ));     # PT A
	$ltr5_tri->addPt($right, $top+($height * 0.5) );     # PT B
	$ltr5_tri->addPt($left, $bottom-($height * 0.05 ));  # PT C
	$ltr5_tri->addPt($left, $top+($height * 0.05  ));    # PT A
	$gd->filledPolygon($ltr5_tri,$self->ltr_arrow_bg_color);
	$gd->polygon($ltr5_tri,$self->ltr_arrow_fg_color);


    };


    #-----------------------------+
    # PRIMER BINDING SITE         |
    #-----------------------------+
    if ( $self->feature->type =~ /primer_binding_site/i ) {

	$self->filled_box($gd,
			  $left,
			  $top+($height * 0.05 ),
			  $right,
			  $bottom-($height * 0.05 ),
			  $self->pbs_bg_color,
			  $self->pbs_fg_color);

    }

    #-----------------------------+
    # POLYPURINE TRACT            |
    #-----------------------------+
    if ( $self->feature->type =~ /RR_tract/i ) {
	
	$self->filled_box($gd,
			  $left,
			  $top+($height * 0.05 ),
			  $right,
			  $bottom-($height * 0.05 ),
			  $self->pbs_bg_color,
			  $self->pbs_fg_color);

    }

    #-----------------------------+
    # Target Site Duplication     |
    #-----------------------------+
    if ( $self->feature->type =~ /target_site_duplication/i ) {
	
	$self->filled_box($gd,
			  $left,
			  $top+($height * 0.05 ),
			  $right,
			  $bottom-($height * 0.05 ),
			  $self->tsd_bg_color,
			  $self->tsd_fg_color);

    }

    #-----------------------------+
    # MATURE PROTEIN REGION       |
    #-----------------------------+
    # Trying to allow for printing of mutlple protein types
    if ( $self->feature->type =~ /mature_protein_region/i ) {
	
	my $prot_bg_color;
	my $prot_fg_color;
	
	$self->filled_box($gd,
			  $left,
			  $top+($height * 0.05 ),
			  $right,
			  $bottom-($height * 0.05 ),
			  $self->protein_bg_color,
			  $self->protein_fg_color);

	# Label features if desired
	# TO DO: Add check that there is enough space
	if ($self->label_feat) {
	    my $height = $bottom - $top;
	    my $name = $self->feature->name;
	    my $char_height = $feat_font->height;
	    my $word_width = length($name) * $feat_font->width;
	    my $part_width = $right - $left;
 	    my $feat_start_x = $left + 
 		( ( ($right - $left)/2 ) - 
 		  ( (($feat_font->width * length($name))/2) 
 		  )
 		) ;
	    my $feat_start_y =  $top + ($height * 0.5) - ($char_height/2);

	    # Print label if there is room
	    if ($word_width < $part_width) {
		$gd->string($feat_font, 
			    $feat_start_x,
			    $feat_start_y,
			    $name,
			    $self->seq_con_color);
	    }

	}

    }

    #-----------------------------+
    # TRANSPOSABLE ELEMENT PROTEIN|
    #-----------------------------+
    # SO:0000111
    # The current SO for TE proteins is not well defined
    # Trying to allow for printing of mutlple protein types
    # A simple copy and paste of mature_protein_region for now
    if ( $self->feature->type =~ /transposable_element_protein/i ) {
	
	my $prot_bg_color;
	my $prot_fg_color;
	
	$self->filled_box($gd,
			  $left,
			  $top+($height * 0.05 ),
			  $right,
			  $bottom-($height * 0.05 ),
			  $self->protein_bg_color,
			  $self->protein_fg_color);

	# Label features if desired
	# TO DO: Add check that there is enough space
	if ($self->label_feat) {
	    my $height = $bottom - $top;
	    my $name = $self->feature->name;
	    my $char_height = $feat_font->height;
	    my $word_width = length($name) * $feat_font->width;
	    my $part_width = $right - $left;
 	    my $feat_start_x = $left + 
 		( ( ($right - $left)/2 ) - 
 		  ( (($feat_font->width * length($name))/2) 
 		  )
 		) ;
	    my $feat_start_y =  $top + ($height * 0.5) - ($char_height/2);

	    # Print label if there is room
	    if ($word_width < $part_width) {
		$gd->string($feat_font, 
			    $feat_start_x,
			    $feat_start_y,
			    $name,
			    $self->seq_con_color);
	    }

	}

    }


    #-----------------------------+
    # POLYPEPTIDE DOMAIN          |
    #-----------------------------+
    if ( $self->feature->type =~ /polypeptide_domain/i ) {
	
	my $prot_bg_color;
	my $prot_fg_color;
	
	$self->filled_box($gd,
			  $left,
			  $top+($height * 0.05 ),
			  $right,
			  $bottom-($height * 0.05 ),
			  $self->protein_bg_color,
			  $self->protein_fg_color);

	# Label features if desired
	# TO DO: Add check that there is enough space
	if ($self->label_feat) {
	    my $height = $bottom - $top;
	    my $name = $self->feature->name;
	    my $char_height = $feat_font->height;
	    my $word_width = length($name) * $feat_font->width;
	    my $part_width = $right - $left;
 	    my $feat_start_x = $left + 
 		( ( ($right - $left)/2 ) - 
 		  ( (($feat_font->width * length($name))/2) 
 		  )
 		) ;
	    my $feat_start_y =  $top + ($height * 0.5) - ($char_height/2);

	    # Print label if there is room
	    if ($word_width < $part_width) {
		$gd->string($feat_font, 
			    $feat_start_x,
			    $feat_start_y,
			    $name,
			    $self->seq_con_color);
	    }

	}

    }


    #-----------------------------+
    # 5' SEQUENCE CONTEXT         |
    #-----------------------------+
#    if ($self->seq_con_ltr5) {
#
#	my $char_width = $con_font->width;
#	my $char_height = $con_font->height;
#	my $seq_width = length($self->seq_con_ltr5) * $char_width;
#
#	my $seq_val = $self->seq_con_ltr5;
#	my $con_start_y =  $top + ($height * 0.5) - ($char_height/2);
#
#	my $con_start_x = $left - $seq_width - 1;
#
#	if ($self->tsd_len) {
#	    my $offset = $self->tsd_len*$self->scale;
#	    $offset = roundup($offset);
#	    $con_start_x = $con_start_x - $offset;
#	}
#	
#	$gd->string($con_font, 
#		    $con_start_x, 
#		    $con_start_y, 
#		    $seq_val, 
#		    $self->seq_con_color);
#    }
#
#    #-----------------------------+
#    # 3' SEQUENCE CONTEXT         |
#    #-----------------------------+
#    if ($self->seq_con_ltr3) {
#	
#	my $char_width = $con_font->width;
#	my $char_height = $con_font->height;
#	my $seq_width = length($self->seq_con_ltr3) * $char_width;
#
#	# Need to switch the following to not use pad_right since pad_right
#	# will need to include the text string
#	my $con_start_y =  $top + ($height * 0.5) - ($char_height/2);
#	my $seq_val = $self->seq_con_ltr3;
#
#	my $con_start_x = $right + 1;
#	
#	if ($self->tsd_len) {
#	    my $offset = $self->tsd_len*$self->scale;
#	    $offset = roundup($offset);
#	    $con_start_x = $con_start_x + $offset;
#	}
#
#	$gd->string($con_font, 
#		    $con_start_x, 
#		    $con_start_y, 
#		    $seq_val, 
#		    $self->seq_con_color);
#    }

}
 
1;

#-----------------------------------------------------------+
# HISTORY
#-----------------------------------------------------------+
# 03/13/2008 -- 03/14/2008
#  - Main object started, includes
#    - LTR by start stopewith triangle
#    - LTR by Length
#    - TSD
#    - GAG
#    - Pol
#    - Pro
#    - INT
#    - RH
#    - RT
#    - PBS
#    - PPT (using RR to follow song complience)
#  - Single object that only would be interpreted correctly
#    on the plus strand for a full length LTR retro
#
# 03/14/2008
# - Added boolean for drawing LTR arrows, default is true 1
# - Added half circles for the tg/ca based at either end of 
#   the LTR
# - Adding color features as variables
#    
# 03/16/2008
# - Added color features as variables for 
#    Pol, Gag, RR, PBS, INT, RH, PRO, LTR Arrows
# - Started genomic context 
#
# 03/17/2008
# - Continuing work on genomic context
#    This buffers sequences by one pixel on either side to 
#    avoid overlap with drawn objects, 
# - Modifed left and right padding to take into account all
#   sequence features that can exist outside of the drawn 
#   object. This includes a 1 pixel offset at all times.
# - Added tgca foreground and background colors as variables
# - Cleaned up code, removed STDERR print statments
# - Added span bgcolor and fg color
# - Added Lablels for Gag, Pol, AP, INT, RT, RH
# 03/26/2008
# - Added chromo domain with
#    -chromo_start,chromo_end,chromo_bg_color, chromo_fg_color
#    -This domain is about 60 aa so will not be labeled for now
#
# 03/28/2008
# - Adding env domain
#
# 04/18/2008
# - Working on dealing with value 'N' being passed as a coordinate
#   I will treat this as a null/N value and not draw the feature
#   by returned 0/false when fetching feature locations
# - Added labels to zinc finger cchc, env, and chromo domains 
# 
# 01/28/2010
# - Rewritten to work with SO complient GFF3 in Gbrowse
