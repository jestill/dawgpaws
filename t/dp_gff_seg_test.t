#!/usr/bin/perl -w
#-----------------------------------------------------------+
#                                                           |
# dp_gff_seg_test.pl - Testing gff segmentation program     |
#                                                           |
#-----------------------------------------------------------+
#  AUTHOR: James C. Estill                                  |
# CONTACT: jestill_at_sourceforge.net                       |
# STARTED: 04/15/2009                                       |
# UPDATED: 04/20/2010                                       |
#                                                           |
# SHORT DESCRIPTION:                                        |
#  Test of gff segmentation programs.                       |
#-----------------------------------------------------------+

use Test::More tests => 4;
diag ("Testing GFF segmentation programs ...");

#-----------------------------+
# GFF SIMPLE SEGMENTATION     |
#-----------------------------+
open (EXPECTED, '<data/exp/gff_seg_50x_expected.txt');
my @gffseg_exp = <EXPECTED>;
close (EXPECTED);

my @gffseg_obs = `gff_seg.pl -i data/oligo_count/gff_seg_test_in.gff --thresh 50`;

is_deeply ( \@gffseg_obs, \@gffseg_exp, "gff_seg.pl");

#-----------------------------+
# MATRIX BASED SEGMENTATION   |
#-----------------------------+

# GFF2 output
open (EXPECTED, '<data/exp/gff_segmat_def_gff2.txt');
my @gff2_gffmatrix_exp = <EXPECTED>;
close (EXPECTED);

my @gff2_gffmatrix_obs = `gff_seg_matrix.pl -i data/oligo_count/gff_seg_test_in.gff --gff-ver gff2`;

is_deeply ( \@gff2_gffmatrix_obs, \@gff2_gffmatrix_exp, "gff_seg_matrix.pl GFF2 output");

# GFF3 output
open (EXPECTED, '<data/exp/gff_segmat_def_gff3.txt');
my @gff3_gffmatrix_exp = <EXPECTED>;
close (EXPECTED);

my @gff3_gffmatrix_obs = `gff_seg_matrix.pl -i data/oligo_count/gff_seg_test_in.gff --gff-ver gff3`;

is_deeply ( \@gff3_gffmatrix_obs, \@gff3_gffmatrix_exp, "gff_seg_matrix.pl GFF3 output");


# GFF2 thresh variable
open (EXPECTED, '<data/exp/gff_segmat_newt_gff2.txt');
my @newt_gffmatrix_exp = <EXPECTED>;
close (EXPECTED);

my @newt_gffmatrix_obs = `gff_seg_matrix.pl -i data/oligo_count/gff_seg_test_in.gff --gff-ver gff2 -t 1,2,4,8,16,32,64,128,256,512`;

is_deeply ( \@newt_gffmatrix_obs, \@newt_gffmatrix_exp, "gff_seg_matrix.pl novel threshold values");
