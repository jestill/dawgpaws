#!/usr/bin/perl -w
#-----------------------------------------------------------+
#                                                           |
# dp_cnv_test.pl - Test DAWGPAWS conversion programs.       |
#                                                           |
#-----------------------------------------------------------+
#  AUTHOR: James C. Estill                                  |
# CONTACT: jestill_at_sourceforge.net                       |
# STARTED: 04/14/2009                                       |
# UPDATED: 04/20/2010                                       |
#                                                           |
# SHORT DESCRIPTION:                                        |
#  Test of the DAWGPAWS conversion programs.                |
#-----------------------------------------------------------+

use Test::More tests => 20;

diag ("Testing DAWGPAWS conversion programs ...");

#-----------------------------------------------------------+
# FGENESH                                                   |
#-----------------------------------------------------------+

#-----------------------------+
# FGENESH TXT/ GFF2 OUT       |
#-----------------------------+
open (EXPECTED, '<data/exp/fgenesh_cnv_expected.txt');
my @gff2_fgenesh_exp = <EXPECTED>;
close (EXPECTED);

my @gff2_fgenesh_obs =  `cnv_fgenesh2gff.pl -i data/HEX3045G05/fgenesh/fgenesh.txt -gff-ver gff2`;

# Test the observed fgenesh parser result against the expected result
is_deeply ( \@gff2_fgenesh_obs, \@gff2_fgenesh_exp, 
	    "cnv_fgenesh2gff.pl text conversion GFF2 output" );

#-----------------------------+
# FGENESH HTML/GFF2 OUT       |
#-----------------------------+
open (EXPECTED, '<data/exp/fgenesh_cnv_expected.txt');
my @gff2_fgenesh_html_exp = <EXPECTED>;
close (EXPECTED);

# GFF2 Test
my @gff2_fgenesh_html_obs =  `cnv_fgenesh2gff.pl -i data/HEX3045G05/fgenesh/fgenesh_html.txt --html --gff-ver gff2`;

# Test the observed fgenesh parser result against the expected result
is_deeply ( \@gff2_fgenesh_html_obs, \@gff2_fgenesh_html_exp, 
	    "cnv_fgenesh2gff.pl html conversion GFF2 output" );

#-----------------------------+
# FGENESH TXT/ GFF3 OUT       |
#-----------------------------+
open (EXPECTED, '<data/exp/fgenesh_cnv_expected_gff3.txt');
my @gff3_fgenesh_exp = <EXPECTED>;
close (EXPECTED);

my @gff3_fgenesh_obs =  `cnv_fgenesh2gff.pl -i data/HEX3045G05/fgenesh/fgenesh.txt -gff-ver gff3`;

# Test the observed fgenesh parser result against the expected result
is_deeply ( \@gff3_fgenesh_obs, \@gff3_fgenesh_exp, 
	    "cnv_fgenesh2gff.pl text conversion GFF3 output" );

#-----------------------------+
# FGENESH HTML/GFF3 OUT       |
#-----------------------------+
open (EXPECTED, '<data/exp/fgenesh_cnv_expected_gff3.txt');
my @gff3_fgenesh_html_exp = <EXPECTED>;
close (EXPECTED);

# GFF2 Test
my @gff3_fgenesh_html_obs =  `cnv_fgenesh2gff.pl -i data/HEX3045G05/fgenesh/fgenesh_html.txt --html --gff-ver gff3`;

# Test the observed fgenesh parser result against the expected result
is_deeply ( \@gff3_fgenesh_html_obs, \@gff3_fgenesh_html_exp, 
	    "cnv_fgenesh2gff.pl html conversion GFF3" );

#-----------------------------------------------------------+
# GENEMARK                                                  |
#-----------------------------------------------------------+

#-----------------------------+
# GFF2 OUTPUT TEST            |
#-----------------------------+
open (EXPECTED, '<data/exp/genemark_cnv_expected_gff2.txt');
my @gff2_genemark_exp = <EXPECTED>;
close (EXPECTED);
my @gff2_genemark_obs = `cnv_genemark2gff.pl -i data/HEX3045G05/genemark/HEX3045G05_genemark_ta.out --gff-ver gff2 --seqname HEX3045G05`;

is_deeply ( \@gff2_genemark_obs, \@gff2_genemark_exp, "cnv_genemark2gff.pl GFF2 output");

#-----------------------------+
# GFF2 OUTPUT TEST            |
#-----------------------------+
open (EXPECTED, '<data/exp/genemark_cnv_expected_gff3.txt');
my @gff3_genemark_exp = <EXPECTED>;
close (EXPECTED);
my @gff3_genemark_obs = `cnv_genemark2gff.pl -i data/HEX3045G05/genemark/HEX3045G05_genemark_ta.out --gff-ver gff3 --seqname HEX3045G05`;

is_deeply ( \@gff2_genemark_obs, \@gff2_genemark_exp, "cnv_genemark2gff.pl GFF3 output");

#-----------------------------------------------------------+
# FIND LTR                                                  |
#-----------------------------------------------------------+

#-----------------------------+
# GFF2 OUTPUT TEST            |
#-----------------------------+
open (EXPECTED, '<data/exp/findltr_cnv_expected_gff2.txt');
my @gff2_findltr_exp = <EXPECTED>;
close (EXPECTED);
my @gff2_findltr_obs = `cnv_findltr2gff.pl -i data/HEX3045G05/find_ltr/HEX3045G05_findltr_def.txt --gff-ver gff2`;

is_deeply ( \@gff2_findltr_obs, \@gff2_findltr_exp, "cnv_findltr2gff.pl GFF2 output");

#-----------------------------+
# GFF3 OUTPUT TEST            |
#-----------------------------+
open (EXPECTED, '<data/exp/findltr_cnv_expected_gff3.txt');
my @gff3_findltr_exp = <EXPECTED>;
close (EXPECTED);
my @gff3_findltr_obs = `cnv_findltr2gff.pl -i data/HEX3045G05/find_ltr/HEX3045G05_findltr_def.txt --gff-ver gff3`;

is_deeply ( \@gff3_findltr_obs, \@gff3_findltr_exp, "cnv_findltr2gff.pl GFF3 output");

#-----------------------------------------------------------+
# LTR FINDER                                                |
#-----------------------------------------------------------+

#-----------------------------+
# GFF2 OUTPUT TEST            |
#-----------------------------+
open (EXPECTED, '<data/exp/ltrfinder_cnv_exp_gff2.txt');
my @gff2_ltrfinder_exp = <EXPECTED>;
close (EXPECTED);
my @gff2_ltrfinder_obs = `cnv_ltrfinder2gff.pl -i data/HEX3045G05/ltr_finder/HEX3045G05_def_ltr_finder.txt --gff-ver gff2`;

is_deeply ( \@gff2_ltrfinder_obs, \@gff2_ltrfinder_exp, "cnv_ltrfinder2gff.pl GFF2 output");

#-----------------------------+
# GFF3 OUTPUT TEST            |
#-----------------------------+

open (EXPECTED, '<data/exp/ltrfinder_cnv_exp_gff3.txt');
my @gff3_ltrfinder_exp = <EXPECTED>;
close (EXPECTED);
my @gff3_ltrfinder_obs = `cnv_ltrfinder2gff.pl -i data/HEX3045G05/ltr_finder/HEX3045G05_def_ltr_finder.txt --gff-ver gff3`;

is_deeply ( \@gff3_ltrfinder_obs, \@gff3_ltrfinder_exp, "cnv_ltrfinder2gff.pl GFF3 output");

#-----------------------------------------------------------+
# LTR_seq CONVERSION                                        |
#-----------------------------------------------------------+

#-----------------------------+
# GFF2 OUTPUT TEST            |
#-----------------------------+
open (EXPECTED, '<data/exp/ltrseq_cnv_expected_gff2.txt');
my @gff2_ltrseq_exp = <EXPECTED>;
close (EXPECTED);
my @gff2_ltrseq_obs = `cnv_ltrseq2gff.pl -i data/HEX3045G05/ltr_seq/HEX3045G05_ltrseq.txt -s HEX3045G05 --gff-ver gff2`;

is_deeply ( \@gff2_ltrseq_obs, \@gff2_ltrseq_exp, "cnv_ltrseq2gff.pl GFF2 output" );


#-----------------------------+
# GFF3 OUTPUT TEST            |
#-----------------------------+
open (EXPECTED, '<data/exp/ltrseq_cnv_expected_gff3.txt');
my @gff3_ltrseq_exp = <EXPECTED>;
close (EXPECTED);

# GFF2 Test
my @gff3_ltrseq_obs = `cnv_ltrseq2gff.pl -i data/HEX3045G05/ltr_seq/HEX3045G05_ltrseq.txt -s HEX3045G05 --gff-ver gff3`;

is_deeply ( \@gff3_ltrseq_obs, \@gff3_ltrseq_exp, "cnv_ltrseq2gff.pl GFF3 output" );

#-----------------------------------------------------------+
# REPSEEK CONVERSION                                        |
#-----------------------------------------------------------+

#-----------------------------+
# GFF2 OUTPUT TEST            |
#-----------------------------+
open (EXPECTED, '<data/exp/repseek_cnv_expected_gff2.txt');
my @gff2_repseek_exp = <EXPECTED>;
close (EXPECTED);

my @gff2_repseek_obs = `cnv_repseek2gff.pl -i data/HEX3045G05/repseek/repseek_l20.txt --gff-ver gff2`;

is_deeply ( \@gff2_repseek_obs, \@gff2_repseek_exp, "cnv_repseek2gff.pl GFF2 output");

#-----------------------------+
# GFF3 OUTPUT TEST            |
#-----------------------------+
open (EXPECTED, '<data/exp/repseek_cnv_expected_gff3.txt');
my @gff3_repseek_exp = <EXPECTED>;
close (EXPECTED);

my @gff3_repseek_obs = `cnv_repseek2gff.pl -i data/HEX3045G05/repseek/repseek_l20.txt --gff-ver gff3`;

is_deeply ( \@gff3_repseek_obs, \@gff3_repseek_exp, "cnv_repseek2gff.pl GFF3 output");


#-----------------------------------------------------------+
# BLAST CONVERSION                                          |
#-----------------------------------------------------------+

#-----------------------------+
# GFF2 OUTPUT TEST            |
#-----------------------------+

open (EXPECTED, '<data/exp/blast_cnv_AcTA_gff2.txt');
my @gff2_blast_exp = <EXPECTED>;
close (EXPECTED);
my @gff2_blast_obs = `cnv_blast2gff.pl -i data/HEX3045G05/blast/HEX3045G05_AcTA_1.bln --gff-ver gff2`;
is_deeply ( \@gff2_blast_obs, \@gff2_blast_exp, "cnv_blast2gff.pl with GFF2 output" );

#-----------------------------+
# GFF3 OUTPUT TEST            |
#-----------------------------+

open (EXPECTED, '<data/exp/blast_cnv_AcTA_gff3.txt');
my @gff3_blast_exp = <EXPECTED>;
close (EXPECTED);
my @gff3_blast_obs = `cnv_blast2gff.pl -i data/HEX3045G05/blast/HEX3045G05_AcTA_1.bln --gff-ver gff3`;
is_deeply ( \@gff3_blast_obs, \@gff3_blast_exp, "cnv_blast2gff.pl with GFF3 output" );

#-----------------------------------------------------------+
# REPEATMASKER CONVERSION                                   |
#-----------------------------------------------------------+

#-----------------------------+
# GFF2 OUTPUT TEST            |
#-----------------------------+

open (EXPECTED, '<data/exp/repeatmasker_cnv_expected_gff2.txt');
my @gff2_repmask_exp = <EXPECTED>;
close (EXPECTED);
my @gff2_repmask_obs = `cnv_repmask2gff.pl -i data/HEX3045G05/rm/HEX3045G05_TREP9.rm.out --gff-ver gff2`;
is_deeply ( \@gff2_repmask_obs, \@gff2_repmask_exp, "cnv_repmask2gff.pl with GFF2 output" );

#-----------------------------+
# GFF3 OUTPUT TEST            |
#-----------------------------+
open (EXPECTED, '<data/exp/repeatmasker_cnv_expected_gff3.txt');
my @gff3_repmask_exp = <EXPECTED>;
close (EXPECTED);
my @gff3_repmask_obs = `cnv_repmask2gff.pl -i data/HEX3045G05/rm/HEX3045G05_TREP9.rm.out --gff-ver gff3`;

is_deeply ( \@gff3_repmask_obs, \@gff3_repmask_exp, "cnv_repmask2gff.pl with GFF3 output" );


#-----------------------------------------------------------+
# TE NEST CONVERSION                                        |
#-----------------------------------------------------------+

#-----------------------------+
# GFF2 OUTPUT TEST            |
#-----------------------------+
open (EXPECTED, '<data/exp/tenest_cnv_expected_gff2.txt');
my @gff2_tenest_exp = <EXPECTED>;
close (EXPECTED);
my @gff2_tenest_obs = `cnv_tenest2gff.pl -i data/HEX3045G05/tenest/HEX3045G05.LTR --gff-ver gff2`;

is_deeply ( \@gff2_tenest_obs, \@gff2_tenest_exp, "cnv_tenest2gff.pl with GFF2 output" );

#-----------------------------+
# GFF3 OUTPUT TEST            |
#-----------------------------+
open (EXPECTED, '<data/exp/tenest_cnv_expected_gff3.txt');
my @gff3_tenest_exp = <EXPECTED>;
close (EXPECTED);
my @gff3_tenest_obs = `cnv_tenest2gff.pl -i data/HEX3045G05/tenest/HEX3045G05.LTR --gff-ver gff3`;
is_deeply ( \@gff3_tenest_obs, \@gff3_tenest_exp, "cnv_tenest2gff.pl with GFF3 output" );
