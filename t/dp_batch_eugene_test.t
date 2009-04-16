#!/usr/bin/perl -w
#-----------------------------------------------------------+
#                                                           |
# dp_batch_eugene_test.pl - Test batch_eugene.pl            |
#                                                           |
#-----------------------------------------------------------+
#  AUTHOR: James C. Estill                                  |
# CONTACT: jestill_at_sourceforge.net                       |
# STARTED: 04/16/2009                                       |
# UPDATED: 04/16/2009                                       |
#                                                           |
# SHORT DESCRIPTION:                                        |
#  Test the batch_eugene.pl program.                        |
#                                                           |
#-----------------------------------------------------------+

use Test::More tests => 1;

#-----------------------------+
# TESTING BINARY PATH         |
#-----------------------------+
my $eugene_bin = "eugene";
my $get_sites =  "egn_getsites4eugene.pl";
my $param_file = "hex_eugene.par";

# Basic command to see that eugene is installed
my $basic_eug_cmd = "eugene";
ok ( system($basic_eug_cmd), "Eugene is Installed" ) ||
    diag ("You must have eugene installed, and in your PATH.");



# Use the parameters file
#my $eugene_cmd = $eugene_bin." -p g -A $param_file ".
#    " -O $eugene_dir".
#    " $eug_fasta";

exit;

#-----------------------------+
# TEST RUNNING THE PROGRAM    |
#-----------------------------+
diag ("Testing running the batch_eugene.pl program, this is slow ...");
my $out_dir = $ENV{HOME}."/dp_temp_test/";
my $eugene_cmd = "batch_eugene.pl -i data/fasta/ -o $out_dir".
    " -p ";

# This will exit zero when things work
ok ( system($findgaps_cmd)==0 , "batch_eugene.pl");

#-----------------------------+
# TEST PROGRAM OUTPUT         |
#-----------------------------+
diag ("Testing gap results ..");
my $exp_file = "data/exp/HEX2903P03_gaps.gff";

open (EXPECTED, "<".$exp_file);
my @gaps_exp = <EXPECTED>;
close (EXPECTED);

# TEST THE THE GFF OUTPUT FILE EXISTS
my $obs_file = $out_dir."HEX2903P03/gff/HEX2903P03_gaps.gff";
ok ( (-e $obs_file) , "Gaps GFF test files appears to exist") ||
    diag("I expected to see the file\n $obs_file");

open (OBSERVED, "<".$obs_file);
my @gaps_obs = <OBSERVED>;
close (OBSERVED);

is_deeply ( \@gaps_obs, \@gaps_exp,
	    "Gaps GFF file contains correct data");
