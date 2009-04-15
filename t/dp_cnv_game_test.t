#!/usr/bin/perl -w
#-----------------------------------------------------------+
#                                                           |
# dp_gene_parse_test.pl - Test gene prediction converters   |
#                                                           |
#-----------------------------------------------------------+
#  AUTHOR: James C. Estill                                  |
# CONTACT: jestill_at_sourceforge.net                       |
# STARTED: 04/14/2009                                       |
# UPDATED: 04/14/2009                                       |
#                                                           |
# SHORT DESCRIPTION:                                        |
#  Test of the converts in the dawpaws program.             |
#-----------------------------------------------------------+

use Test::More tests => 3;

print "Testing Apollo dependent game.xml converters ...\n";

#-----------------------------+
# TESTING APOLLO BIN DEFINED  |
#-----------------------------+
my $ap_bin;
ok ( $ap_bin = $ENV{DP_APOLLO_BIN} || "apollo", 
     "Apollo Binary Path Defined as $ap_bin" );

#-----------------------------+
# APOLLO TEST COMMAND         |
#-----------------------------+
my $ap_test_cmd = "$ap_bin --help";
ok ( my $apollo_result = `$ap_test_cmd`,
    "Apollo test command");

#-----------------------------+
# cnv_gff2game.pl             |
#-----------------------------+
print STDERR "Testing apollo mediated conversion ..\n";
open (EXPECTED, '<data/exp/genscan_gamexml_expected.txt');
my @game_xml_exp;
my $in_analysis = 0;
while (<EXPECTED>) {

    if ($in_analysis) {
	push (@game_xml_exp, $_);
    }
    
    if (m/\<computational_analysis\>/) {
	$in_analysis = 1;
    }
    elsif (m/\<\/computational_analysis\>/) {
	$in_analysis = 0;
    }
    
}
close (EXPECTED);

my $num_lines = @game_xml_expected;
print STDERR "Lines $num_lines\n";

# OBSERVED
my $tmp_result = "~/tmp_gamexml_observed.txt";
my $cnv_cmd = "cnv_gff2game.pl".
    " -i data/HEX3045G05/HEX3045G05_TREP9.masked.fasta".
    " -g data/HEX3045G05/gff/HEX3045G05.genscan.gff".
    " -o $tmp_result";
system ($cnv_cmd);

my $obs_file = $ENV{HOME}."/tmp_gamexml_observed.txt";
open (OBSERVED, '<'.$obs_file) ||
    die "Can not open the result at $tmp_result\n";
my @game_xml_obs;
while (<OBSERVED>) {

    if ($in_analysis) {
	push (@game_xml_obs, $_);
    }
    
    if (m/\<computational_analysis\>/) {
	$in_analysis = 1;
    }
    elsif (m/\<\/computational_analysis\>/) {
	$in_analysis = 0;
    }
    
}
close (OBSERVED);

is_deeply ( \@game_xml_obs, \@game_xml_exp,
	    "cnv_gff2game.pl");
