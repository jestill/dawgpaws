#!/usr/bin/perl -w
#-----------------------------------------------------------+
#                                                           |
# dp_test_all.pl - Run all DAWGPAWS tests.                  |
#                                                           |
#-----------------------------------------------------------+
#  AUTHOR: James C. Estill                                  |
# CONTACT: jestill_at_sourceforge.net                       |
# STARTED: 04/15/2009                                       |
# UPDATED: 04/15/2009                                       |
#                                                           |
# SHORT DESCRIPTION:                                        |
#  Run all of the test scripts for the DAWGPAWS programs.   |
#-----------------------------------------------------------+

use TAP::Harness;

print STDERR "Test\n";

my %args = (verbosity =>1, color=>1);

my $harness = TAP::Harness->new( \%args );

$harness->runtests( 
    ['dp_module_test.t', 'TESTING DAWGPAWS REQUIRED MODULES'],
    ['dp_env_test.t', 'TESTING DAWGPAWS ENVIRONMENT'],
    ['dp_cnv_test.t', 'TESTING DAWGPAWS CONVERTER SCRIPTS'],
    ['dp_cnv_game_test.t', 'TESTING GFF2GAME CONVERTER'],
    ['dp_venn_test.t', 'TESTING VENN PROGRAM'],
    ['dp_gff_seg_test.t', 'TESTING GFF SEGMENTATION'],
    );
