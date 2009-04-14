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
#  Test to see if modules required by DAWG-PAWS are         |
#  present.                                                 |
#-----------------------------------------------------------+

use Test::More tests => 1;

print "Testing ab initio gene parsers ...\n";

#-----------------------------+
# FGENESH                     |
#-----------------------------+

