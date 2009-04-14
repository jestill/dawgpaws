#!/usr/bin/perl -w
#-----------------------------------------------------------+
#                                                           |
# dp_env_test.t - Test variables in the user environment    |
#                                                           |
#-----------------------------------------------------------+
#  AUTHOR: James C. Estill                                  |
# CONTACT: jestill_at_sourceforge.net                       |
# STARTED: 12/18/2007                                       |
# UPDATED: 12/18/2007                                       |
#                                                           |
# SHORT DESCRIPTION:                                        |
#  Test to see if variables defined in the user environment |
#  are okay.                                                |
#                                                           |
#-----------------------------------------------------------+

#-----------------------------+
# INCLUDES                    |
#-----------------------------+
use Test::More tests => 11;

print "Testing user environment ...\n";

#-----------------------------+
# VENNMASTER                  |
#-----------------------------+
ok ( $ENV{VMASTER_DIR}, 
     "VennMaster Dir Defined"  );
ok ( $ENV{VMASTER_JAVA_BIN}, 
     "VennMaster Java Binary Defined" );

#-----------------------------+
# FIND LTR                    |
#-----------------------------+
ok ( $ENV{FIND_LTR_PATH}, 
     "FIND LTR PATH Defined" );

#-----------------------------+
# APOLLO                      |
#-----------------------------+
ok ( $ENV{DP_APOLLO_BIN}, 
     "Apollo Binary Path Defined" );

#-----------------------------+
# GENEMARK                    |
#-----------------------------+
ok ($ENV{GM_BIN_DIR},
    "GENEMARK Binary Defined");
ok ($ENV{GM_LIB_DIR},
    "GENEMARK Library defined");

#-----------------------------+
# GENSCAN                     |
#-----------------------------+
ok ($ENV{DP_GENSCAN_BIN},
    "Genscan Binary Defined");
ok ($ENV{DP_GENSCAN_LIB},
    "Genscan libary directory defined");

#-----------------------------+
# LTR FINDER                  |
#-----------------------------+
ok ($ENV{TRNA_DB},
    "LTR FINDER TRNA_DB Defined");
ok ($ENV{PROSITE_DIR}, 
    "LTR FINDER Prosite Dir Defined");
ok ($ENV{LTR_FINDER},
    "LTR Finder binary path Defind");
