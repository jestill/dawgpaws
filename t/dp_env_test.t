#!/usr/bin/perl -w
#-----------------------------------------------------------+
#                                                           |
# dp_env_test.t - Test variables in the user environment    |
#                                                           |
#-----------------------------------------------------------+
#  AUTHOR: James C. Estill                                  |
# CONTACT: jestill_at_sourceforge.net                       |
# STARTED: 12/18/2007                                       |
# UPDATED: 04/15/2009                                       |
#                                                           |
# SHORT DESCRIPTION:                                        |
#  Test to see if variables defined in the user environment |
#  are okay.                                                |
#                                                           |
#-----------------------------------------------------------+

#-----------------------------+
# INCLUDES                    |
#-----------------------------+
use Test::More tests => 19;

diag ("-----------------------------------------------------------");
diag ("Testing user environment ...");
diag ("-----------------------------------------------------------");

#-----------------------------+
# VENNMASTER                  |
#-----------------------------+
# Dropped these and added to the dp_venn_test.t script
#ok ( $ENV{VMASTER_DIR}, 
#     "VennMaster Dir Defined"  );
#ok ( $ENV{VMASTER_JAVA_BIN}, 
#     "VennMaster Java Binary Defined" );

#-----------------------------+
# FIND LTR                    |
#-----------------------------+
ok ( $ENV{FIND_LTR_PATH}, 
     "FIND LTR path defined in the environment" ) ||
    diag ("The find_ltr path can be defined with FIND_LTR_PATH");

#-----------------------------+
# APOLLO                      |
#-----------------------------+
ok ( $ENV{DP_APOLLO_BIN}, 
     "Apollo binary path defined in the environment" ) ||
    diag("The path to the Apollo binary can be defined with DP_APOLLO_BIN");

#-----------------------------+
# GENEMARK                    |
#-----------------------------+
ok ($ENV{GM_BIN_DIR},
    "GENEMARK binary directory defined in the environment") ||
    diag("The path to the genemark directory can be defined with GM_BIN_DIR");
ok ($ENV{GM_LIB_DIR},
    "GENEMARK library directory defined in the environment") ||
    diag ("The path to the genemark libary directory can be defined".
	  " with GM_LIB_DIR");

#-----------------------------+
# GENSCAN                     |
#-----------------------------+
ok ($ENV{DP_GENSCAN_BIN},
    "Path to GENSCAN binary defined in the environment") ||
    diag("The path to the Genscan library can be defined with DP_GENSCAN_BIN");
ok ($ENV{DP_GENSCAN_LIB},
    "Genscan libary directory defined in the environment") ||
    diag("The path to the genscan library directory can be defined with".
	 " DP_GENSCAN_LIB");

#-----------------------------+
# LTR FINDER                  |
#-----------------------------+
ok ($ENV{TRNA_DB},
    "LTR FINDER TRNA_DB defined in environment") ||
    diag("The path to the TRNA DB can be defined with TRNA_DB");
ok ($ENV{PROSITE_DIR}, 
    "LTR FINDER Prosite directory defined in environment") ||
    diag("The path to prosite dir can be defined with PROSITE_DIR");
ok ($ENV{LTR_FINDER},
    "LTR Finder binary path defined in environment") ||
    diag("The path to the LTR finder binary can be defined with LTR_FINDER");

#-----------------------------+
# GFF VERSION                 |
#-----------------------------+
ok ($ENV{DP_GFF},
    "Default GFF version defined in environment") ||
    diag("The default GFF version can be defined with DP_GFF");

#-----------------------------+
# SNAP BINARY PATH            |
#-----------------------------+
ok ($ENV{SNAP_PATH},
    "SNAP binary path defined in environment") ||
    diag("The path to the SNAP binary can be defined with SNAP_PATH");

#-----------------------------+
# REPEATMASKER                |
#-----------------------------+
ok ($ENV{DP_RM_BIN},
    "RepatMasker binary path defined in environment") ||
    diag("The path to the RepeatMasker binary can be defined with DP_RM_BIN");

ok ($ENV{DP_RM_DIR},
    "RepatMasker directory defined in environment") ||
    diag("The path to the RepeatMasker directory can be defined with DP_RM_DIR");

#-----------------------------+
# LTRSEQ
#-----------------------------+
ok ($ENV{LTR_SEQ_BIN},
    "Path to LTRSeq binary defined in the environment") ||
    diag("The path to the LTRSeq library can be defined with LTR_SEQ_BIN");
ok ($ENV{LTR_SEQ_DIR},
    "LTRSeq  directory defined in the environment") ||
    diag("The path to the genscan library directory can be defined with".
	 " LTR_SEQ_DIR");

#-----------------------------+
# GENEID                      |
#-----------------------------+
ok ($ENV{GENEID_BIN},
    "Path to geneid binary defined in the environment") ||
    diag("The path to the geneid binary can be defined with GENEID_BIN");

#-----------------------------+
# FINDMITE                    |
#-----------------------------+
ok ($ENV{FINDMITE_BIN},
    "Path to FindMITE binary defined in the environment") ||
    diag("The path to the geneid binary can be defined with FINDMITE_BIN");

#-----------------------------+
# FGENESH                     |
#-----------------------------+
ok ($ENV{FGENESH_PATH},
    "Path to FGENESH binary defined in the environment") ||
    diag("The path to the FGENESH library can be defined with FGENESH_PATH");

#-----------------------------+
# BLAST                       |
#-----------------------------+
ok ($ENV{DP_BLAST_BIN},
    "Path to NCBI-BLAST binary defined in the environment") ||
    diag("The path to the NCBI-BLAST library can be defined with DP_BLAST_BIN");
