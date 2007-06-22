#!/usr/bin/perl -w
#-----------------------------------------------------------+
#
# GFF2Game.pl
#
#-----------------------------------------------------------+
#
#  AUTHOR: James C. Estill
# STARTED: 02/09/2007
# UPDATED: 02/09/2007
# DESCRIPTION:
#  Converts gff data tracks from gff format to the game
#  xml format for use in the Apollo Genome Annotation 
#  Curation program.
# 
# DEPENDENCIES:
# *Apollo
#    Requires that apollo be installed on the local machine
#    since apollo is being used as the engine to do the
#    conversion between formats.
#-----------------------------------------------------------+ 
# TO DO:
# The following three should be variables that can be passed
# at the command line if desired.

my $InFastaFile = "/home/jestill/projects/wheat_annotation/wed/Wed/".
    "DQ537336/DQ537336.fasta";

my $InputRoot = "/home/jestill/projects/wheat_annotation/wed/Wed/".
    "DQ537336/JamesEstill17103038/";

my $OutputRoot = "/home/jestill/projects/wheat_annotation/wed/Wed/".
    "DQ537336/game/";

# The following are typical for the
# TriAnnot Pipeline
#my $FileRoot = "1blastnRpb";       # RepBase

#
# TO DO:
# * Load this array using the command to list the gff files in the
#   directory
# * Also use the ls command to get the fasta file to use with this
# 

my @Files2Convert = (
		     "1blastnTIGRRep",    # TIGR Repeat Database
		     "1blastnTrp",        # Dont Know
		     "1blastnTtot",       # TREP total
		     "1blastxTpr",        # TREP protein
		     "1trf",              # Trandem Repeat Finder
		     "2fGh",              # Fgenesh - Gene Model Prediction
		     "2gID",              # GeneID - Gene Model Prediction
		     "2gmHv",             # GeneMark Hordeum vulgare
		     "2gmOs",             # GeneMark Oryza sative
		     "2gmTa",             # GeneMark Triticum aestevum
		     "2gmZm",             # GeneMark Zea mays 
		     "3blastnOsU",        # Blastn Oryza sativa Unigenes
		     "3blastnTaU",        # Blastn Triticum aestivm Unigenes
		     "3blastnZmU",        # Blastn Zea mays Unigenes
		     "3blastxSprot",      # Blastx
		     "3blastxTrembl",     # Blastx
 		     "3tblastxTAta"       # tBLASTx
		     );

# Can do a while $FileRoot 
# across the array to conver here

# !!! APOLLO IS HAVING TROUBLE PARSING THE SCORE FIELD
#     FROM PROPERLY FORMATTED GFF FILES THIS NEEDS TO BE FIXED

foreach my $FileRoot (@Files2Convert)
{
    
    my $InputTest = $InputRoot.$FileRoot.".gff";
    my $OutputTest = $OutputRoot.$FileRoot.".game.xml";
    
    # Convert each gff file game xml using the ApolloConvert command
    &ApolloConvert ($InputTest, "gff", $OutputTest, 
		    "game", $InFastaFile, "NULL");
}

# Add bioperl dependent code here to concatenate the game xml files
# to a single XML file name for the root name in the fasta file or
# a root name given by the user.

exit;

#-----------------------------------------------------------+
# SUBFUNCTIONS                                              |
#-----------------------------------------------------------+

sub ApolloConvert
{
#-----------------------------+
# CONVERT AMONG FILE FORMATS  |
# USING THE APOLLO PROGRAM    |
#-----------------------------+
# Converts among the various data formats that can be used 
# from the command line in tbe Apollo program. For example
# can convert GFF format files into the game XML format.
# NOTES:
#  - Currently assumes that the input file is in the correct
#    coordinate system.
#  - GFF files will require a sequence file
#  - ChadoDB format will require a db password


    # ApPath - the path of dir with the Apollo binary
    #          Specifying the path will allow for cases
    #          where the program is not in the PATHS
    # ApCmd  - the apollo commands to run

    my $InFile = $_[0];        # Input file path
    my $InForm = $_[1];        # Output file format:
                               # game|gff|gb|chadoxml|backup
    my $OutFile = $_[2];       # Output file path
    my $OutForm = $_[3];       # Ouput file foramt
                               # chadoDB|game|chadoxml|genbank|gff|backup
    my $SeqFile = $_[4];       # The path of the sequence file
                               # This is only required for GFF foramt files
                               # When not required this can be passed as na
    my $DbPass = $_[5];        # Database password for logging on to the 
                               # chado database for reading or writing.
    my ( $ApPath, $ApCmd );

    $ApPath = "/home/jestill/Apps/Apollo/";

    # Set the base command line. More may need to be added for different
    # formats. For example, GFF in requires a sequence file and the CHADO
    # format will require a database password.
    $ApCmd = $ApPath."Apollo -i ".$InForm." -f ".$InFile.
	" -o ".$OutForm." -w ".$OutFile;

    # Make sure that that input output formats are in lowercase
    # may need to add something here to avoid converting chadoDB
    $InForm = lc($InForm);
    $OutForm = lc($OutForm);
    
    # Determine the proper command to use based on the input format
    # since GFF file also require a sequence file
    if ($InForm =~ "gff" )
    {
	$ApCmd = $ApCmd." -s ".$SeqFile;
    }
    
    if ($InForm =~ "chadodb")
    {
	$ApCmd = $ApCmd." -D ".$DbPass;
    }

    # Do the apollo command
    system ( $ApCmd );

}


# NOTES CONCERNING GFF IN APOLLO
# GFF in apollo is completed screwy, and does not
# use the standard definition of the GFF format when
# importing GFF data
# Typics

# This is a complete pain in the ass so we need
# to convert from standard GFF-1 to apollo formatted
# GFF before we can use the GFF import funcion in apoll0
# APOLLO "GFF" IS
# The following is assuming that the AutoMask.pl program
# puts everything in the correct place
# [0] - The name you want to group your file by
# [1] - SOURCE
# [2] - FEATURE
# [3] - START
# [4] - END
# [5] - SCORE
# [6] - STRAND
# [7] - FRAME
# HOWEVER FOLLOWING THIS FORMAT DOES NOT APPEAR TO PARSE
# THE SCORE, STILL GET ERROR
