#!/usr/bin/perl -w
#-------------------------------------------------------+
#                                                       |
# REPEAT MASK PARSER PIPELINE                           |
#                                                       |
#-------------------------------------------------------+
#  AUTHOR: James C. Estill                              |
# CONTACT: jestill_at_sourceforge.net                   |
# STARTED: 4/10/2006                                    |
# UPDATED: 09/28/2006                                   |
#                                                       |
# DESCRIPTION:                                          |
#  Runs the RepeatMasker program for a set of input     |
#  FASTA files against a set of repeat library files &  |
#  then converts the repeat masker *.out file into the  |
#  GFF format and then to the game XML format for       |
#  visualization by the Apollo genome anotation program.|
#                                                       |
# USAGE:                                                |
#  RepMaskParse.pl                                      |
#                                                       |
# REQUIRED SOFTWARE:                                    |
#  -Apollo (Genome Annotation Curation Tool)            |
#   http://www.fruitfly.org/annot/apollo/               |
#  -grep                                                |
#  -awk                                                 | 
#                                                       |
# NOTE:                                                 |
#  The sequence name used in the FASTA file header must |
#  be short enough to make it through to the            | 
#  RepeatMasker *.out file.                             |
#-------------------------------------------------------+


print "The RepeatMask Parser program has started\n";

#-----------------------------+
# INCLUDES                    |
#-----------------------------+
use File::Copy;

#-----------------------------+
# LOCAL VARIABLES             |
#-----------------------------+
# TempWorkDir - A directory to temporarily store the data that will later
#               be transfered to a more appropriate long term storage
#               location. This could be a place to put the original fasta file
#               and associated intermediate output.
# FileToMask - The fasta file that will be masked using RepeatMasker
# GFFOutfile - The path to write the gff outfile
# RepMaskOutfile - The path of the *.out file produced by RepeatMasker
# SeqName - The name of the sequence file in the original [name].fasta
#           file that was used as input into RepeatMasker
# RepDb - The short name of the repeat mask DB that was used

# Many of the following variables will need to be changed to set to 
# array elements, but these are hard coded for now to allow for testing and
# getting the program up and running.

# MakeGffDBCmd - Command to make a GFF file with the data grouped by the 
#                the repeat database that was used for assignment.
# MakeGffElCmd - Command to make a GFF file with the data grouped by the
#                type of element that was identified.


# OLD VARIABLES
#my $bac_data_dir = "/home/jestill/hmpr_msll/bacs/zm_testbacs/";
#my $bac_parent_dir = "/home/jestill/hmpr_msll/bacs/";
#my $local_cp_dir = "/home/jestill/hmpr_msll/bacs/zm_testbacs/";  
#my $LogFile = "/home/jestill/projects/wheat_annotation/RepMaskLog.txt";


# FROM USE ON THE CLUSTER
#my $LogFile = "/scratch/jestill/wheat/RepMaskLog.txt";
#my $bac_data_dir = "/scratch/jestill/wheat/bac2proc/";
#my $bac_parent_dir = "/scratch/jestill/wheat/";  
#my $TmpDir = "/scratch/jestill/wheat/temp/";
#my $local_cp_dir = "/scratch/jestill/wheat/copy/"; 
#my $LocalCpDir = $local_cp_dir;

# FROM THE LOCAL
my $LogFile = "/scratch/jestill/wheat/RepMaskLog.txt";
my $bac_data_dir = "/scratch/jestill/wheat/bac2proc/";
my $bac_parent_dir = "/scratch/jestill/wheat/";  
my $TmpDir = "/scratch/jestill/wheat/temp/";
my $local_cp_dir = "/scratch/jestill/wheat/copy/"; 
my $LocalCpDir = $local_cp_dir;

my ( $LibData , $RepMaskCmd, $MakeGffDbCmd, $MakeGffElCmd );
my ( @RepLibs, @FastaFiles );
my $ProcNum = 1;


#-----------------------------+
# OPEN THE LOG FILE           |
#-----------------------------+
open ( LOG, ">>$LogFile" );    # Open file for appending
$TheTime = time;
print LOG "==================================\n";
print LOG "  JOB: $TheTime\n";
print LOG "==================================\n";

#-----------------------------+
# Get the FASTA files from the|
# directory provided by the   |
# var $bac_data_dir           |
#-----------------------------+
opendir( DIR, $bac_data_dir ) || 
    die "Can't open $bac_data_dir"; 
# Read the dir ignoring . and ..
my @FastaFiles = grep !/^\.\.?$/, readdir DIR ;
closedir( DIR );

#-----------------------------+
# REPEAT LIBRARY INFORMATION  |
#-----------------------------+
# This is an array of arrays:
# The nested array should contain the following info:
# [0] - Name of the repeat library
#       This will be used to name the output files from the analysis
#       and to name the data tracks that will be used by Apollo.
# [1] - Path of the repeat library
# The number of records in the fasta file shown in brackets
# afte the description of the db
@RepLibs = (

	    #-----------------------------+
	    # SANMIGUEL REPEAT DATABASE   |
	    # 493                         |
	    #-----------------------------+
	    ["SanMiguel",     
	     "/scratch/jestill/repmask/Philipretros.fasta"],
	    
	    #-----------------------------+
	    # WESSLER REPEAT LIBRARY FROM |
	    # MAGI WEBSITE                |
	    # 571                         |
	    #-----------------------------+
#	    ["WesRep",     
#	     "/scratch/jestill/repmask/magi_Wes_Rep_db.fasta"],

	    #-----------------------------+
	    # MAGI NONREDUNDANT           |
	    # STATISTICALY DEFINED REPEATS|
	    # 13,564                     |
	    #-----------------------------+
#	    ["MAGI_NR",    
#	     "/scratch/jestill/repmask/magi_zmSDRv3_1.fasta"],

	    #-----------------------------+
	    # MAIZE TIGR 4 UNCHARACTERIZED| 
	    # REPEATS                     |
	    # 21,133                      |
	    #-----------------------------+
#	    ["TIGRUnchar", 
#	     "/scratch/jestill/repmask/tigr_uncharacterized_02202004.fasta"],

	    #-----------------------------+
	    # MAIZE TIGR 4 CHARACTERIZED  |
	    # REPEATS                     |
	    # 623                         |
	    #-----------------------------+
#	    ["TIGRChar",   
#	     "/scratch/jestill/repmask/tigr_characterized_02202004.fasta"],

	    #-----------------------------+
	    # MAIZE TIGR 4 RECON PREDICTED|
	    # 5,035                       |
	    #-----------------------------+
#	    ["TIGRRECON",  
#	     "/scratch/jestill/repmask/tigr_RECON_prediction_02202004.fasta"],

	    #-----------------------------+
	    # TIGR Brassicaceae           |
	    # 775                         |
	    #-----------------------------+
	    ["TIGRBras_2",  
	     "/scratch/jestill/repmask/TIGR_Brassicaceae_Repeats.v2.rm.fasta"],

	    #-----------------------------+
	    # TIGR Fabaceae               |
	    # 308                         |
	    #-----------------------------+
	    ["TIGRFab_2",  
	     "/scratch/jestill/repmask/TIGR_Fabaceae_Repeats.v2.rm.fasta"],

	    #-----------------------------+
	    # TIGR Gramineae              |
	    # 4,475                       |
	    #-----------------------------+
	    ["TIGRGram_2",  
	     "/scratch/jestill/repmask/TIGR_Gramineae_Repeats.v3.1.rm.fasta"],

	    #-----------------------------+
	    # TIGR Solanaceae             |
	    # 252                         |
	    #-----------------------------+
	    ["TIGRSol_2",  
	     "/scratch/jestill/repmask/TIGR_Solanaceae_Repeats.v2.rm.fasta"],

	    #-----------------------------+
	    # TIGR ARABIDOPSIS GSS REPEATS|
	    # 8,018                       |
	    #-----------------------------+
	    ["TIGR_Arab_GSS_2",  
	     "/scratch/jestill/repmask/TIGR_Arabidopsis_GSS_Repeats.v2.rm.fasta"],

	    #-----------------------------+
	    # TIGR BRASSICA GSS REPEATS   |
	    # 34,314                      |
	    #-----------------------------+
	    ["TIGR_Bras_GSS_2",  
	     "/scratch/jestill/repmask/TIGR_Brassicaceae_GSS_Repeats.v2.rm.fasta"],
	    
	    #-----------------------------+
	    # TIGR GLYCINE GSS REPEATS    |
	    # 149                         |
	    #-----------------------------+
	    ["TIGR_Gly_GSS_2",  
	     "/scratch/jestill/repmask/TIGR_Glycine_GSS_Repeats.v2.rm.fasta"],

	    #-----------------------------+
	    # TIGR LOTUS GSS REPEATS      |
	    # 479                         |
	    #-----------------------------+
	    ["TIGR_Gly_GSS_2",  
	     "/scratch/jestill/repmask/TIGR_Lotus_GSS_Repeats.v2.rm.fasta"],

	    #-----------------------------+
	    # TIGR LYCOPERSICON GSS REPEAT|
	    # 343                         |
	    #-----------------------------+
	    ["TIGR_Lyc_GSS_2",  
	     "/scratch/jestill/repmask/TIGR_Lycopersicon_GSS_Repeats.v2.rm.fasta"],

	    #-----------------------------+
	    # TIGR MEDICAGO GSS REPEATS   |
	    # 15                          |
	    #-----------------------------+
	    ["TIGR_Med_GSS_2",  
	     "/scratch/jestill/repmask/TIGR_Medicago_GSS_Repeats.v2.rm.fasta"],

	    #-----------------------------+
	    # TIGR ORYZA GSS REPEATS      |
	    # 4,627                       |
	    #-----------------------------+
	    ["TIGR_Ory_GSS_2",  
	     "/scratch/jestill/repmask/TIGR_Oryza_GSS_Repeats.v2.rm.fasta"],

	    #-----------------------------+
	    # TIGR SORGHUM GSS REPEATS    |
	    # 1,755                       |
	    #-----------------------------+
	    ["TIGR_Sor_GSS_2",  
	     "/scratch/jestill/repmask/TIGR_Sorghum_GSS_Repeats.v2.rm.fasta"],

	    #-----------------------------+
	    # TIGR TRITICUM GSS REPEATS   |
	    # 17                          |
	    #-----------------------------+
	    ["TIGR_Trit_GSS_2",  
	     "/scratch/jestill/repmask/TIGR_Triticum_GSS_Repeats.v2.rm.fasta"],

	    #-----------------------------+
	    # TIGR ZEA GSS REPEATS        |
	    # 124,724                     |
	    #-----------------------------+
	    ["TIGR_Zea_GSS_2",  
	     "/scratch/jestill/repmask/TIGR_Zea_GSS_Repeats.v2.rm.fasta"],

	    # This should remain at the bottom since we should
	    # always run the TREP database with the rice annotation.
	    #-----------------------------+
	    # TREP v 8                    |
	    # 350                         |
	    #-----------------------------+
	    ["TREP8",  
	     "/scratch/jestill/repmask/TREP_8_20050802.nr.fasta"]

	    );


my $NumLibs = @RepLibs;
my $NumFiles = @FastaFiles;
my $NumProcTotal = $NumLibs * $NumFiles;

#-----------------------------+
# RUN REPEAT MAKSER AND PARSE |
# RESULTS FOR EACH SEQ IN THE |
# FastaFiles ARRAY FOR EACH   |
# REPEAT LIBRARY IN THE       |
# RepLibs ARRAY               |
#-----------------------------+

for $SeqName (@FastaFiles)
{

    $FileToMask = $bac_data_dir.$SeqName;
    $RepMaskOutfile = $FileToMask.".out";
    $RepMaskCatFile = $FileToMask.".cat";
    $RepMaskTblFile = $FileToMask.".tbl";
    $RepMaskMaskedFile = $FileToMask.".masked";

    print LOG "Process $ProcNum of $NumProcTotal.\n";
    $ProcNum++;
    print "Process $ProcNum of $NumProcTotal.\n";


    my $bac_out_dir = $bac_parent_dir.$SeqName;
    # make the bac output dir if it does not exist    
    mkdir $bac_out_dir, 0777 unless (-e $bac_out_dir); 
    
    for $LibData (@RepLibs)
    {
	
	$RepDbName = @$LibData[0];
	$RepDbPath = @$LibData[1];
	
	#-----------------------------+
	# REMOVE ZM PREFIX FROM NAME  |
	#-----------------------------+
	# This search name will still work with
	# names that do not have unique information
	# in the first three characters
	# 11/21/2006 - TO make this work for the wheat assembled
	# BACS I will just use $SearhcName = $SeqName and see
	# how well this works
	#$SearchName = substr ($SeqName, 3);
	$SearchName = $SeqName; 
	
	$GffElOut = $bac_data_dir.$RepDbName."_".$SeqName."_EL.gff";
	$XmlElOut = $bac_data_dir.$RepDbName."_".$SeqName."_EL.game.xml"; 
 	$GffAllDbOut = $bac_data_dir."ALLDB_".$SeqName.".gff";
	$XmlAllDbOut = $bac_data_dir."ALLDB_".$SeqName."game.xml";

	#$RepMaskCmd = "RepeatMasker -fixed -lib ".$RepDbPath.
	#    " -xsmall $FileToMask";
	# USE THE FOLLOWING TO RUN ON TWO PROCESSORS
	$RepMaskCmd = "RepeatMasker -lib ".$RepDbPath.
	    " -xsmall -pa 3 $FileToMask";

# THE FOLLOWING DOES NOT APPEAR TO WORK ON ALTIX -- 09/28/2006
# USING THE -w flag will force this to use wublast
#	$RepMaskCmd = "RepeatMasker -lib ".$RepDbPath.
#	    " -xsmall -pa 4 -w $FileToMask";



	#-----------------------------+
	# COMMAND TO CONVERT OUTPUT TO|
	# GFF FORMAT AND APPEND TO    |
	# A SINGLE FILE FOR THE REPEAT|
	# LIBRARY                     |
	#-----------------------------+
	# I did not use gff out from repeatmasker because it
	# did not appear to work properly with Apollo
	$MakeGffElCmd = "grep ".$SearchName." ".$RepMaskOutfile.
	    " | awk 'BEGIN{OFS=\"\\t\"}{print \$11,  \"RepeatMasker: ".
	    $RepDbName."\", \$11, ".
	    "\$6,\$7,\$1, \".\", \".\"}' > ".$GffElOut;
	
	#-----------------------------+
	# COMMAND TO CONVERT OUTPUT TO|
	# GFF FORMAT AND APPEND TO    |
	# A SINGLE FILE FOR ALL REPEAT|
	# LIBRARIES                   |
	#-----------------------------+
	print "SEARCH: ".$SearchName."\n";
	print "OUTFILE: ".$RepMaskOutfile."\n";
	print "DB-NAME: ".$RepDbName."\n";
	print "GFF: ".$GffAllDbOut."\n";
	#exit;


	$MakeGffAllDbCmd = "grep ".$SearchName." ".$RepMaskOutfile.
	    " | awk 'BEGIN{OFS=\"\\t\"}{print \$11,  \"RepeatMasker: ".
	    $RepDbName."\", \$11, ".
	    "\$6,\$7,\$1, \".\", \".\"}' >> ".$GffAllDbOut;
	
	#-----------------------------+
	# SHOW THE USER THE COMMANDS  | 
	# THAT WILL BE USED           |
	#-----------------------------+
	print "Lib Name: ".$RepDbName."\n";
	print "Lib Path: ".$RepDbPath."\n";
	print "EL Out:   ".$GffElOut."\n";
	print "RepCmd:   ".$RepMaskCmd."\n";
	#print "GffDbCmd:\n\t".$MakeGffDbCmd."\n";
	print "GffElCmd:\n\t".$MakeGffElCmd."\n";
	print "\n\n";

	#-----------------------------+
	# PRINT INFO TO LOG FILE      | 
	#-----------------------------+
	print LOG "\tLib Name: ".$RepDbName."\n";
	print LOG "\tLib Path: ".$RepDbPath."\n";
	print LOG "\tEL Out:   ".$GffElOut."\n";
	print LOG "\tRepCmd:   ".$RepMaskCmd."\n";
	#print LOG "\tGffDbCmd:\n\t".$MakeGffDbCmd."\n";
	print LOG "\tGffElCmd:\n\t".$MakeGffElCmd."\n";
	print LOG "\n\n";

	system ( $RepMaskCmd );
	system ( $MakeGffElCmd );
	system ( $MakeGffAllDbCmd );
	
	#-----------------------------+
	# CONVERT THE FILES FROM GFF  | 
	# FORMAT TO A MORE USABLE     |
	# FORMAT FOR APOLLO           |
	#-----------------------------+

	# APOLLO FUNCTION WILL NOT WORK ON ALTIX
	# SO THIS HAS BEEN COMMENTED OUT
	#print "\n\n\n\nCONVERTING\n\n\n";
#	&ApolloConvert ( $GffElOut, "gff", $XmlElOut, "game", 
#			 $FileToMask, "none" );  
	
	#-----------------------------+
	# COPY THE RM OUTPUT FILES TO | 
	# THE RM (REPEATMASK) FOLDER  |
	# MAKE THE DIR IF NEEDED      |
	#-----------------------------+
	$bac_out_dir = $bac_parent_dir.$SeqName."/";
	$bac_rep_out_dir = "$bac_out_dir"."rm/";
	mkdir $bac_rep_out_dir, 0777 unless (-e $bac_rep_out_dir); 

	$RepMaskOutCp = $bac_rep_out_dir.$RepDbName."_".$SeqName.".out";
	$RepMaskCatCp = $bac_rep_out_dir.$RepDbName."_".$SeqName.".cat";
	$RepMaskTblCp = $bac_rep_out_dir.$RepDbName."_".$SeqName.".tbl";
	$RepMaskMaskedCp = $bac_rep_out_dir.$RepDbName."_".$SeqName.".masked";
	$RepMaskLocalCp = $LocalCpDir.$RepDbName."_".$SeqName.".masked";
	$RepMaskElCp = $bac_rep_out_dir.$RepDbName."_".$SeqName."_EL.gff";
	$RepMaskXmlElCp = $bac_rep_out_dir.$RepDbName."_".$SeqName.
	    "_EL.game.xml"; 

	
	# THE FOLLOWING ADDED 09/28/2006
	$RepMaskLog = $bac_data_dir.$RepDbName."_".$SeqName.".log";
	$RepMaskLogCp = $bac_rep_out_dir.$RepDbName."_".$SeqName.".log";
	move ( $RepMaskLog, $RepMaskLogCp);

	#-----------------------------+
	# MAKE A COPY OF THE MASKED   |
	# FASTA FILE TO A SINGLE DIR  |
	#-----------------------------+
	copy ( $RepMaskMaskedFile, $RepMaskLocalCp  ) ||
	    die "Can not copy ".$RepMaskMaskedFile." to\n".$RepMaskLocalCp;

	#-----------------------------+
	# MOVE THE RM OUTPUT FILES TO |
	# THE TARGET DIR	      |
	#-----------------------------+
	move ( $RepMaskOutfile, $RepMaskOutCp) ||
	    die "Can not move ".$RepMaskOutfile." to\n ".$RepMaskOutCp."\n";
	move ( $RepMaskCatFile, $RepMaskCatCp ) ||
	    die "Can not move ".$RepMaskCatFile."\n";
# THE TABLE FILE DID NOT GET MADE ON THE ATLITX MACHINE
#	move ( $RepMaskTblFile, $RepMaskTblCp ) ||
#	    die "Can not move ".$RepMaskTblFile."\n";
	move ( $RepMaskMaskedFile, $RepMaskMaskedCp ) ||
	    die "Can not move ".$RepMaskMaskedFile."\n";
	move ( $GffElOut , $RepMaskElCp ) ||
	    die "Can not move ".$GffElOut."\n";
# THE GAME XML FILE NOT MADE ON ALTIX
#	move ( $XmlElOut, $RepMaskXmlElCp ) ||
#	    die "Can not move ".$XmlElOut."\n";

    }
    
    #-----------------------------+ 
    # THE APOLLO CONVERT FOR THE  |
    # ENTIRE SET OF REPEAT LIBS   |
    # FOR A GIVEN SEQUENCE FILE   |
    # THAT IS BEING MASKED.       |
    #-----------------------------+
#    &ApolloConvert ( $GffAllDbOut, "gff", $XmlAllDbOut , "game", 
#    		     $FileToMask, "none" );  
    
    $RepMaskALL_GFFCp = $bac_rep_out_dir."ALLDB_".$SeqName.".gff";
    $RepMaskAll_XMLCp = $bac_rep_out_dir."ALLDB_".$SeqName."game.xml";
    move ( $GffAllDbOut, $RepMaskALL_GFFCp ) ||
	die "Can not move ".$GffAllDbOut."\n";
#    move ( $XmlAllDbOut, $RepMaskAll_XMLCp ) ||
#	die "Can not move ".$XmlAllDbOut."\n";
    
    #exit; # TEMP EXIT FOR DEBUG, WIll JUST RUN FIRST FILE TO BE MASKED


} # End of for each file in the input folder

close LOG;
exit;

#-----------------------------------------------------------+
# SUBFUNCTIONS                                              |
#-----------------------------------------------------------+

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
sub ApolloConvert
{
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

=head1 HISTORY
The program history
=cut

#-------------------------------------------------------+
# HISTORY                                               |
#-------------------------------------------------------+
#
# 4/10/2006
# - Program started
# - Parsing of RepeatMasker *.out file to an Apollo 
#   compatible *.GFF file done. Working copy for both
#   a single track for each repeat class, and a data
#   track with all of the results from a Repeat Library
#   grouped toether.
# - Two dimensional array containing the repeat library
#   database informaiton.
# 4/11/2006
# - Additional code comments and reformat
# - Code to cycle through a set of FASTA files and store
#   the output data in a separate named folder for each
#   of the FASTA flies.
# - Added the -fixed to the command to run RepeatMasker
#   to make sure that the text format is the same for all
# - Since Apollo does not allow additional GFF files to be
#   added later, I added code to create a GFF file that has
#   the outcome for all of the RepeatMask runs in one 
#   database.
# 4/12/2006
# - Added the AllRepeats to the dataset to make sure that 
#   that the repeat characterization is cumulative
# - ConvertGFF2Chado subfunction
# - Adding Smith Waterman score from RepeatMasker to the
#   GFF output file. This should be the first column in the 
#   *.out file and the 6th column in the GFF file.
# 4/16/2005
# - Added array for list of fasta files to process and testing
#   with small set of sequences.
# - The fasta file was manually edited to use just the GB
#   ID in the FASTA header. This could be done using the 
#   the seq object PERL module from bioperl.
# - The fasta files should be read from an input directory
#   and then the output should be copied to an output dir.
# 4/17/2006
# - Working out variables names
# 4/25/2006
# - Adding ability to get the input set of FASTA files
#   from a directory 
# 4/26/2006
# - Working out the copy of the relevant output files to
#   a central directory for the FASTA file (ie. all 
#   programatic output from the FASTA file goes to 
#   a dir named for the file.) 
# - Added use of File::Copy module for cp commands
# 5/24/2006
# - Added process log to the program. This will allow
#   me to launch the program at the lab and then
#   monitor the process at home.
# 5/31/2006
# - Made local copy dir a variable that can be set at the top
# 09/26/2006
# - A few changes made to the base code to make this
#   work on the altix. Using databases at
#   /scratch/jestill/repmask/
# 09/28/2006
# - Changes make to get this to work on the altix
# - Changed the format of the repeat databases list
#   This is currently a two-d array

=head1 To DO
The to do list
=cut
#-----------------------------------------------------------+
# TO DO
#-----------------------------------------------------------+
#
# - Make the results compatable for an upload to a chado
#   database.
# - This currently requires that the fasta header be the
#   exact same name as the, this should be fixed such
#   that the program can juse read the name. There also
#   appears to be a name trimming problen with the contigs.
# - MUST USE SHORT NAMES
