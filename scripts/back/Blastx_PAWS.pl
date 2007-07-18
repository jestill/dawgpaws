#!/usr/bin/perl -w
#-----------------------------------------------------------+
#                                                           |
# BLAST PAWS                                                |
#                                                           |
#-----------------------------------------------------------+
#  AUTHOR: James C. Estill                                  |
# CONTACT: jestill_at_sourceforge.net                       |
# STARTED: 05/10/2006                                       |
# UPDATED: 02/13/2007                                       |
#                                                           |
# DESCRIPTION:                                              |
#  Run NCBI blast for the asgr sequences against a set of   |
#  repeast databases of interest.                           |
#                                                           |
# USAGE:                                                    |
#                                                           |
#-----------------------------------------------------------+

# Can run this with a test flag to just run through
# the set to make sure that all of the databases and input
# files exist.

print "The program has started\n";

#-----------------------------+
# INCLUDES                    |
#-----------------------------+
use Getopt::Std;               # Allows to get options from the command line
use Text::Wrap;                # Allows word wrapping and hanging indents
                               # for more readable output for long strings.

#-----------------------------+
# LOCAL VARIABLES             |
#-----------------------------+
my $QDir;                     # Base dir for qry sequences  
my $DbDir;                    # Base dir for BLAST databases
my $BlSuf;                    # Blast suffix
my $BlProg;                   # Blast program to use 
my ( $BlastCmd, $QryPath, $DbPath, $OutPath ); 

$QDir = "/home/jestill/projects/wheat_annotation/".
    "wed/Wed/Bac2Pro/"; 
$DbDir = "/home/jestill/blast/paws/";
# E is 10^-5
# Run on both processors
# Give tab delim output
# -u Uses lowercase filtering
#$BlSuf = "-e 0.00001 -a 2 -m 8 -U";
# Full text blast output
#$BlSuf = "-e 0.00001 -a 2 -U";
# Full text an not MASKED 
$BlSuf = "-e 0.00001 -a 2 -U";

$BlProg = "tblastx";
my $ProcNum = 0;

#-----------------------------+
# COMMAND LINE VARIABLES      |
#-----------------------------+
# Options hash to hold the options from the command line
my %Options;                  
getopts('t', \%Options);      # Get the options from the command line
my $test;
$test = $Options{t};

#-----------------------------+
# BASE NAMES OF FILES TO USE  | 
# AS QUERY SETS FOR BLAST     |
#-----------------------------+
# TODO Have this search the dir for fasta files
# At present use seq name without fasta extension
@Qry = ( "DQ537335",
	 "DQ537336",
	 "DQ537337"
	 );
  
#-----------------------------+
# BLAST QUERY DATABASES       |
#-----------------------------+
@Db = ( 
        #---------------------+
	# TIGR GENE INDICES   |
	#---------------------+
#	"TaGI_10",
#	"AtGI_13",
#	"ZmGI_17",
#	"SbGI_8",
#	"OsGI_17",
#	"HvGI_9",
	#---------------------+
	# CDS MODEL SEQS      |
#       #---------------------+
	"ATH1_cds",            # Arabidopsis
	"Os_5_cds"             # Rice
	#---------------------+
	# REPEAT DATABASES    |
	#---------------------+
#	"mips_REdat_4_3",
#	"RB_pln",
#	"SanMiguel_200610",
#	"TIGR_Gram_3_3",
#	"TREP9_nr",
#	"TREP9_total",
#	"Wessler"
	);

my $LenQry = @Qry;
my $LenDb =  @Db;
my $NumProc = $LenQry * $LenDb;

for $IndQry (@Qry)
{

    for $IndDb (@Db)
    {

	$ProcNum++; # Increment the process number
	print "BLAST Process ".$ProcNum." of ".$NumProc."\n";

	$QryPath = $QDir.$IndQry."/".$IndQry.".fasta";
	$DbPath = $DbDir.$IndDb;
	$TestFile = $DbPath.".nhr";

#	$OutPath = $QDir.$IndQry."/".$IndQry."_".$IndDb.".blo";
# FULL TEXT BLAST OUTPUT
#	$OutPath = $QDir.$IndQry."/".$IndQry."_".$IndDb.".full.blo";
# FULL TEXT NOT MASKED
	$OutPath = $QDir.$IndQry."/".$IndQry."_".$IndDb.".full.blx";

	#-----------------------------+
	# DOES THE BLAST DB EXIST     |
	#-----------------------------+
	# CHANGE THIS TO AN UNLESS STATEMENT
	if (-e $TestFile ) 
	{
	    #print "DB: $IndDb exists\n";
	}
	else 
	{die "Can not find database:\n$IndDb\n"; }
	
	#-----------------------------+
	# DOES THE BLAST QRY EXIST    |
	#-----------------------------+
	# CHANGE THIS TO AN UNLESS STATEMENT
	if (-e $QryPath)
	{
	    #print "QRY: $QryPath exists\n";
	}
	else
	{die "Can not find qry file:\n$QryPath\n";}

	#------------------------------+
	# PRINT THE BLAST COMMAND      |
	#------------------------------+
	$BlastCmd = "blastall -p ".$BlProg." -i $QryPath -d $DbPath -o $OutPath $BlSuf";
	print wrap("\t", "\t", $BlastCmd );
	print "\n";
	#------------------------------+
	# RUN THE BLAST COMMAND        |
	#------------------------------+
	if (! $test)  # If this is not a test then run the BlastCmd
	{
	    system ($BlastCmd);
	}
	


    } # End of for each database loop

} # End of for each query seq loop


exit;


#-----------------------------------------------------------+
# SUBFUNCTIONS                                              |
#-----------------------------------------------------------+


#-----------------------------------------------------------+
# HISTORY                                                   |
#-----------------------------------------------------------+
#
# 05/10/2006
# - Program was started to create an easy way to blast all
#   of the query datasets of interest by all of the repeat
#   databases of interest. The BLAST reports should then
#   be parsed and uploaded to the dbASGR database.
# 02/13/2007
# - Existing program was modified for the PAWS program
