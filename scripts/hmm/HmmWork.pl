#!/usr/bin/perl -w
#-----------------------------------------------------------+
#                                                           |
# HMMER REPEATS                                             |
#                                                           |
#-----------------------------------------------------------+
#  AUTHOR: James C. Estill                                  |
# CONTACT: jestill_at_sourceforge.net                       |
# STARTED: 08/07/2006                                       |
# UPDATED: 08/07/2006                                       |
#                                                           |
# DESCRIPTION:                                              |
# Run HMMMER program to identify repeats. This will         |
# initially just use the models from the Bureau lab rice    |
# HMM models.                                               |
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
# DIR FOR ALL ASGR BLAST ANALYSIS QUERY SEQUENCES
$QDir = "/home/jestill/projects/asgr/hmm_mite/seqs/";  # Base dir for query sequences
$ModDir = "/home/jestill/projects/asgr/hmm_mite/mite_models/";  # Base dir for HMM models
$ProcNum = 0;                                          # Process number starts at zero

my ( $HmmCmd, $QryPath, $DbPath, $OutPath );  # Declare scope for varaiables used later

#-----------------------------+
# COMMAND LINE VARIABLES      |
#-----------------------------+
my %Options;                  # Options hash to hold the options from the command line
getopts('t', \%Options);      # Get the options from the command line
my $test;
$test = $Options{t};

#-----------------------------+
# BASE NAMES OF FILES TO USE  | 
# AS QUERY SEQS FOR HMMER     |
#-----------------------------+
@Qry = ( "PPAN058CON001"
	 );
		  
#-----------------------------+
# HMMER MODELS                |
#-----------------------------+

opendir( MODDIR, $ModDir );
# Read the directory ignoring the . and ..
my @Mod = grep !/^\.\.?$/, readdir MODDIR ;
closedir( MODDIR );    


# Determine the total of BLAST queries that will be run
my $LenQry = @Qry;
my $LenMod =  @Mod;
my $NumProc = $LenQry * $LenMod;
#exit;

for $IndQry (@Qry)
{

    for $IndMod (@Mod)
    {

	$ProcNum++; # Increment the process number
	print "HMM Process ".$ProcNum." of ".$NumProc."\n";

	$QryPath = $QDir.$IndQry.".fasta";
	$ModPath = $ModDir.$IndMod;
	$ModFile = $ModPath;
	$OutPath = $QDir.$IndQry."/".$IndQry."_".$IndMod.".hmmout";

	#-----------------------------+
	# DOES THE HMM MODEL EXIST    |
	#-----------------------------+
	if (-e $ModFile ) 
	{
	    #print "DB: $IndDb exists\n";
	}
	else 
	{die "Can not find model:\n$IndMod\n"; }
	
	#-----------------------------+
	# DOES THE HMM QRY SEQ EXIST  |
	#-----------------------------+
	if (-e $QryPath)
	{
	    #print "QRY: $QryPath exists\n";
	}
	else
	{die "Can not find qry file:\n$QryPath\n";}

	#------------------------------+
	# PRINT THE BLAST COMMAND      |
	#------------------------------+
	$HmmCmd = "hmmsearch " . 
	    "--domT 2 $ModFile $QryPath >$OutPath";
	print wrap("\t", "\t", $HmmCmd );
	print "\n";

	#------------------------------+
	# RUN THE BLAST COMMAND        |
	#------------------------------+
	if (! $test)  # If this is not a test then run the BlastCmd
	{
	    system ($HmmCmd);
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
