#!/usr/bin/perl -w
#-----------------------------------------------------------+
#                                                           |
# ApPushFiles                                               |
#                                                           |
#-----------------------------------------------------------+
#                                                           |
#  AUTHOR: James C. Estill                                  |
# STARTED: 02/11/2007                                       |
# UPDATED: 02/11/2007                                       |
# DESCRIPTION:                                              |
#  Pushes the files needed for Apollo to the computers where|
#  the annotation team is working (Currently C116).         |
#                                                           |
# DEPENDENCIES:                                             |
#  * Net::SFTP may need to be installed from CPAN           | 
#-----------------------------------------------------------+ 

#-----------------------------+
# INCLUDES                    |
#-----------------------------+
use Getopt::Std;               # Options flags at command line
use Net::SFTP;                 # Secure FTP connects
use warnings;
#use strict;

#-----------------------------+
# SET VARIABLE SCOPE          |
#-----------------------------+

#-----------------------------+
# GET VARS FROM COMMAND LINE  |
#-----------------------------+
my %Options;
getopts('u:q', \%Options);
my $Usage = "ApPushFiles.lp -u UserName -q\n";
my $UserName = $Options{u} ||
    die "You must provide a user name for the SFTO connection.\n$Usage\n";
my $Quiet = $Options{q};

#-----------------------------+
# GET PASSWORD FROM USER      |
#-----------------------------+
print "\nPassword for $UserName\n";
system('stty', '-echo') == 0 or die "can't turn off echo: $?";
my $UserPassword = <STDIN>;
system('stty', 'echo') == 0 or die "can't turn on echo: $?";
chomp $UserPassword;

#-----------------------------+
# ARRAY VARIABLES             |
#-----------------------------+
# The list of target IP addresses for pushing
# files to.
my @Targets = ( 
#		"128.192.158.11",
#		"128.192.158.12",
#		"128.192.158.13",
#		"128.192.158.14",
#		"128.192.158.15",
#		"128.192.158.16",
#		"128.192.158.17"
		"128.192.158.12" # SINGLE MACHINE CHECK
		);

#-----------------------------+
# FILE TRANSFER DATA          |
#-----------------------------+
# This is a two-dimensional array with 
# each row containing the information:
# (SrcPath, DestinationPath)
my @Tr = (
	  ["/home/jestill/.apollo/wheat.tiers" ,
	   ".apollo/wheat_test_1.tiers"]
	  );
print "Test transfer file".$Tr[0][0]."\n";
print "Test transfer file".$Tr[0][1]."\n";

$TrLen = @Tr;
#print "NUMBER OF FILES TO TRANSFER".$TrLen."\n";
#$MaxTr = $TrLen - 1; 
#print "MAX array position is ".$MaxTr."\n";
#exit; 

#my $LocFile = "/home/jestill/.apollo/wheat.tiers";
#my $RemFile = ".apollo/wheat.tiers";

# Get the number of files to transfer 
# in the array @Tr

#-----------------------------+
# SOURCE FILE VALIDATION      |
#-----------------------------+
# Validate the existence of the source files
# Throw die error for file not -e

# For each individual target machine
# in the Targets array
foreach my $IndTarget (@Targets)
{
    my $sftp = Net::SFTP->new( $IndTarget, 
			       user=> "ateam",
			       password=> "annytait2",
			       debug=>"true") 
	|| die "Net::SFTP not working\n";
    

    # ATTEMPTING TO REFERENCE THE ARRAY ABOVE
    $sftp->put($Tr[0][0], $Tr[0][1]) ||
	die "Transfer of $Tr[0][0] not working\n";
    
    print "Completed push to: ".$IndTarget."\n";
    
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


