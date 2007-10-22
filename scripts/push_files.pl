#!/usr/bin/perl -w
#-----------------------------------------------------------+
#                                                           |
# push_files.pl                                             |
#                                                           |
#-----------------------------------------------------------+
#                                                           |
#  AUTHOR: James C. Estill                                  |
# CONTACT: jestill_atthehere_uga.edu                        |
# STARTED: 02/11/2007                                       |
# UPDATED: 02/28/2007                                       |
#                                                           |
# DESCRIPTION:                                              |
#  Pushes the files needed for Apollo to the computers      |
#  where the annotation team is working (Currently C116).   |
#  This will be used to automatically update the necessary  |
#  files in the home dir for these machines.                |
#                                                           |
# DEPENDENCIES:                                             |
#  * Net::SFTP may need to be installed from CPAN           |
#                                                           |  
#-----------------------------------------------------------+ 

#-----------------------------+
# INCLUDES                    |
#-----------------------------+
use Getopt::Std;               # Options flags at command line
use Net::SFTP;                 # Secure FTP connects
use warnings;
use utf8;                      # Use utf8 encoding to try to 
                               # fix log-on problem  
#use strict;

#-----------------------------+
# SET VARIABLE SCOPE          |
#-----------------------------+
my $DebugVal;                 # True to use debug on SFTP, false otherwise
my $sftp;                     # SFTP connection object

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
# Since the password is coming from the command line
# this may not be in utf8 (ie. may be ansi)
print "\nPassword for $UserName\n";
system('stty', '-echo') == 0 or die "can't turn off echo: $?";
my $UserPassword = <STDIN>;
system('stty', 'echo') == 0 or die "can't turn on echo: $?";
chomp $UserPassword;           # Remove newline character
# Convert to utf8 to avoid problem with DES.pm module
utf8::encode ($UserPassword); 

#-----------------------------+
# ARRAY VARIABLES             |
#-----------------------------+
# The list of target IP addresses for pushing
# files to.
# This should be changed to an external text file.
# Replace the number symbols below with IP addresses
my @Targets = ( 
		"###.###.###.##",
		"###.###.###.##",
		);

#-----------------------------+
# FILE TRANSFER DATA          |
#-----------------------------+
# This is a two-dimensional array with 
# each row containing the information:
# [SrcPath, DestinationPath1,DestinationPath2]
# The use of multiple destination paths allows me
# to use the same script to pus files to multiple IP
# sets.
my @Tr = (
	  #-----------------------------+
	  # WHEAT TIERS FILE            |
	  #-----------------------------+
#	  ["/home/jestill/.apollo/wheat.tiers",   # Src file
#	   ".apollo/wheat.tiers",                 # C116 Destination
#	   "ApolloFiles/wheat.tiers"],            # Local Ateam
	  #-----------------------------+
	  # FLY.STYLE FILE              |
	  #-----------------------------+
#	  ["/home/jestill/.apollo/fly.style",     # Src file
#	   ".apollo/fly.style",                   # C116 Destination 
#	   "ApolloFiles/fly.style"]               # Local Ateam
	  #-----------------------------+
          # BAC TO ANNOTATE             |
	  #-----------------------------+
#	  ["/home/jestill/projects/wheat_annotation/200702/DirtyDozen/HEX1057D07.game.xml",     # Src file
#	   "HEX1057D07.game.xml",                   # C116 Destination 
#	   "DirtyDozen/HEX1057D07.style"]               # Local Ateam
	  # ANOTHER BAC TO ANNOTATE
#	  ["/home/jestill/projects/wheat_annotation/200702/DirtyDozen/HEX1057D14.game.xml",     # Src file
#	   "HEX1057D14.game.xml",                          # C116 Destination 
#	   "DirtyDozen/HEX1057D14.game.xml"]               # Local Ateam

#	  ["/home/jestill/projects/wheat_annotation/200702/DirtyDozen/HEX1057C05.game.xml",     # Src file
#	   "HEX1057C05.game.xml",                          # C116 Destination 
#	   "DirtyDozen/HEX1057C05.game.xml"]               # Local Ateam
#
#	  ["/home/jestill/projects/wheat_annotation/200702/DirtyDozen/HEX1057B03.game.xml",     # Src file
#	   "HEX1057B03.game.xml",                          # C116 Destination 
#	   "DirtyDozen/HEX1057B03.game.xml"]               # Local Ateam

#	  ["/home/jestill/projects/wheat_annotation/200702/DirtyDozen/HEX0961H02.game.xml",     # Src file
#	   "HEX0961H02.game.xml",                          # C116 Destination 
#	   "DirtyDozen/HEX0961H02.game.xml"]               # Local Ateam

#	  ["/home/jestill/projects/wheat_annotation/200702/DirtyDozen/HEX0961G12.game.xml",     # Src file
#	   "HEX0961G12.game.xml",                          # C116 Destination 
#	   "DirtyDozen/HEX0961G12.game.xml"]               # Local Ateam

#	  ["/home/jestill/projects/wheat_annotation/200702/DirtyDozen/HEX0961F17.game.xml",     # Src file
#	   "HEX0961F17.game.xml",                          # C116 Destination 
#	   "DirtyDozen/HEX0961F17game.xml"]               # Local Ateam

#	  ["/home/jestill/projects/wheat_annotation/200702/DirtyDozen/HEX0961F07.game.xml",     # Src file
#	   "HEX0961F07.game.xml",                          # C116 Destination 
#	   "DirtyDozen/HEX0961F07game.xml"]               # Local Ateam

#	  ["/home/jestill/projects/wheat_annotation/200702/DirtyDozen/HEX0961E20.game.xml",     # Src file
#	   "HEX0961E20.game.xml",                          # C116 Destination 
#	   "DirtyDozen/HEX0961E20.game.xml"]               # Local Ateam

#	  ["/home/jestill/projects/wheat_annotation/200702/DirtyDozen/HEX0593O07.game.xml",     # Src file
#	   "HEX0593O07.game.xml",                          # C116 Destination 
#	   "DirtyDozen/HEX0593O07.game.xml"]               # Local Ateam

#	  ["/home/jestill/projects/wheat_annotation/200702/DirtyDozen/HEX0075G21.game.xml",     # Src file
#	   "HEX0075G21.game.xml",                          # C116 Destination 
#	   "DirtyDozen/HEX0075G21.game.xml"]               # Local Ateam

 	  ["/home/jestill/projects/wheat_annotation/200704/HEX0014K09/HEX0014K09.game.xml",     # Src file
 	   "HEX0014K09.game.xml",                          # C116 Destination 
 	   "HEX0014K09.game.xml"]               # Local Ateam

	  );

$TrLen = @Tr;                  # Number of files to transfer
$MaxTr = $TrLen - 1;           # Where to stop in array value

#-----------------------------+
# SOURCE FILE VALIDATION      |
#-----------------------------+
# Validate the existence of the source files
# Throw exit for file not -e
print "Running local file check ...";
for ($i=0; $i<=$MaxTr ; $i++)
{
    unless (-e $Tr[$i][0])
    {
	"\nThe source file could not be located at".$Tr[$i][0]."\n";
	exit;
    }
}
print " ok\n";

print "Number of files to transfer: ".$TrLen."\n\n";

#print "Debug Val".$DebugVal."\n";
#exit;


#-----------------------------+
# C116 Computers              |
#-----------------------------+
# For each individual target machine
# in the Targets array
foreach my $IndTarget (@Targets)
{
    
    #-----------------------------+
    # CONNECT, USE DEBUG UNLESS   |
    # THE QUIET FLAG WAS USED     |
    #-----------------------------+
    # For some reason passing the password at the command line
    # throws an error (input must be 8 bytes long). Perhaps
    # this is an ANSI thing
    if ($Quiet)
    {
	print "Connecting to $IndTarget\n"; 
	$sftp = Net::SFTP->new( $IndTarget, 
				user=> $UserName,
				password=> $UserPassword)
	    || die "Net::SFTP not working\n";
    }else{
	$sftp = Net::SFTP->new( $IndTarget, 
				user=> $UserName,
				password=> $UserPassword,
				debug=>"true") 
	    || die "Net::SFTP not working\n";
    }    

    #-----------------------------+
    # TRANSFER  EVERY FILE IN @Tr |  
    #-----------------------------+
    for ($i=0; $i<=$MaxTr ; $i++)
    {
	# Transfer working but giving Couldn't fsetstat
	# error. It appears that is is transfering the
	# file even though it is throwing an error
	# This appears to be for files that already
	# exist on the target machine
	print "\tTransfering:\n\t".$Tr[$i][0]."\n"; 
	$sftp->put($Tr[$i][0], $Tr[$i][1]) ||
	    die "Transfer of $Tr[$i][0] not working\n";
    }
    
    print "Completed push to: ".$IndTarget."\n\n";
    
}

# STOP HERE TO JUST PUSH C116 FILES
#exit;

#-----------------------------+
# JLB10 COMPUTER              | 
#-----------------------------+
# REDEFINE THE TARGETS ARAY
@Targets = ("jlb10.genetics.uga.edu");


# COPY AND PASTE OF ABOVE CODE
# This could be written better to serve as a 
# subfunction taking variables:
#  @ Targets
#  @ Tr
#    - This would be a 2 d array
#      with a simple Src,Target relationship  
foreach my $IndTarget (@Targets)
{
    
    #-----------------------------+
    # CONNECT, USE DEBUG UNLESS   |
    # THE QUIET FLAG WAS USED     |
    #-----------------------------+
    # For some reason passing the password at the command line
    # throws an error (input must be 8 bytes long). Perhaps
    # this is an ANSI thing
    if ($Quiet)
    {
	print "Connecting to $IndTarget\n"; 
	$sftp = Net::SFTP->new( $IndTarget, 
				user=> $UserName,
				password=> $UserPassword)
	    || die "Net::SFTP not working\n";
    }else{
	$sftp = Net::SFTP->new( $IndTarget, 
				user=> $UserName,
				password=> $UserPassword,
				debug=>"true") 
	    || die "Net::SFTP not working\n";
    }    

    #-----------------------------+
    # TRANSFER  EVERY FILE IN @Tr |  
    #-----------------------------+
    for ($i=0; $i<=$MaxTr ; $i++)
    {
	print "\tTransfering:\n\t".$Tr[$i][0]."\n"; 
	$sftp->put($Tr[$i][0], $Tr[$i][2]) ||
	    die "Transfer of $Tr[$i][0] not working\n";
    }
    
    print "Completed push to: ".$IndTarget."\n\n";
    
}


#

exit;

