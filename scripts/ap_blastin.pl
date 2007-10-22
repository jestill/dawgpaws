#!/usr/bin/perl -w
#-----------------------------------------------------------+
#                                                           |
# blast2ap.pl - Load BLAST alignments to Apollo             |
#                                                           |
#-----------------------------------------------------------+
#                                                           |
#  AUTHOR: James C. Estill                                  |
# CONTACT: JamesEstill_at_gmail.com                         |
# STARTED: 06/26/2007                                       |
# UPDATED: 07/02/2007                                       |
#                                                           |
# DESCRIPTION:                                              | 
#  Import full BLAST output information into an existing    |
#  game xml file using the Apollo genome annotation         |
#  curation program.                                        |
#                                                           |
# LICENSE:                                                  |
#  GNU General Public License, Version 3                    |
#  http://www.gnu.org/licenses/gpl.html                     |  
#                                                           |
#-----------------------------------------------------------+
#
# WARNING: THIS SOFTWARE IS CURRENTLY NOT A FUNCTIONAL 
#
# TODO: REGEXP Test to see if multiple windows are open that could
#       take the focus HEX---
#

=head1 NAME

ap_blastin.pl - Add BLAST alignments to a game.xml file.

=head1 SYNOPSIS

  Usage: ap_blastin.pl -i Infile.game.xml -b BlastDir -o Outfile.game.xml

=head1 DESCRIPTION

Import full BLAST output information into an existing
game xml file using the Apollo genome annotation
curation program. This uses the X11::GUITest PERL module
to run the Apollo program as if the user were entering commands.
This is a very ugly way to get this done, but I can't find
something better at the moment.

=head1 ARGUMENTS

=over

=item -i, --infile

Full path to the intput game.xml file that the blast data will
be added to.

=item -b, --blastdir

Full path to the BLAST directory that contains the BLAST alignment 
output files to parse.

=item -o, --outfile

Full path to the outfile.

=item -q, --quiet

Run the program in quiet mode. Minimal information is printed to STDOUT.

=back

=head1 AUTHOR

James C. Estill E<lt>JamesEstill at gmail.comE<gt>

=cut

package DAWGPAWS;

my $time_start = time;

#-----------------------------+
# INCLUDES                    |
#-----------------------------+
use Cwd;                       # Get the current working directory
use Getopt::Long;
use Carp;                      # Access to the croak function
use strict;
#use X11::GUITest;
use X11::GUITest qw/
    StartApp
    WaitWindowViewable
    SendKeys
    GetInputFocus
    SetInputFocus
    GetWindowName
    WaitWindowClose
    FindWindowLike
    /;


#-----------------------------+
# HARD CODED VARIABLES        |
#-----------------------------+
my $ver = "1.0";               # Program version

#-----------------------------+
# VARIABLE SCOPE              |
#-----------------------------+
# Deafult values
my $ap_path = '/home/jestill/Apps/Apollo_1.6.5/apollo/bin/apollo';

# Command line variables
my $logfile;
my $add_dialog_open;           # boolean, is the add_dialog_open
my $blastdir;                  # Directory containing the blast output files
my $infile;                    # Full path to the input file 
                               # This should be a valid game.xml file
my $outfile;                   # Full path to game.xml that you would
                               # like to create. If this file already
                               # exists, it will be overwritten.
my $bac_name;                  # Name of the BAC, this will be the
                               # expected name of the Apollo window

# BOOLEANS
my $verbose = 0;
my $show_help = 0;
my $show_usage = 0;
my $show_man = 0;
my $show_version = 0;
my $quiet = 0;                 # Boolean: Run the program in quiet mode

# COUNTERS AND INDICES
my $blast_file_count = 0;      # Counter for blast files 

# APOLLO WINDOW IDS
my $win_focus_id;              # The id number of the window that has focus
my $win_focus_prev;            # Previous window id
my @winid_load_data;           # window ids of the load data window
my @winid_main;                # window ids of the main Apollo window
                               # This is the window that shows alignment
my @windid_add_data;           # window ids of the add data Apollo window

my $usage = "blast2ap.pl -i infile -b blastdir -o outfile -n bac_name\n".
    " \n\n".
    "blast2ap.pl --help\n".
    "To print the full help documentation.";

#-----------------------------+
# COMMAND LINE OPTIONS        |
#-----------------------------+
my $ok = GetOptions(# Required options
		    "b|blastdir=s" => \$blastdir,
                    "i|infile=s"   => \$infile,
		    "o|outfile=s"  => \$outfile,
		    "n|name=s"     => \$bac_name,
		    # Additional options
		    "ap-path=s"    => \$ap_path,
		    "logfile=s"    => \$logfile,
		    # Booleans
		    "verbose"      => \$verbose,
		    "usage"        => \$show_usage,
		    "version"      => \$show_version,
		    "man"          => \$show_man,
		    "q|quiet"      => \$quiet,
		    "h|help"       => \$show_help
		    );


#-----------------------------+
# SHOW REQUESTED HELP         |
#-----------------------------+
if ($show_usage) {
    print_help("");
}

if ($show_help || (!$ok) ) {
    print_help("full");
}

if ($show_version) {
    print "\n$0:\nVersion: $ver\n\n";
    exit;
}

if ($show_man) {
    # User perldoc to generate the man documentation.
    system("perldoc $0");
    exit($ok ? 0 : 2);
}

print "$0 Starting...\n" if $verbose;

#-----------------------------+
# PRINT USAGE WHEN REQUIRED   |
# CMD LIEN VARAIBLES MISSING  |
#-----------------------------+
# If the Input, blast, and outfile are not specified 
# at the command line show the user the program usage statement
if ( (!$infile) || (!$blastdir) || (!$outfile) || (!$bac_name)) {
    print "No infile path provided\n" if !$infile;
    print "No blast directory provided\n" if !$blastdir;
    print "No outfile path provided" if !$outfile;
    print "No BAC name provided" if !$bac_name;
    print "$usage\n\n";
    exit;
}


#-----------------------------+
# CHECK FOR SLASH IN DIR      |
# VARIABLES                   |
#-----------------------------+
# If the indir does not end in a slash then append one
# TO DO: Allow for backslash
unless ($blastdir =~ /\/$/ ) {
    $blastdir = $blastdir."/";
}

#-----------------------------+
# SET FULL PATHS AND CHECK    |
# FOR FILE EXISTENCE          |
#-----------------------------+ 
#my $infile = $work_dir.$infile;
#my $blastdir = $work_dir.$blastdir;

unless (-e $infile) {
    die "Infile could not be found:\n$infile\n";
}

unless (-e $blastdir) {
    die "blast dir could not be found:\b$blastdir\n";
}

if ($logfile) {
    open (LOG,">".$logfile ) 
	|| die "Can not open error log:\n$logfile\n";
}

#-----------------------------+
# LOAD BLAST FILES TO ARRAY   |
#-----------------------------+
opendir DIR, $blastdir ||
    die "Could not open the dir:\n$blastdir\n";
my @blast_files = grep /\.blx$/ || /\.blo$/ || /\.bln$/, readdir DIR;
my $num_blast = @blast_files;
@blast_files = sort(@blast_files);
print "$num_blast blast files found\n" if $verbose;

# Print the blast dir info
my $i=0;
foreach my $ind_blast_file (@blast_files) {
    print "\t$ind_blast_file\t" if $verbose;
    print "$blast_files[$i]\n" if $verbose;
    $i++;

    my $in_blast_path = $blastdir.$ind_blast_file;
    unless (-e $in_blast_path) {
	die "Could not find blast file:\n$in_blast_path\n";
    }

    # Can check if the blast files have any hits here,
    # and drop them from the array to input if it is empty
    # this could just be loaded into a new array
}

#-----------------------------+
# START APOLLO                |
#-----------------------------+
StartApp($ap_path);

#-----------------------------+
# OPEN INPUT GAME XML FILE    |
#-----------------------------+
# Wait for the window to be available before proceeding
# and only continue if the expected window has the focus
if (WaitWindowViewable('Apollo: load data', undef, 10)) {
   
    @winid_load_data = FindWindowLike ('Apollo: load data', undef,) ||
	die "Could not find the load data window.\n";
    
    # Get the name of the window that has the current focus
    my $win_focus = GetWindowName(GetInputFocus());
    
    print "Focus win Name: $win_focus\n" if $verbose;       
    
    if ($win_focus =~ 'Apollo: load data') { 
	SendKeys('add{TAB}^(a){BAC}'.$infile.'{TAB 2}{SPA}') ||
	    die "Could not send keys to open the file:\n$infile";
    } 
    else {
	die "I expected Apollo: load data to be in focus for\n$infile";
    } 
}


#-----------------------------+
# ADD BLAST ALIGNMENTS        |
#-----------------------------+

foreach my $ind_blast_file (@blast_files) {

    # SET THE BLAST PATH
    my $in_blast_path = $blastdir.$ind_blast_file;
    $blast_file_count++;

    print LOG "\n\nProcessing:\n$in_blast_path\n" if $logfile;
    $add_dialog_open = 0;

    print "PROCESSING:\n\t$in_blast_path\n" if $verbose;
    print "NUM\t$blast_file_count\n" if $verbose;
    print "TOTAL\t\n$num_blast\n" if $verbose;

    #-----------------------------+
    # OPEN THE IMPORT BLAST       |
    # COMMAND WINDOW              |
    #-----------------------------+
    # Wait for the window to be viewable before continuing
    # then make sure that it has the focus before trying to send 
    # the commands to open the import blast dialog. This is required
    # because the window may be viewable even if it does not have the focus
    #if (WaitWindowViewable($bac_name, undef, 10)) {
    if (WaitWindowViewable($bac_name, undef, 120)) {

	my $win_focus = GetWindowName(GetInputFocus());
	
	if ($win_focus =~ $bac_name) {
	    SendKeys ('%(f){DOW SPA}^c') ||
		die "Could not open the dialog window";
	} 
	else {
	    my $err_msg = "I expected the Apollo window $bac_name\n".
		"to have the focus but the focus is currently on:\n".
		"$win_focus\n".
		"BAC:\t$bac_name\n".
		"BLAST:\t$ind_blast_file\n".
		"FNUM:\t$blast_file_count\n".
		"FCOUNT:\t$num_blast\n";
	    die "$err_msg\n";
	}

    } else {
	die "Could not get apollo window viewable for $bac_name:\n";
    }


    #-----------------------------+
    # IMPORT THE BLAST REPORT     |
    #-----------------------------+
    # Change the following to something more straightforward
    if (WaitWindowViewable( 'Apollo: adding data' ,undef, 20)) {
	

	$win_focus_id = GetInputFocus();
	my $win_focus = GetWindowName($win_focus_id);
	print "\n\nFOCUS ID:\t$win_focus_id\n" if $verbose;
	print "FOCUS NAME:\t$win_focus\n" if $verbose;


	#-----------------------------+
	# CHECK FOCUS ID              |
	#-----------------------------+
	# WOrkign change is here
	# IF THE WINDOW FOCUS ID HAS NOT CHANGED
	if ($win_focus_prev == $win_focus_id) {
	    print "The window id did not change.\n";
	    print "PREV:\t$win_focus_prev\n";
	    print "NOW:\t$win_focus_id\n";
	    print "NAME:\t$ind_blast_file\n";
	    print "File $blast_file_count of $num_blast\n";
	    print "";
	    
	    # ALT F4 Will close the current window in apollo
	    # If the window currently in focus is the main widow
	    # this will also be closed.
	    SendKeys('%(F4)');
	    
	    exit;
	} 
	else {
	    print "The focus id has changed.\n" if $verbose;
	    print "PREV:\t$win_focus_prev\n" if $verbose;
	    print "NOW:\t$win_focus_id\n" if $verbose;

	    print LOG "The focus id has changed.\n" if $logfile;
	    print LOG "PREV:\t$win_focus_prev\n" if $logfile;
	    print LOG "NOW:\t$win_focus_id\n" if $logfile;
	}	
#	exit;
	
	
	if ($win_focus =~ 'Apollo: adding data') {
	    # Reset window to one of the d options
	    
	    SendKeys('d') ||
		die "ERROR: Problem with SendKeys d\n";
	    sleep 1;
	    #sleep 0.3;
	    
	    # Select chado xml first
	    SendKeys('c') ||
		die "ERROR: Problem with SendKeys c\n";
	    sleep 1;
	    #sleep 0.3;
	    
	    # Select computational results
	    SendKeys('c') ||
		die "ERROR: Problem with SendKeys c\n";
	    sleep 1;
	    #sleep 0.3;
	    
	    # The following send keys can probably
	    # be sent as a package
            # Press TAB 2 times
	    SendKeys('{TAB 2}') ||
		die "ERROR: Problem with SendKeys TAB 2\n"; 

	    # Ctrl-a to select all
	    SendKeys('^(a)') ||
		die "ERROR: Problem with SendKeys Ctrl-a\n";
	    
	    # Backspace to delete existing information
	    SendKeys('{BAC}') ||
		die "ERROR: Problem with SendKeys BACKSPACE\n";   

	    # Type in the path to the blast file
	    SendKeys($in_blast_path) ||
		die "ERROR: Problem with SendKeys BLAST path\n"; 

	    # Tab to the OK button
	    SendKeys('{TAB 27}') ||
		die "ERROR: Problem with SendKeys TAB 27\n";  
	    
	    # Spacebar to select OK
	    SendKeys('{SPA}') ||
		die "ERROR: Problem with SendKeys SPACE";
 
	} # End of if adding data dialog has focus

	
    } else {
	die "Could not open Apollo: add dialog.\n";
    }
    
    # Reset focus ID
    $win_focus_prev = $win_focus_id; 

#    exit;
    
#    # If the import data window is open, wait for it to close
    if (GetWindowName(GetInputFocus()) =~ 'Apollo: adding data') {
	# To close the window we will need to use the numeric 
	# argument from above, we will wait 30 seconds
	# before moving on
	#$win_focus_id 
	print "Trying to close the adding data window\n";
	if (WaitWindowClose ($win_focus_id, 30)) {
	    print "The add data window is now closed.\n" if $verbose;
	}
	else {
	    die "The add data window did not close.\n";
	}
    } 
    else {
	print "The adding data widow appears to be closed.\n" if $verbose;
    }

    # The alternative to the above is to always sleep
#    sleep 5;

#    # ADDED A WAIT HERE FOR THE WINDOW TO CLOSE
# The WaitWindowClose expects a numeric arugment
# This probably refers to the 
#
#    WaitWindowClose ( $window_focus , 10) || 
#	die "Could not wait for the window close\n"; 
    
} # End of for each ind blast file




print "$0 complete.\n";

my $time_end = time;
my $time_total = $time_end - $time_start;
print "\n\n\n\n\nTOTAL TIME:$time_total\n";

exit;

#-----------------------------------------------------------+
# SUBFUNCTIONS                                              |
#-----------------------------------------------------------+

sub print_help {

    # Print requested help or exit.
    # Options are to just print the full 
    my ($opt) = @_;

    my $usage = "USAGE:\n". 
	"ap_blastin.pl -i Infile.game.xml -b BlastDir -o Outfile.game.xml";
	#"ap_blastin.pl -i InFile -o OutFile";
    my $args = "REQUIRED ARGUMENTS:\n".
	"  --infile       # Path to the input game.xml file\n".
	"  --outfile      # Path to the output game.xml file\n".
	"  --blastdir     # Directory containing BLAST files to import\n".
	"\n".
	"ADDITIONAL OPTIONS:\n".
	"".
	"  --version      # Show the program version\n".     
	"  --usage        # Show program usage\n".
	"  --help         # Show this help message\n".
	"  --man          # Open full program manual\n".
	"  --quiet        # Run program with minimal output\n";
	
    if ($opt =~ "full") {
	print "\n$usage\n\n";
	print "$args\n\n";
    }
    else {
	print "\n$usage\n\n";
    }
    
    exit;
}

=head1 HISTORY

Started: 06/26/2007

Updated: 07/27/2007

=cut



#-----------------------------------------------------------+
# TO DO LIST                                                |
#-----------------------------------------------------------+
# Working with the following later to make sure
# that an apollo window is not already currently in play
#    my $num_load_data = @winid_load_data;
#    # If the number of load data windows open is just one
#    print $num_load_data." load data windows available \n";
#
#    if ($num_load_data == 1) {
#	print "The winid is ".$winid_load_data[0]."\n";
#    } 


#-----------------------------------------------------------+
# HISTORY                                                   |
#-----------------------------------------------------------+
#
# 06/26/2007 
#  -Started program
#
# 06/27/2007
#  - Working with getting this to run without crashing
#    or skipping over records
# 
# 06/28/2007
# - Still working on making this run smoothly  
# - Added tests to make sure that the expected window acutall
#   has the focus before proceeding
#
# 06/29/2007
# - Continuing to add error checks to make this run smoothly
#
# 07/02/2007
# - Working to make this die on errors
