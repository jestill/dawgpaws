#!/usr/bin/perl -w
#-----------------------------------------------------------+
#                                                           |
# merge_ammned.pl - Merge recent annotation an ammendments  |
#                                                           |
#-----------------------------------------------------------+
#                                                           |
#  AUTHOR: James C. Estill                                  |
# CONTACT: JamesEstill_at_gmail.com                         |
# STARTED: 11/09/2007                                       |
# UPDATED: 11/09/2007                                       |
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
my $ammend_file;               # Path to game xml file with ammendments

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

my $usage = "merge_ammend -i orig_file -a ammend_file -o outfile\n";

#-----------------------------+
# COMMAND LINE OPTIONS        |
#-----------------------------+
my $ok = GetOptions(# Required options
                    "i|infile=s"   => \$infile,
		    "a|ammend=s"   => \$ammend_file,
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
		    "h|help"       => \$show_help,
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
if ( (!$infile) || (!$ammend_file) || (!$outfile) ) {
    print "No infile path provided\n" if !$infile;
    print "No outfile path provided" if !$outfile;
    print "No ammended file provided" if !$bac_name;
    print "$usage\n\n";
    exit;
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

unless (-e $ammend_file) {
    die "Ammend could not be found:\n$ammend_file\n";
}


if ($logfile) {
    open (LOG,">".$logfile ) 
	|| die "Can not open error log:\n$logfile\n";
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
# OPEN THE ADD DATA DIALOG    |
#-----------------------------+

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
	    "BAC:\t$bac_name\n";
	die "$err_msg\n";
    }   
}
else {
    die "Could not get apollo window viewable for $bac_name:\n";
}

#-----------------------------+
# ADD THE AMMENDED DATA       |
#-----------------------------+
if (WaitWindowViewable( 'Apollo: adding data' ,undef, 20)) {
    
    SendKeys('add{TAB}^(a){BAC}'.$ammend_file.'{TAB 2}{SPA}') ||
	die "Could not send keys to open the file:\n$infile";
    
}


#-----------------------------+
# OPEN SAVE DATA DIALOG       |
#-----------------------------+
sleep 5;
if (WaitWindowViewable($bac_name, undef, 120)) {

    my $win_focus = GetWindowName(GetInputFocus());
    
    if ($win_focus =~ $bac_name) {
#	SendKeys ('%(f){DOW SPA}^c') ||
	SendKeys ('^(a)') ||
	    die "Could not open the dialog window";
    } 
    else {
	my $err_msg = "I expected the Apollo window $bac_name\n".
	    "to have the focus but the focus is currently on:\n".
	    "$win_focus\n".
	    "BAC:\t$bac_name\n";
	die "$err_msg\n";
    }   
    
}

#-----------------------------+
# SAVE THE MERGE FILE         |
#-----------------------------+
if (WaitWindowViewable( 'Write data' ,undef, 20)) {
    
    SendKeys('{TAB}^(a){BAC}'.$outfile.'{TAB 4}{SPA}') ||
	die "Could not send keys to save the file:\n$outfile";
    
}


# CAN QUIT HERE

exit;

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
