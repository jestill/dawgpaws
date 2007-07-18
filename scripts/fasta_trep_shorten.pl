#!/usr/bin/perl -w
#-----------------------------------------------------------+
#                                                           |
# fasta_short.pl - Give fasta file shorter names            |
#                                                           |
#-----------------------------------------------------------+
#                                                           |
#  AUTHOR: James C. Estill                                  |
# CONTACT: JamesEstill_@_gmail.com                          |
# STARTED: 07/17/2007                                       |
# UPDATED: 07/17/2007                                       |
#                                                           |
# DESCRIPTION:                                              |
#   Change headers in a fasta file to give shorter names    |
#   for trep files. This will reduce the fasta headers to   |
#   just the TREP ID.                                       |
#                                                           |
# USAGE:                                                    |
#  ShortFasta Infile.fasta Outfile.fasta                    |
#                                                           |
# VERSION:                                                  |
# $Id::                                                  $: |
#                                                           |
#-----------------------------------------------------------+

=head1 NAME

fasta_shorten.pl - Change headers in a fasta file to give shorter names.

=head1 VERSION

This documentation refers to program version 0.1

=head1 SYNOPSIS

 Usage:
  fasta_shorten.pl -i InDir -o OutDir

=head1 DESCRIPTION

This program will take all of the fasta files in a input directory 
and will shorten the name in the fasta header file. This is primarily
developed for the wheat project. The name used in the fasta header
is the name used by apollo as well as RepeatMaksker.

=head1 REQUIRED ARGUMENTS

=over 2

=item -i,--infile

Path of the input file.

=item -o,--outfile

Path of the output file.

=back

=head1 OPTIONS

=over 2

=item --usage

Short overview of how to use program from command line.

=item --help

Show program usage with summary of options.

=item --version

Show program version.

=item --man

Show the full program manual. This uses the perldoc command to print the 
POD documentation for the program.

=item -q,--quiet

Run the program with minimal output.

=back

=head1 DIAGNOSTICS

The list of error messages that can be generated,
explanation of the problem
one or more causes
suggested remedies
list exit status associated with each error

=head1 CONFIGURATION AND ENVIRONMENT

Names and locations of config files
environmental variables
or properties that can be set.

=head1 DEPENDENCIES

Other modules or software that the program is dependent on.

=head1 BUGS AND LIMITATIONS

Any known bugs and limitations will be listed here.

=head1 LICENSE

GNU LESSER GENERAL PUBLIC LICENSE

http://www.gnu.org/licenses/lgpl.html

=head1 AUTHOR

James C. Estill E<lt>JamesEstill at gmail.comE<gt>

=cut

package JPERL;

#-----------------------------+
# INCLUDES                    |
#-----------------------------+
use strict;
use Getopt::Long;

#-----------------------------+
# PROGRAM VARIABLES           |
#-----------------------------+
my $ver = "0.1";

#-----------------------------+
# VARIABLE SCOPE              |
#-----------------------------+
my $indir;                     # Input dir of fasta files to shorten
my $outdir;                    # Output dir for shortened fasta files
my $new_len = 20;             # New length of the header
my $head_new;                  # New, shorter header
my $test_len;                  # Test length of the header
my $cur_len;                   # Current length of the header

# Booleans
my $help = 0;
my $quiet = 0;
my $show_help = 0;
my $show_usage = 0;
my $show_man = 0;
my $show_version = 0;
my $test = 0;
my $uppercase = 0;            # Convert sequence strigs to uppercase

#-----------------------------+
# COMMAND LINE OPTIONS        |
#-----------------------------+
my $ok = GetOptions(
		    # Required options
		    "i|infile=s"  => \$indir,
                    "o|outfile=s" => \$outdir,
		    # Booleans
		    "uppercase"   => \$uppercase,
		    "usage"       => \$show_usage,
		    "version"     => \$show_version,
		    "man"         => \$show_man,
		    "h|help"      => \$show_help,
		    "test"        => \$test,
		    "q|quiet"     => \$quiet,);

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

#-----------------------------------------------------------+
# MAIN PROGRAM BODY                                         |
#-----------------------------------------------------------+
print "\n" unless $quiet;

#-----------------------------+
# CHECK FOR SLASH IN DIR      |
# VARIABLES                   |
#-----------------------------+
# If the indir does not end in a slash then append one
# TO DO: Allow for backslash
unless ($indir =~ /\/$/ ) {
    $indir = $indir."/";
}

unless ($outdir =~ /\/$/ ) {
    $outdir = $outdir."/";
}

#-----------------------------+
# Get the FASTA files from the|
# directory provided by the   |
# var $indir                  |
#-----------------------------+
opendir( DIR, $indir ) || 
    die "Can't open directory:\n$indir"; 
my @fasta_files = grep /\.fasta$|\.fa$/, readdir DIR ;
closedir( DIR );

#-----------------------------+
# CREATE THE OUT DIR          |
# IF IT DOES NOT EXIST        |
#-----------------------------+
print "Creating output dir ...\n" unless $quiet;
unless (-e $outdir) {
    mkdir $outdir ||
	die "Could not create the output directory:\n$outdir";
}


#-----------------------------+
# PROCESS EACH INDIVIDUAL     |
# FILE IN THE DIRECTORY       |
#-----------------------------+
for my $ind_file (@fasta_files) {

    my $infile = $indir.$ind_file;
    my $outfile = $outdir.$ind_file;

    open (IN, $infile) ||
	die "Can not open infile:\n$infile\n";
    
    open (OUT, ">".$outfile) ||
	die "Can not open outfile:\n$outfile\n";

    while (<IN>) {
	
	chomp;
	# If the first characeter is a >
	if (/^\>/) {
	    #-----------------------------+
	    # FASTA HEADER PROCESSING     |
	    #-----------------------------+
	    
	    print " processing $_\n";

	    my @header_data = split;
	    print "Full_ID:".$header_data[0]."\n";
	    
	    my @full_trep_id = split(/\|/,$header_data[0]);
	    print "TREP:".$full_trep_id[2]."\n";

	    print OUT ">".$full_trep_id[2]."\n";

	}
	
	else {
	    #-----------------------------+
	    # SEQUENCE STRING PROCESSING  |
	    #-----------------------------+
	    # print the string to the outfile unchanged if this is
	    # not a fasta header line
	    # T
	    #if ($)
	    if ($uppercase) {
		# Convert the sequence string to uppercase
		# I will try this as a bare regexp since 
		# we are workign with $_
		tr/a-z/A-Z/;
		print OUT "$_\n";

	    }
	    else {
		# If not converting the sequence string to uppercase
		# print the output file unchanged.
		print OUT "$_\n";
	    }
	}


    } # End of while IN

} # End of for each file in the directory

exit;

#-----------------------------------------------------------+ 
# SUBFUNCTIONS                                              |
#-----------------------------------------------------------+

sub print_help {

    # Print requested help or exit.
    # Options are to just print the full 
    my ($opt) = @_;
    
    my $usage = "USAGE:\n". 
	"fasta_shorten.pl -i InFile -o OutFile";
    
    my $args = "REQUIRED ARGUMENTS:\n".
	"  --indir        # Path to the input directory\n".
	"  --outdir       # Path to the output directory\n".
	"                 # Will be created if needed".
	"\n".
	"OPTIONS:\n".
	"  --length       # Length to shorten the name to\n".
	"                 # default is 20\n".
	"  --uppercase    # All sequence strings are converted to uppercase\n".
	"  --test         # Run a program test\n".
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

STARTED: 07/17/2007

UPDATED: 07/17/2007

VERSION: $Id: script_template.pl 68 2007-07-15 00:25:18Z JamesEstill $

=cut

#-----------------------------------------------------------+
# HISTORY                                                   |
#-----------------------------------------------------------+
#
# 07/17/2007
# - Program started.
