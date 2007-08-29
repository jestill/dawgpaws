#!/usr/bin/perl -w
#-----------------------------------------------------------+
#                                                           |
# fetch_tenest.pl - Fetch TE Nest LTR and SVG files         |
#                                                           |
#-----------------------------------------------------------+
#                                                           |
#  AUTHOR: James C. Estill                                  |
# CONTACT: JamesEstill_@_gmail.com                          |
# STARTED: 08/29/2007                                       |
# UPDATED: 08/29/2007                                       |
#                                                           |
# DESCRIPTION:                                              |
#  Workflow to fetch the tenest program results and covert  |
#  the output to gff format for use in the DAWG-PAWS        |
#  pipeline.                                                |
#                                                           |
# VERSION:                                                  |
# $Id:: fetch_tenest.pl 85 2007-08-29 14:29:27Z JamesEst#$: |
#                                                           |
#-----------------------------------------------------------+

=head1 NAME

fetch_tenest.pl - Fetch TE Nest LTR and SVG files

=head1 VERSION

This documentation refers to program version 1.0

=head1 SYNOPSIS

  USAGE:
    fetch_tenest.pl -i SeqList.txt -o OutDir

=head1 DESCRIPTION

Given a tab delimited input file of theFetches the TE Nest results using wget.
The results are converted to 

=head1 REQUIRED ARGUMENTS

=over 2

=item -i,--infile

Path of the input file. This is a tab delimited input file that must be
in the following format:

    HEX0350E24     1188334554_18873
    HEX0358K19     1188334578_19208
    HEX0411K13     1188334604_19821

Where the first column is the name of the sequence that was analyzed
by TE Nest and the second column is the id of the TE Nest job. This job 
id is returned by TE Nest when the job is initially submitted.

=item -o,--outdir

Path of the output directory. This is the base dir that the analysis 
results are stored in. 

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

package DAWGPAWS;

#-----------------------------+
# INCLUDES                    |
#-----------------------------+
use strict;
use Getopt::Long;

#-----------------------------+
# PROGRAM VARIABLES           |
#-----------------------------+
my ($VERSION) = q$Rev$ =~ /(\d+)/;

#-----------------------------+
# VARIABLE SCOPE              |
#-----------------------------+

# my 

my $infile;
my $outfile;

# Booleans
my $help = 0;
my $quiet = 0;
my $show_help = 0;
my $show_usage = 0;
my $show_man = 0;
my $show_version = 0;

#-----------------------------+
# COMMAND LINE OPTIONS        |
#-----------------------------+
my $ok = GetOptions("i|infile=s"  => \$infile,
                    "o|outfile=s" => \$outfile,
		    "usage"       => \$show_usage,
		    "version"     => \$show_version,
		    "man"         => \$show_man,
		    "h|help"      => \$show_help,
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
    print "\n$0:\nVersion: $VERSION\n\n";
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

#-----------------------------+
# LOAD DATA TO FETCH FROM     |
# TAB DELIM TEXT FILE         |
#-----------------------------+
open (INFILE, "<$infile")
    || die "Can not open the input file:\n$infile\n";

while (<INFILE>) {

}

close INFILE;

#-----------------------------+
# FETCH DATA TO APPROPRIATE
# DIR, 

# Create tenest dir if it does not exist

# COPY FILES TO NEW NAME
# Use copy instead of move to keep track of what the original
# file id was

# 

# Convert te_nest output to gff format

# copy gff file to gff dir


exit;

#-----------------------------------------------------------+ 
# SUBFUNCTIONS                                              |
#-----------------------------------------------------------+

sub print_help {

    # Print requested help or exit.
    # Options are to just print the full 
    my ($opt) = @_;

    my $usage = "USAGE:\n". 
	"fetch_tenest.pl -i InFile -o OutFile";
    my $args = "REQUIRED ARGUMENTS:\n".
	"  --infile       # Path to the input file\n".
	"  --outfile      # Path to the output file\n".
	"\n".
	"OPTIONS::\n".
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

STARTED:

UPDATED:

VERSION: $Id: fetch_tenest.pl 85 2007-08-29 14:29:27Z JamesEstill $

=cut

#-----------------------------------------------------------+
# HISTORY                                                   |
#-----------------------------------------------------------+
#
# 08/29/2007
# - Program started
