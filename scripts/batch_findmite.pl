#!/usr/bin/perl -w
#-----------------------------------------------------------+
#                                                           |
# batch_findmite.pl - Batch run the findmite program        |
#                                                           |
#-----------------------------------------------------------+
#                                                           |
#  AUTHOR: James C. Estill                                  |
# CONTACT: JamesEstill_@_gmail.com                          |
# STARTED: 08/30/2007                                       |
# UPDATED: 08/30/2007                                       |
#                                                           |
# DESCRIPTION:                                              |
#  Short Program Description                                |
#                                                           |
# USAGE:                                                    |
#  ShortFasta Infile.fasta Outfile.fasta                    |
#-----------------------------------------------------------+

=head1 NAME

batch_findmite.pl - Run the findmite program in batch mode.

=head1 VERSION

This documentation refers to $Rev$

=head1 SYNOPSIS

  USAGE:
    batch_findminte.pl -i InDir -o OutDir

=head1 DESCRIPTION



=head1 COMMAND LINE ARGUMENTS

=head2 Required Argumens

=over 2

=item -i,--indir

Path of the input directory. This is a directory that contains the fasta
files to process.

=item -o,--outdir

Path of the base output directory.

=back

=head2 Additional Options

=over 2

=item -q,--quiet

Run the program with minimal output. Does not require user interaction.

=item --verbose

Run the program with maximal output.

=back

=head2 Additional Information

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
my ($VERSION) = q$Rev$ =~ /(\d+)/;

#-----------------------------+
# VARIABLE SCOPE              |
#-----------------------------+
my $indir;
my $outdir;

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
my $ok = GetOptions(# Required Arguments
		    "i|indir=s"   => \$indir,
                    "o|outdir=s"  => \$outdir,
		    # Additional options
		    "q|quiet"     => \$quiet,
		    # Additional information
		    "usage"       => \$show_usage,
		    "version"     => \$show_version,
		    "man"         => \$show_man,
		    "h|help"      => \$show_help,);


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


exit;

#-----------------------------------------------------------+ 
# SUBFUNCTIONS                                              |
#-----------------------------------------------------------+

sub print_help {

    # Print requested help or exit.
    # Options are to just print the full 
    my ($opt) = @_;

    my $usage = "USAGE:\n". 
	"batch_findmite.pl -i InDir -o OutDir";
    my $args = "REQUIRED ARGUMENTS:\n".
	"  --infdir       # Directory of fasta files to process\n".
	"  --outdir       # Path to the output direcotry\n".
	"\n".
	"OPTIONS:\n".
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

STARTED: 08/30/2007

UPDATED: 08/30/2007

VERSION: $Id: batch_findmite.pl 95 2007-08-30 15:32:36Z JamesEstill $

=cut

#-----------------------------------------------------------+
# HISTORY                                                   |
#-----------------------------------------------------------+
#
# 08/30/2007
# - Program started
