#!/usr/bin/perl -w
#-----------------------------------------------------------+
#                                                           |
# cnv_gff2game.pl - Convert a gff file to game xml          |
#                                                           |
#-----------------------------------------------------------+
#                                                           |
#  AUTHOR: James C. Estill                                  |
# CONTACT: JamesEstill at gmail.com 
# STARTED: 02/09/2007                                       |
# UPDATED: 08/07/2007                                       |
# DESCRIPTION:                                              |
#  Converts gff data tracks from gff format to the game     |
#  xml format for use in the Apollo Genome Annotation       |
#  Curation program.                                        |
#                                                           |
#-----------------------------------------------------------+ 

=head1 NAME

cnv_gff2game.pl - Convert a gff file to game xml

=head1 VERSION

This documentation refers to program version 1.0

=head1 SYNOPSIS

 Usage:
  cnv_gff2game.pl -i InFile.fasta -g GffFile.gff -o OutFile.game.xml

=head1 DESCRIPTION

This is what the program does

=head1 REQUIRED ARGUMENTS

=over 2

=item -i,--infile

Path of the input fasta filefile.

=item -g, --gff

Path of the input gff file

=item -o,--outfile

Path of the output game xml file.

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

=over 2 

=item Apollo

Requires that apollo be installed on the local machine
since apollo is being used as the engine to do the
conversion between formats.

=back

=head1 BUGS AND LIMITATIONS

Any known bugs and limitations will be listed here.

=head1 LICENSE

GNU LESSER GENERAL PUBLIC LICENSE

http://www.gnu.org/licenses/lgpl.html

=head1 AUTHOR

James C. Estill E<lt>JamesEstill at gmail.comE<gt>

=cut


#-----------------------------+
# INCLUDES                    |
#-----------------------------+
use strict;
use Getopt::Long;

#-----------------------------+
# VARIABLES
#-----------------------------+

# CONSTANTS
my $VERSION = "1.0";

my $infile_fasta;
my $infile_gff;
my $outfile;

# Options with default values
my $ap_path = "/home/jestill/Apps/Apollo_1.6.5/apollo/bin/apollo";

# Booleans
my $quiet = 0;
my $show_help = 0;
my $show_usage = 0;
my $show_man = 0;
my $show_version = 0;

#-----------------------------+
# COMMAND LINE OPTIONS        |
#-----------------------------+
my $ok = GetOptions(# REQUIRED
		    "i|infile=s"  => \$infile_fasta,
		    "g|gff=s"     => \$infile_gff,
                    "o|outfile=s" => \$outfile,
		    # ADDITIONAL OPTIONS
		    "ap-path=s"   => \$ap_path,
		    # BOOLEANS
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

#-----------------------------+
# CHECK FOR REQUIRED VARS     |
#-----------------------------+
if ( (!$infile_fasta) || (!$infile_gff) || (!$outfile) ) {
    print_help("full");
}

# Convert each the gff file to game xml using the apollo_convert command
&apollo_convert ($infile_gff, "gff", $outfile, 
		 "game", $infile_fasta, "NULL");

exit;

#-----------------------------------------------------------+
# SUBFUNCTIONS                                              |
#-----------------------------------------------------------+

sub apollo_convert
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

    my ($InFile,$InForm,$OutFile,$OutForm,$SeqFile,$DbPass) = @_;

    #InFile = $_[0];        # Input file path
    #InForm = $_[1];        # Output file format:
    #                         # game|gff|gb|chadoxml|backup
    #$OutFile = $_[2];       # Output file path
    #$OutForm = $_[3];       # Ouput file foramt
    #                           # chadoDB|game|chadoxml|genbank|gff|backup
    #$SeqFile = $_[4];       # The path of the sequence file
    #                           # This is only required for GFF foramt files
    #                           # When not required this can be passed as na
    #$DbPass = $_[5];        # Database password for logging on to the 
    #                           # chado database for reading or writing.
    my $ApCmd;

    # Set the base command line. More may need to be added for different
    # formats. For example, GFF in requires a sequence file and the CHADO
    # format will require a database password.
    $ApCmd = $ap_path." -i ".$InForm." -f ".$InFile.
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

STARTED: 02/09/2007

UPDATED: 08/07/2007

VERSION: $Id: script_template.pl 68 2007-07-15 00:25:18Z JamesEstill $

=cut

#-----------------------------------------------------------+
# HISTORY                                                   |
#-----------------------------------------------------------+
#
