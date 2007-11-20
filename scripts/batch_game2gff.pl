#!/usr/bin/perl -w
#-----------------------------------------------------------+
#                                                           |
# batch_game2gff.pl - Batch convert gamexml to gff format   |
#                                                           |
#-----------------------------------------------------------+
#                                                           |
#  AUTHOR: James C. Estill                                  |
# CONTACT: JamesEstill_@_gmail.com                          |
# STARTED: 11/01/2007                                       |
# UPDATED: 11/01/2007                                       |
#                                                           |
# DESCRIPTION:                                              |
#  Convert game.xml annotated feature file format to        |
#  Goal is to only export the curated results.              | 
#                                                           |
# REQUIREMENTS:                                             |
#  This requires the apollo genome annotation curation      |
#  program.                                                 |
#                                                           |
# VERSION: $Rev$                                            |
#                                                           |
#                                                           |
# LICENSE:                                                  |
#  GNU General Public License, Version 3                    |
#  http://www.gnu.org/licenses/gpl.html                     |  
#                                                           |
#-----------------------------------------------------------+

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

my $indir;
my $outdir;

# Booleans
my $quiet = 0;
my $verbose = 0;
my $show_help = 0;
my $show_usage = 0;
my $show_man = 0;
my $show_version = 0;

# Options with default values
my $ap_path = "/home/jestill/Apps/Apollo_1.6.5/apollo/bin/apollo";

#-----------------------------+
# COMMAND LINE OPTIONS        |
#-----------------------------+
my $ok = GetOptions(# REQUIRED OPTIONS
		    "i|infile=s"  => \$indir,
                    "o|outfile=s" => \$outdir,
		    # ADDITIONAL OPTIONS
		    "ap-path=s"   => \$ap_path,
		    "q|quiet"     => \$quiet,
		    "verbose"     => \$verbose,
		    # ADDITIONAL INFORMATION
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

# Throw error when required options not in command line
if ( (!$indir) || (!$outdir) ) {
    print "\a";
    print "ERROR: Input directory must be specified\n" if !$indir;
    print "ERROR: Output directory must be specified\n" if !$outdir;
    print_help("full");

}

#-----------------------------+
# MAIN PROGRAM BODY           |
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
# Get the seq files from the  |
# directory provided by the   |
# var $indir                  |
#-----------------------------+
opendir( DIR, $indir ) || 
    die "Can't open directory:\n$indir"; 
my @seq_files = grep /\.game.xml$|\.xml$/, readdir DIR ;
closedir( DIR );

my $count_files = @seq_files;


if ($count_files == 0) {
    print "\a";
    print "No game.xml files were found in the intput directory\n";
    print "$indir\n";
}

print "$count_files to process\n";

for my $ind_file (@seq_files)
{
    print "Processing: $ind_file\n";

    my $in_seq_path = $indir.$ind_file;
    my $out_gff_path = $outdir.$ind_file.".gff";
    my $out_fasta_path = $outdir.$ind_file.".fasta";
 
    #  Convert each the game xml to gff using the apollo_convert command
    &apollo_convert ($in_seq_path, "game", $out_gff_path, 
		     "gff", "NULL", "NULL");
    

}


exit;

#-----------------------------------------------------------+ 
# SUBFUNCTIONS                                              |
#-----------------------------------------------------------+

sub print_help {

    # Print requested help or exit.
    # Options are to just print the full 
    my ($opt) = @_;

    my $usage = "USAGE:\n". 
	"MyProg.pl -i InFile -o OutFile";
    my $args = "REQUIRED ARGUMENTS:\n".
	"  --indir       # Path to the input file\n".
	"  --outdir      # Path to the output file\n".
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


=head1 NAME

batch_game2gff.pl - Convert game.xml annotations to gff format

=head1 VERSION

This documentation refers to program version $Rev$

=head1 SYNOPSIS

  USAGE:
    batch_game2gff.pl -i InDir -o OutDir

    --indir         # Path to the input directory
    --outdir        # Path to the output directory

=head1 DESCRIPTION

This is what the program does

=head1 COMMAND LINE ARGUMENTS

=head2 Required Arguments

=over 2

=item -i,--infile

Path of the input file.

=item -o,--outfile

Path of the output file.

=back

=head1 Additional Options

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

GNU General Public License, Version 3

L<http://www.gnu.org/licenses/gpl.html>

=head1 AUTHOR

James C. Estill E<lt>JamesEstill at gmail.comE<gt>

=head1 HISTORY

STARTED:

UPDATED:

VERSION: $Rev$

=cut

#-----------------------------------------------------------+
# HISTORY                                                   |
#-----------------------------------------------------------+
#
