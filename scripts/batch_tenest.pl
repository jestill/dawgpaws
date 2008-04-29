#!/usr/bin/perl -w
#-----------------------------------------------------------+
#                                                           |
# batch_tenest.pl - Run TeNest in batch mode                |
#                                                           |
#-----------------------------------------------------------+
#                                                           |
#  AUTHOR: James C. Estill                                  |
# CONTACT: JamesEstill_@_gmail.com                          |
# STARTED: 04/29/2008                                       |
# UPDATED: 04/29/2008                                       |
#                                                           |
# DESCRIPTION:                                              |
#  Run TeNest in batch mode.
#                                                           |
# USAGE:                                                    |
#  ShortFasta Infile.fasta Outfile.fasta                    |
#                                                           |
# VERSION: $Rev$                                            |
#                                                           |
#                                                           |
# LICENSE:                                                  |
#  GNU General Public License, Version 3                    |
#  http://www.gnu.org/licenses/gpl.html                     |  
#                                                           |
#-----------------------------------------------------------+
# NOTE: TENEST Assumes that directories do not end with /

package DAWGPAWS;

#-----------------------------+
# INCLUDES                    |
#-----------------------------+
use strict;
use Getopt::Long;
# The following needed for printing help
use Pod::Select;               # Print subsections of POD documentation
use Pod::Text;                 # Print POD doc as formatted text file
use IO::Scalar;                # For print_help subfunction
use IO::Pipe;                  # Pipe for STDIN, STDOUT for POD docs
use File::Spec;                # Convert a relative path to an abosolute path

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

# TE NEST VARIABLES WITH DEFAULT VALUES
my $organism_db = "maize";
my $wublast_dir = "/usr/local/genome/wu_blast/";
my $home_dir = $ENV{'HOME'};
my $te_db_dir = "$home_dir/apps/te_nest";

#-----------------------------+
# COMMAND LINE OPTIONS        |
#-----------------------------+
my $ok = GetOptions(# REQUIRED OPTIONS
		    "i|indir=s"    => \$indir,
                    "o|outdir=s"   => \$outdir,
		    # ADDITIONAL OPTIONS
		    "org=s"        => \$organism_db,
		    "blast-dir=s"  => \$wublast_dir,
		    "tenest-dir=s" => \$te_db_dir,
		    "q|quiet"      => \$quiet,
		    "verbose"      => \$verbose,
		    # ADDITIONAL INFORMATION
		    "usage"        => \$show_usage,
		    "version"      => \$show_version,
		    "man"          => \$show_man,
		    "h|help"       => \$show_help,);

#-----------------------------+
# SHOW REQUESTED HELP         |
#-----------------------------+
if ( ($show_usage) ) {
#    print_help ("usage", File::Spec->rel2abs($0) );
    print_help ("usage", $0 );
}

if ( ($show_help) || (!$ok) ) {
#    print_help ("help",  File::Spec->rel2abs($0) );
    print_help ("help",  $0 );
}

if ($show_man) {
    # User perldoc to generate the man documentation.
    system ("perldoc $0");
    exit($ok ? 0 : 2);
}

if ($show_version) {
    print "\batch_genmark.pl:\n".
        "Version: $VERSION\n\n";
    exit;
}


#-----------------------------+
# CHECK REQUIRED ARGS         |
#-----------------------------+
if ( (!$indir) || (!$outdir) ) {
    print "\a";
    print STDERR "\n";
    print STDERR "ERROR: An input directory was not specified at the".
        " command line\n" if (!$indir);
    print STDERR "ERROR: An output directory was specified at the".
        " command line\n" if (!$outdir);
    print_help("full");
}



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

#unless ($te_db_dir =~ /\/$/ ) {
#    $te_db_dir = $te_db_dir."/";
#}

#-----------------------------+
# CHECK FOR REQUIRED PROGRAMS |
#-----------------------------+
unless (-e $te_db_dir."/TEnest.pl") {
    die "TENest.pl does not exist at:\n$te_db_dir\n"
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
my $num_files = @fasta_files;

#-----------------------------+
# SHOW ERROR IF NO FILES      |
# WERE FOUND IN THE INPUT DIR |
#-----------------------------+
if ($num_files == 0) {
    print STDERR "\a";
    print STDERR "\nERROR: No fasta files were found in the input directory\n".
        "$indir\n".
        "Fasta files must have the fasta or fa extension.\n\n";
    exit;
}

#-----------------------------+
# CREATE THE OUT DIR          |
# IF IT DOES NOT EXIST        |
#-----------------------------+
unless (-e $outdir) {
    print STDERR "Creating output dir ...\n" unless $quiet;
    mkdir $outdir ||
        die "Could not create the output directory:\n$outdir";
}


#-----------------------------------------------------------+ 
# MAIN PROGRAM BODY                                         |
#-----------------------------------------------------------+ 
my $proc_num = 0;
my $file_num =0;

for my $ind_file (@fasta_files)
{

    my $name_root;

    $proc_num++;
    $file_num++;

    #-----------------------------+
    # GET THE ROOT NAME OF THE    |
    # FASTA FILE                  |
    #-----------------------------+
    if ($ind_file =~ m/(.*)\.hard\.fasta$/) {
        # file ends in .hard.fasta
        # This is hard masked fasta files
        $name_root = "$1";
    }
    elsif ($ind_file =~ m/(.*)\.masked\.fasta$/) {
        # file ends in .masked.fasta
        $name_root = "$1";
    }
    elsif ($ind_file =~ m/(.*)\.fasta$/ ) {
        # file ends in .fasta
        $name_root = "$1";
    }
    elsif ($ind_file =~ m/(.*)\.fa$/ ) {
        # file ends in .fa
        $name_root = "$1";
    }
    else {
        $name_root = $ind_file;
    }


    #-----------------------------+
    # CREATE ROOT NAME DIR        |
    #-----------------------------+
    my $name_root_dir = $outdir.$name_root;
    unless (-e $name_root_dir) {
        mkdir $name_root_dir ||
            die "Could not create dir:\n$name_root_dir\n"
	}

    #-----------------------------+
    # CREATE TENEST OUTDIR        |
    #-----------------------------+
    # Dir to hold gene prediction output from local software
    my $tenest_out_dir = $outdir.$name_root."/tenest";
    unless (-e $tenest_out_dir) {
        mkdir $tenest_out_dir ||
            die "Could not create tenest out dir:\n$tenest_out_dir\n";
    }

    #-----------------------------+
    # CREATE GFF OUTDIR           |
    #-----------------------------+
    # Dir to hold gene prediction output from local software
    my $gff_out_dir = $outdir.$name_root."/gff/";
    #print STDERR "$gff_out_dir\n";
    unless (-e $gff_out_dir) {
        mkdir $gff_out_dir ||
            die "Could not create gff out dir:\n$gff_out_dir\n";
    }

    #-----------------------------+
    # RUN TE NEST                 |
    #-----------------------------+

    # TE NEST COMMAND
    my $te_nest_cmd = $te_db_dir."/TEnest.pl ".$indir.$ind_file.
	" --output ".$tenest_out_dir;

# Trying usage as program OPTIONS INFILE
#    my $te_nest_cmd = $te_db_dir."/TEnest.pl ".
#	" --output ".$tenest_out_dir.
#	$indir.$ind_file;
#	" --org ".$organism_db;
#	" --blast ".$wublast_dir.
#	" --current ".$te_db_dir;

    print STDERR "CMD: $te_nest_cmd\n" if $verbose;
    print STDERR "\n" if $verbose;
	
    # RUN THE TE NEST COMMAND
    system ($te_nest_cmd);
    
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


=head1 NAME

batch_tenest.pl - Run TENest in batch mode

=head1 VERSION

This documentation refers to program version $Rev$

=head1 SYNOPSIS

  USAGE:
    batch_tenest.pl -i indir -o outdir

    -i,--indir         # Dir containing fasta files to process 
    -o,--outfie        # Dir to place the output in

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

STARTED: 09/29/2008

UPDATED: 09/29/2008

VERSION: $Rev$

=cut

#-----------------------------------------------------------+
# HISTORY                                                   |
#-----------------------------------------------------------+
#
# 09/29/2008
# - Program started
# - Basic idea is to run the tenest program from 
