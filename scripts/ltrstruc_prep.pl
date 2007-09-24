#!/usr/bin/perl -w
#-----------------------------------------------------------+
#                                                           |
# ltrstruc_prep.pl - Creates files needed for LTR_SRUC      |
#                                                           |
#-----------------------------------------------------------+
#                                                           |
#  AUTHOR: James C. Estill                                  |
# CONTACT: JamesEstill_@_gmail.com                          |
# STARTED: 09/24/2007                                       |
# UPDATED: 09/24/2007                                       |
#                                                           |
# DESCRIPTION:                                              |
#  Given a directory of fasta files with UNIX line endings  |
#  this will create DOS formated files with the txt         |
#  extension. This will also create the flist.txt file.     | 
#                                                           |
# USAGE:                                                    |
#  ltrstruc_prep.pl -i InDir -o OutDir                      |
#                                                           |
# VERSION: $Rev$                                            |
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
my $indir;                      # Input directory, has the fasta files
my $outdir;                     # Output directory, the txt files go here
my $name_root;                  # Root name of the file

# Booleans
my $quiet = 0;
my $verbose = 0;
my $show_help = 0;
my $show_usage = 0;
my $show_man = 0;
my $show_version = 0;

# Counters
my $file_num =0;

#-----------------------------+
# COMMAND LINE OPTIONS        |
#-----------------------------+
my $ok = GetOptions(# REQUIRED OPTIONS
		    "i|indir=s"   => \$indir,
                    "o|outdir=s"  => \$outdir,
		    # ADDITIONAL OPTIONS
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

#-----------------------------------------------------------+
# MAIN PROGRAM BODY                                         |
#-----------------------------------------------------------+

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
# CREATE THE OUT DIR          |
# IF IT DOES NOT EXIST        |
#-----------------------------+
unless (-e $outdir) {
    print "Creating output dir ...\n" if $verbose;
    mkdir $outdir ||
	die "Could not create the output directory:\n$outdir";
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

if ($num_files == 0) {
    print "\a";
    print "No fasta files were found in the input direcotry:\n";
    print "$indir\n";
    print "Fasta file MUST have the fasta or fa extension to be".
	" recognized as fasta files\n";
    exit;
}


my $flist_path = $outdir."flist.txt";
open (FLISTOUT, ">$flist_path")
    || die "Can not open flist for output at:\n$flist_path\n";


for my $ind_file (@fasta_files) {

    $file_num++;   
    
    #-----------------------------+
    # GET ROOT FILE NAME          |
    #-----------------------------+
    if ($ind_file =~ m/(.*)\.fasta$/ ) {	    
	$name_root = "$1";
    }  
    elsif ($ind_file =~ m/(.*)\.fasta$/ ) {	    
	$name_root = "$1";
    } 
    else {
	$name_root = "UNDEFINED";
    }

    # File paths
    my $unix_file_in = $indir.$ind_file;
    my $dos_file_out = $outdir.$name_root.".txt";
    my $dos_file_name = $name_root.".txt";

    print FLISTOUT $dos_file_name."\n";

    unix2dos($unix_file_in, $dos_file_out);


}

close (FLISTOUT);


exit;

#-----------------------------------------------------------+ 
# SUBFUNCTIONS                                              |
#-----------------------------------------------------------+

sub unix2dos {
    
    # Path of the file in UNIX Format
    my $file_in = $_[0];
    my $file_out = $_[1];

    # Original dos line
    my $dos = "\012\015";
    my $unix = "\n";
    # Start lineNo at zero
    my $line_num = 0;
    my $line;

    open( IN, $file_in)  || die ("\n Can not read $file_in \n");
    open( OUT, ">".$file_out) || die ("\n Can not write $file_out \n");

    while( <IN>){
       $line_num++;
       $line = $_;
       # Use the regular expression search operator to replace unix
       $line =~ s/$unix/$dos/g;
       print OUT "$line";
       print "$line_num \n" if $verbose;
    }
    close(IN);
    close(OUT);

    print "Unix2DOX Conversion is complete. \n" if $verbose;

}

sub print_help {

    # Print requested help or exit.
    # Options are to just print the full 
    my ($opt) = @_;

    my $usage = "USAGE:\n". 
	"ltrstruc_prep.pl -i InDir -o OutDir";
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

ltrstruc_prep.pl - Prepare files needed to do a run in LTR_STRUC

=head1 VERSION

This documentation refers to program version 0.1

=head1 SYNOPSIS

  USAGE:
    ltrstruc_prep.pl -i InDir -o OutDir

    --infile        # Path to the input file
    --outfie        # Path to the output file

=head1 DESCRIPTION

This is what the program does

=head1 COMMAND LINE ARGUMENTS

=head 2 Required Arguments

=over 2

=item -i,--infile

Path of the input file.

=item -o,--outfile

Path of the output file.

=back

=head1 Additional Options

=over

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

=head1 HISTORY

STARTED: 09/24/2007

UPDATED: 09/24/2007

VERSION: $Rev$

=cut

#-----------------------------------------------------------+
# HISTORY                                                   |
#-----------------------------------------------------------+
#
