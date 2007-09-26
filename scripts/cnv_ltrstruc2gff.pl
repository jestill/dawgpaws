#!/usr/bin/perl -w
#-----------------------------------------------------------+
#                                                           |
# cnv_ltrstruc2gff.pl - Convert ltrstruc rpt to gff format  |
#                                                           |
#-----------------------------------------------------------+
#                                                           |
#  AUTHOR: James C. Estill                                  |
# CONTACT: JamesEstill_@_gmail.com                          |
# STARTED: 09/25/2007                                       |
# UPDATED: 09/25/2007                                       |
#                                                           |
# DESCRIPTION:                                              |
#  Convert *rpt.txt output files from LTR_STRUC to          |
#  gff formatted annotation coordinates.                    |
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
my $indir;
my $outdir;
my $repdir;
my $flist;                     # Flist is the file list used by LTR_Struc

# Booleans
my $quiet = 0;
my $verbose = 0;
my $show_help = 0;
my $show_usage = 0;
my $show_man = 0;
my $show_version = 0;


#-----------------------------+
# COMMAND LINE OPTIONS        |
#-----------------------------+
my $ok = GetOptions(# REQUIRED OPTIONS
		    "i|indir=s"     => \$indir,
                    "o|outdir=s"    => \$outdir,
		    "r|repdir=s"    => \$repdir,
		    "f|flist=s"     => \$flist,
		    # ADDITIONAL OPTIONS
		    "q|quiet"       => \$quiet,
		    "verbose"       => \$verbose,
		    # ADDITIONAL INFORMATION
		    "usage"         => \$show_usage,
		    "version"       => \$show_version,
		    "man"           => \$show_man,
		    "h|help"        => \$show_help,);

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

# THROW ERROR IF REQRUIED VARIABLES NOT SPECIFIED
if ( (!$indir) || (!$outdir) || (!$repdir) || (!$flist) ) {
    print "ERROR: Input directory must be specified\n" if !$indir;
    print "ERROR: Output directory must be specified\n" if !$outdir;
    print "ERROR: Reports directory must be specified\n" if !$repdir;
    print "\n";
    print_help("full");
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

unless ($repdir =~ /\/$/ ) {
    $repdir = $repdir."/";
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
opendir( DIR, $indir ) 
    || die "Can't open directory:\n$indir"; 
my @fasta_files = grep /\.fasta$|\.fa$/, readdir DIR ;
closedir( DIR );
my $num_fasta_files = @fasta_files;

if ($num_fasta_files == 0) {
    print "\a";
    print "No fasta files were found in the input direcotry:\n";
    print "$indir\n";
    print "Fasta file MUST have the fasta or fa extension to be".
	" recognized as fasta files\n";
    exit;
}


#-----------------------------+
# Get the LTR_STRUC report    |
# files                       |
#-----------------------------+
opendir(REPDIR,$repdir) 
    || die "Can not open results direcotry:\n$repdir\n";
my @report_files = grep /rprt\.txt$/, readdir REPDIR;
@report_files = sort(@report_files);
closedir (REPDIR); 
# Sort the array of reports

my $num_report_files = @report_files;

print "\n-----------------------------------------------\n";
print "Report files to process: $num_report_files\n";
print "-----------------------------------------------\n\n";

my $report_num = 0;
for my $ind_report (@report_files) {
    
    $report_num++;

    print "Processing report $report_num of $num_report_files: $ind_report\n";
    
#    my $fasta_file = $indir.

#    ltrstruc2gff ()
    
    if ($report_num = 1) {exit;}
    
}


# Test of using the ltrstruc2gff subfunction on a single sequence


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


sub ltrstruc2gff {
    
    # VARS PASSED TO THE SUBFUNCTION
    my ($fasta_in, $report_in, $gff_out, $gff_append) = @_;

    # FASTA RELATED VARS
    my $qry_seq;

    # LTR STRUC VARS
    my $ls_score;
    my $ls_contig_len;
    my $ls_orientation;
    my $ls_longest_orf;
    my $ls_rt_frame1;
    my $ls_rt_frame2;
    my $ls_rt_frame3;
    my $ls_5ltr_len;
    my $ls_3ltr_len;
    my $ls_ltr_homology;
    my $ls_5dinuc;
    my $ls_3dinuc;
    my $ls_5flank_seq;
    my $ls_3flank_seq;
    my $ls_pbs_seq;
    my $ls_ppt_seq;
    my $ls_5id_seq;
    my $ls_3id_seq;
    my $ls_5ltr_seq;
    my $ls_3ltr_seq;
    my $ls_full_retro_seq;

    # BOOLEANS
    my $in_rt_frame_1 = 0;
    my $in_rt_frame_2 = 0;
    my $in_rt_frame_3 = 0;
    my $in_5id_seq = 0;
    my $in_3id_seq = 0;
    my $in_5ltr = 0;
    my $in_3ltr = 0;
    my $in_complete_seq = 0;

    #-----------------------------+
    # OPEN GFF OUTPUT FILE        |
    #-----------------------------+
    if ($gff_append) {
	open (GFFOUT, ">>gff_out") 
	    || die "ERROR: Could not open gff output file:\b$gff_out\n"
    }
    else {
	open (GFFOUT, ">$gff_out") 
	    || die "ERROR: Could not open gff output file:\n$gff_out\n";
    }
    
    #-----------------------------+
    # LOAD FASTA SEQUENCE         |
    #-----------------------------+
    open (INFASTA, "$fasta_in") 
	|| die "Can not open fasta input file:\n$fasta_in\n";
    while (<INFASTA>) {
	chomp;
	unless(m/^\>/) {
	    $qry_seq = $qry_seq.$_;
	}
    }
    close (INFASTA);
    
    #-----------------------------+
    # GET DATA FROM REPORT FILE   |
    #-----------------------------+
    open (REPIN, "<$report_in")
	|| die "ERROR: Can not open report file:\n$report_in\n";

    while (<REPIN>) {
	chomp;
	if ($in_rt_frame_1) {

	}
	elsif ($in_rt_frame_2) {

	}
	elsif ($in_rt_frame_3) {

	}
	elsif ($in_5id_seq) {

	}
	elsif ($in_3id_seq) {

	}
	elsif ($in_5ltr) {
	    
	}
	elsif ($in_3ltr) {
	    # Jump out on the first empty line encountered
	}
	elsif ($in_complete_seq) {
	    $ls_full_retro_seq = $ls_full_retro_seq.$_;
	}
    }
    close (REPIN);

    close (GFFOUT);

}

=head1 NAME

cnv_ltrstruc2gff.pl - Convert LTR_STRUC report output files to gff

=head1 VERSION

This documentation refers to program version $Rev$

=head1 SYNOPSIS
    
  USAGE:
    cnv_ltrstruc2gff.pl -i InDir -o OutDir -f FastaDir
    
    --indir         # Directory with the fasta files
    --outdir        # Directory for the base output dir
    --results       # Directory containing the LTR_STRUC results

=head1 DESCRIPTION

Given a directory containing the output from LTR_STRUC and a directory
containing the fasta files that the structural predictions were made
from,

=head1 COMMAND LINE ARGUMENTS

=head 2 Required Arguments

=over 2

=item -i,--indir

Path of the intput directory containing the fasta files that were
analyzed by LTR_STRUC.

=item -o,--outdir

Path of the output directory that will serve as the base for the
output from the conversion to gff.

=item -r, --results

Path of the directory containing the results from LTR_STRUC. It is 
expected that these file names will end with rprt.txt.

=back

=head2 Additional Options

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

STARTED:

UPDATED:

VERSION: $Rev$

=cut

#-----------------------------------------------------------+
# HISTORY                                                   |
#-----------------------------------------------------------+
#
