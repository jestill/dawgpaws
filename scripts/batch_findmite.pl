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

package DAWGPAWS;

#-----------------------------+
# INCLUDES                    |
#-----------------------------+
use strict;
use Getopt::Long;

# Array to hold the findmite parameters
my @fm_params = ();

#-----------------------------+
# PROGRAM VARIABLES           |
#-----------------------------+
my ($VERSION) = q$Rev$ =~ /(\d+)/;

#-----------------------------+
# VARIABLE SCOPE              |
#-----------------------------+
my $indir;                     # Input directory
my $outdir;                    # Output directory
my $parfile;                   # Parameters file
my $name_root;                 # Root name of the input file

# Counters/Index Vals
my $i = 0;                     # Array index val
my $file_num = 0;              # File number
my $proc_num = 0;              # Process number

# Booleans
my $verbose = 0;
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
		    "i|indir=s"       => \$indir,
                    "o|outdir=s"      => \$outdir,
		    "p|params=s"      => \$parfile,
		    # Additional options
		    "q|quiet"         => \$quiet,
		    "verbose"         => \$verbose,
		    # Additional information
		    "usage"           => \$show_usage,
		    "version"         => \$show_version,
		    "man"             => \$show_man,
		    "h|help"          => \$show_help,);

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
# LOAD FINDMITE PARAMATERS    |
#-----------------------------+
open (PARFILE, "<$parfile")
    || die "Can not open the paramter file:\n$parfile\n";

while (<PARFILE>) {
    chomp;

#    # Use the following if an error crops up later
#    print "IN: $_\n";
    
    # Ignores comment lines starting with #
    unless (m/^\#/) {
	#my @in_line = split(/\t/, $_);
	my @in_line = split;           # Implicit split of $_ by whitespace
	my $num_in_line = @in_line; 
	$fm_params[$i][0] = $in_line[0] || "NULL"; 
	$fm_params[$i][1] = $in_line[1] || "NULL";
	$fm_params[$i][2] = $in_line[2] || "NULL";
	$fm_params[$i][3] = $in_line[3] || "NULL";
	$fm_params[$i][4] = $in_line[4] || "NULL";
	$fm_params[$i][5] = $in_line[5] || "NULL";
	$fm_params[$i][6] = $in_line[6] || "NULL";
	$fm_params[$i][7] = $in_line[7] || "NULL";
	$fm_params[$i][8] = $in_line[8] || "NULL";
	$fm_params[$i][9] = $in_line[9] || "NULL";
	
#        # Use the following if a error crops up later
#	if ($verbose) {
#	    print STDERR "\tProcessed as:\n";
#	    print STDERR "\t".$fm_params[$i][0]."\n";
#	    print STDERR "\t".$fm_params[$i][1]."\n";
#	    print STDERR "\t".$fm_params[$i][2]."\n";
#	    print STDERR "\t".$fm_params[$i][3]."\n";
#	    print STDERR "\t".$fm_params[$i][4]."\n";
#	    print STDERR "\t".$fm_params[$i][5]."\n";
#	    print STDERR "\t".$fm_params[$i][6]."\n";
#	    print STDERR "\t".$fm_params[$i][7]."\n";
#	    print STDERR "\t".$fm_params[$i][8]."\n";
#	    print STDERR "\t".$fm_params[$i][9]."\n";
#	}
	
	$i++;
    }

} # End of while INFILE

my $num_par_sets = $i;
my $max_i = $i-1;

close PARFILE;


# TO DO: CHECK PARAM FILE VALIDITY


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

my $num_proc_total = $num_files * $num_par_sets;

#-----------------------------------------------------------+
#-----------------------------------------------------------+
#-----------------------------------------------------------+
#-----------------------------------------------------------+

# FINDMITE REQUIRES INPUT FOR THE FOLLOWING QUESTIONS
# Maximal fragment length : ie. 14000
# Direct repeast : ie. TA
#     From Mobile DNA II:
#               TA 
#               TWA = TTA or TAA  
#               TWW = TTA or TAA or TTA or TTT
#               WWW = TTA or TAA or TTA or TTT or 
#                     ATA or AAA or ATA or ATT
#               TNA = TAT or TGT or TCT or TTT
#
# Length of TIR (3-15): ie 11
# Number of mismatches(0-3) : ie 1
# Filtering A/T strings (y/n) : y                 -- Set at cmd line
# Filtering C/G strings (y/n) : y                 -- set at cmd line
# Filtering AT/TA repeats (y/n) : y               -- Set at cmd line
# Filtering 2 Base (0 -- don't filter) : 85   [Range 0 to 100]
# Minimum distance : 30
# Maximum distance : 700

print STDERR "$num_proc_total findmite runs to process\n";

for my $ind_file (@fasta_files) {

    $file_num++;
    
    # Get root file name
    if ($ind_file =~ m/(.*)\.fasta$/ ) {	    
	$name_root = "$1";
    }  
    elsif ($ind_file =~ m/(.*)\.fa$/ ) {	    
	$name_root = "$1";
    } 
    else {
	$name_root = "UNDEFINED";
    }

    #-----------------------------+
    # CREATE ROOT NAME DIR        |
    #-----------------------------+
    my $name_root_dir = $outdir.$name_root."/";
    unless (-e $name_root_dir) {
	mkdir $name_root_dir ||
	    die "Could not create dir:\n$name_root_dir\n"
    }

    #-----------------------------+
    # CREATE GENSCAN OUTDIR       |
    #-----------------------------+
    # Dir to hold gene prediction output from local software
    my $findmite_dir = $name_root_dir."findmite/";
    unless (-e $findmite_dir) {
	mkdir $findmite_dir ||
	    die "Could not create genscan out dir:\n$findmite_dir\n";
    }

    #-----------------------------------------------------------+
    # RUN FINDMITE FOR EACH SET OF PARAM VALUES                 |
    #-----------------------------------------------------------+
    for ($i=0; $i<$num_par_sets; $i++) {
	$proc_num++;
	
	# Load array vals to usefull short names
	my $mite_name    = $fm_params[$i][0];
	my $rep          = $fm_params[$i][1];
	my $tir_len      = $fm_params[$i][2];
	my $num_mismatch = $fm_params[$i][3];
	my $filt_at      = $fm_params[$i][4];
	my $filt_cg      = $fm_params[$i][5];
	my $filt_atta    = $fm_params[$i][6];
	my $filt_2_base  = $fm_params[$i][7];
	my $min_dist     = $fm_params[$i][8];
	my $max_dist     = $fm_params[$i][9];

	if ($verbose) {
	    print STDERR "\n\n========================================\n";
	    print STDERR " FindMITE Process $proc_num of $num_proc_total\n";
	    print STDERR " ".$name_root."_".$mite_name."\n";
	    print STDERR "========================================\n";
	}


	my $fm_out = $findmite_dir.$name_root."_".$mite_name.".mite.txt";
	my $fm_cmd = "FINDMITE1New.bin $ind_file $fm_out";
	print "CMD: $fm_cmd\n" if $verbose;

	# OPEN THE FINDMITE PROGRAM 
	# Then should be made a subfunction ..
	open(FINDMITE,"|$fm_cmd");
	print FINDMITE "150000\n";
	print FINDMITE "$rep\n";
	print FINDMITE "$tir_len\n";
	print FINDMITE "$num_mismatch\n";
	print FINDMITE "$filt_at\n";
	print FINDMITE "$filt_cg\n";
	print FINDMITE "$filt_atta\n";
	print FINDMITE "$filt_2_base\n";
	print FINDMITE "$min_dist\n";
	print FINDMITE "$max_dist\n";
	close FINDMITE;

	# print return to clean up the command line
	print "\n";

    } # End of run find mite for each set of param vals
    
    

} # End of for each $ind_file


exit;

#-----------------------------------------------------------+ 
# SUBFUNCTIONS                                              |
#-----------------------------------------------------------+

sub sub_run_findmite {
    # Copied from sub_run in the bioperl 

}

sub print_help {

    # Print requested help or exit.
    # Options are to just print the full 
    my ($opt) = @_;

    my $usage = "USAGE:\n". 
	"batch_findmite.pl -i InDir -o OutDir";
    my $args = "REQUIRED ARGUMENTS:\n".
	"  --indir        # Directory of fasta files to process\n".
	"  --outdir       # Path to the output direcotry\n".
	"  --params       # Path to the params file\n".
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
