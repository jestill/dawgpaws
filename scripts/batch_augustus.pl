#!/usr/bin/perl -w
#-----------------------------------------------------------+
#                                                           |
# batch_augustus.pl - Run Augustus gene prediction in batch |
#                                                           |
#-----------------------------------------------------------+
#                                                           |
#  AUTHOR: James C. Estill                                  |
# CONTACT: JamesEstill_@_gmail.com                          |
# STARTED: 06/01/2011                                       |
# UPDATED: 07/20/2011                                       |
#                                                           |
# DESCRIPTION:                                              |
#  Run Augustus gene prediction program in batch mode for   |
#  a directory of contigs. The configuration file uses      |
#  param name to tag sets of parameters to allow for        |
#  multiple parameter sets for each contig.                 |
#                                                           |
# USAGE:                                                    |
#  batch_augustus.pl -i indir -o outdir -c config.cfg       |
#                                                           |
# VERSION: $Rev$                                      |
#                                                           |
# LICENSE:                                                  |
#  GNU General Public License, Version 3                    |
#  http://www.gnu.org/licenses/gpl.html                     |  
#                                                           |
#-----------------------------------------------------------+
# TO DO: Deal with Augustus hanging on large contigs ... may need
#        to have timer and check that last line of outfile 
#        contains the command before killing the job.


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
my $indir;                     # Files to process
my $outdir;                    # Base output dir
my $config_file;               # Configuration file
my @aug_params = ();           # Augustus Parameters
                               # One row for each parameter set
my $gff_ver = "GFF3";

# PATH TO AUGUSTUS BINARY
# Path to the Augustus binary file otherwise
# assumes the program already located in users path
my $aug_path = $ENV{AUGUSTUS_BIN_PATH} ||
    "augustus";                  
my $program = "AUGUSTUS";
my $param;

# BOOLEANS
my $quiet = 0;
my $verbose = 0;
my $show_help = 0;
my $show_usage = 0;
my $show_man = 0;
my $show_version = 0;
my $do_test = 0;                  # Run the program in test mode
my $apollo_format = 0;            # Apollo can't render CDS so use exon
my $append_gff = 0;               # Append results to an existing GFF 
my $delim = ":";                  # Character delimiter for source/suffix

#-----------------------------+
# COMMAND LINE OPTIONS        |
#-----------------------------+
my $ok = GetOptions(# REQUIRED OPTIONS
		    "i|indir=s"      => \$indir,
                    "o|outdir=s"     => \$outdir,
                    "c|config=s"     => \$config_file,
		    # ADDITIONAL OPTIONS
		    "p|program=s"    => \$program,
		    "augustus-bin=s" => \$aug_path,
		    "gff-ver=s"      => \$gff_ver,
		    "d|delim=s"      => \$delim,
		    "apollo"         => \$apollo_format,
		    "q|quiet"        => \$quiet,
		    "verbose"        => \$verbose,
		    # ADDITIONAL INFORMATION
		    "usage"          => \$show_usage,
		    "test"           => \$do_test,
		    "version"        => \$show_version,
		    "man"            => \$show_man,
		    "h|help"         => \$show_help,);

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
    print "\nbatch_augustus.pl:\n".
	"Version: $VERSION\n\n";
    exit;
}

#-----------------------------+
# STANDARDIZE GFF VERSION     |
#-----------------------------+
unless ($gff_ver =~ "GFF3" || 
	$gff_ver =~ "GFF2") {
    # Attempt to standardize GFF format names
    if ($gff_ver =~ "3") {
	$gff_ver = "GFF3";
    }
    elsif ($gff_ver =~ "2") {
	$gff_ver = "GFF2";
    }
    else {
	print "\a";
	die "The gff-version \'$gff_ver\' is not recognized\n".
	    "The options GFF2 or GFF3 are supported\n";
    }
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
    print_help ("usage", $0 );
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

#-----------------------------+
# Get the FASTA files from the|
# directory provided by the   |
# var $indir                  |
#-----------------------------+
opendir( DIR, $indir ) || 
    die "Can't open directory:\n$indir"; 
my @fasta_files = grep /\.fasta$|\.fa$/, readdir DIR ;
closedir( DIR );

my $count_files = @fasta_files;

#-----------------------------+
# SHOW ERROR IF NO FILES      |
# WERE FOUND IN THE INPUT DIR |
#-----------------------------+
if ($count_files == 0) {
    print STDERR "\a";
    print STDERR "\nERROR: No fasta files were found in the input directory\n".
	"$indir\n".
	"Fasta files must have the fasta or fa extension.\n\n";
    exit;
}

print STDERR "NUMBER OF FILES TO PROCESS: $count_files\n" if $verbose;


#-----------------------------+
# PARSE THE CONFIG FILE       |
#-----------------------------+
my $i=0;
my $config_line_num=0;

open (CONFIG, "<$config_file") ||
    die "ERROR Can not open the config file:\n $config_file";

while (<CONFIG>) {
    chomp;
    $config_line_num++;
    unless (m/^\#/) {
       	my @in_line = split (/\t/);           # Implicit split of $_ by tab
	my $num_in_line = @in_line; 
	
	# Config file must be atleast two
	if ($num_in_line == 2 ||
	    $num_in_line == 3) { 
	    $aug_params[$i][0] = $in_line[0];            # Name
	    $aug_params[$i][1] = $in_line[1];            # Species
	    $aug_params[$i][2] = $in_line[2] || "NULL";  # Options
	    $i++;
	} # End of if $num_in_line is more than three or less than two
	else {
	    print "\a";
	    print STDERR "WARNING: Unexpected number of line in config".
		" file line $config_line_num\n$config_file\n";
	}

   } # End of unless comment line
} # End of while CONFIG file
close CONFIG;

# Number of parameter sets specified in the config file
my $num_par_sets = $i;

if ($num_par_sets == 0) {
    print "\a";
    print STDERR "ERROR: No parameter sets were found in the config file:\n".
	"$config_file\n";
    exit;
}

my $num_proc_total = $count_files * $num_par_sets;

print STDERR "$num_proc_total AUGUSTUS runs to process\n\n" if $verbose;

#-----------------------------------------------------------+
# MAIN PROGRAM BODY                                         |
#-----------------------------------------------------------+
my $file_num = 0;
my $proc_num = 0;

for my $ind_file (@fasta_files) {

    $file_num++;
    my $name_root;
    # Get root file name
    if ($ind_file =~ m/(.*)\.masked\.fasta$/) {
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
    my $name_root_dir = $outdir.$name_root."/";
    unless (-e $name_root_dir) {
	mkdir $name_root_dir ||
	    die "Could not create dir:\n$name_root_dir\n"
    }

    #-----------------------------+
    # CREATE AUGUSTUS OUTDIR      |
    #-----------------------------+
    # Dir to hold gene prediction output from local software
    my $aug_dir = $name_root_dir."augustus/";
    unless (-e $aug_dir) {
	mkdir $aug_dir ||
	    die "Could not create augustus output dir:\n$aug_dir\n";
    }
    
    #-----------------------------+
    # CREATE GFF OUTDIR           |
    #-----------------------------+
    my $gff_dir = $name_root_dir."gff/";
    unless (-e $gff_dir) {
	mkdir $gff_dir ||
	    die "Could not create gff out dir:\n$gff_dir\n";
    }

    #-----------------------------+
    # FOR EACH PARAM SET  IN THE  |
    # CONFIG FILE                 |
    #-----------------------------+
    my $gff_count=0;               # Count of the number of gff files

    for ($i=0; $i<$num_par_sets; $i++) {

	$proc_num++;
	my $aug_param_name = $aug_params[$i][0];
	my $aug_species = $aug_params[$i][1];
	my $aug_opts = chomp($aug_params[$i][2]) 
	    || "NULL";

	#-----------------------------+
	# SET THE AUGUSUTS COMMAND    |
	#-----------------------------+
	my $aug_res = $aug_dir.$name_root."_augustus_"
	    .$aug_param_name.".gff";

	# Be default we will ask for output in GFF3 format,
	# This will be normalized to the standard GFF3 format
	# used for the rest of DAWGPAWS, and will be simplified
	# to GFF2 format.
	my $aug_cmd = $aug_path.
	    " --species="."$aug_species".
	    " --gff3=on".
	    " ".$indir.$ind_file.
	    " --outfile=".$aug_res;
	
	unless ($aug_opts ="NULL") {
	    $aug_cmd = $aug_cmd." ".$aug_opts;
	}

	# PRINT the command
	print STDERR "CMD: ".$aug_cmd."\n"
	    if $verbose;

	# RUN Augustus
	system ($aug_cmd) unless $do_test;

	# CONVERT THE AUGUSTUS OUTPUT TO GFF
	if (-e $aug_res) {

	    my $aug_gff = $gff_dir.$name_root."_augustus_"
		.$aug_param_name.".gff";

	    augustus2gff ( $program,
			   $aug_res,
			   $aug_gff,
			   $name_root,
			   $aug_param_name,
			   $delim,
			   $append_gff,
			   $apollo_format
			  )

	    # vars to send are source
	    # $aug_in in file path
	    # $gffout gff out file path
	    # seq ID of the sequence being analyszed seq_id
	    # src_suffix (this is the parameter name) 
	    # do append -- if we are appending to an existing gff file

	}

    } # End of for each fasta file for each parameter set


}


exit 0;

#-----------------------------------------------------------+ 
# SUBFUNCTIONS                                              |
#-----------------------------------------------------------+


sub augustus2gff {
# Technically Augustus returns data in GFF format, this modifies the 
# Augustus native output to the GFF format that is standard for 
# DAWGPAWS

    
    my ($source, $aug_in, $gffout, $seq_id, $src_suffix, $d,
	$do_append, $do_apollo_format) = @_;
    # Array to hold all snap results for a single contig

#    my @aug_results;

    my $i = -1;
    my $j = -1;
    my $model_num = 0;

    # The following delimiter (:) works for Apollo, but throws problems
    # with GBROWSE
    if ($src_suffix) {
	$source = $source.":".$src_suffix;
    }

    my $attribute;

    #-----------------------------+
    # OPEN INPUT FILE HANDLE      |
    #-----------------------------+
    if ($aug_in) {
	open (INFILE, "<$aug_in") ||
	    die "ERROR: Can not open Augustus result file\n $aug_in\n";
	
    }
    else {
	print STDERR "Expecting input from STDIN\n";
	open (INFILE, "<&STDIN") ||
	    die "Can not accept input from standard input.\n";
    }
    
    #-----------------------------+
    # OPEN GFF OUTFILE            |
    #-----------------------------+
    # Default to STDOUT if no argument given
    if ($gffout) {
	if ($do_append) {
	    open (GFFOUT, ">>$gffout") ||
		die "ERROR: Can not open gff outfile:\n $gffout\n";
	}
	else {
	    open (GFFOUT,">$gffout") ||
		die "ERROR: Can not open gff outfile:\n $gffout\n";
	}
    } 
    else {
	open (GFFOUT, ">&STDOUT") ||
	    die "Can not print to STDOUT\n";
    }

    if ($gff_ver =~ "GFF3") {
	print GFFOUT "##gff-version 3\n";
    }


    #-----------------------------+
    # PROCESS AUGUSTUS OUTPUT
    #-----------------------------+
    my $prev_gene_name = "NULL";

    while (<INFILE>) {
	chomp;
	# Skip comment lines
	next if m/^\#/;

	print STDERR $_."\n"
	    if $verbose;

	my @gff_parts = split;
	my $num_gff_parts = @gff_parts;

	#-----------------------------+
	# MAKE START < END            |
	#-----------------------------+
	my $start;
	my $end;
	if ($gff_parts[3] < $gff_parts[4]) {
	    $start = $gff_parts[3];
	    $end = $gff_parts[4];
	} else {
	    $end = $gff_parts[4];
	    $start = $gff_parts[3];
	}

	my $type = $gff_parts[2];
	my $attribute = $gff_parts[8];

	# WE MAY NEED TO ORDER PARTS OF A SINGLE GENE
	# SEQUENTIALY BY OCCURENCE ON THE CONTIG

	# It may also be necessary to include contig name in result
	# name to that results for multiple contigs can have
	# unique IDs.
	# Load results to the $aug_results hash
#	$i++;
#	$aug_results[$i][0] = $gff_parts[0];
#	$aug_results[$i][1] = $program.":".$src_suffix;
#	$aug_results[$i][2] = $gff_parts[2];
#	$aug_results[$i][3] = $start;
#	$aug_results[$i][4] = $end;
#	$aug_results[$i][5] = $gff_parts[5];
#	$aug_results[$i][6] = $gff_parts[6];
#	$aug_results[$i][7] = $gff_parts[7];
#	$aug_results[$i][8] = $gff_parts[8];


	#/////////////////
	# Kludge alert
	#/////////////////
	# For GFF2 will only report the CDS
	if ($gff_ver =~ "GFF2") {
	    next unless $type =~ "CDS";

	    # Replace longer GFF3 type name
	    if ($attribute =~ m/Parent\=(.*)/) {
		print STDERR "YUP it is ".$1."\n";
		$attribute = $1;
	    }
	    
	    # Apollo can't render CDS features so 
	    # converting type to exon for Apollo
	    if ($do_apollo_format) {
		$type = "exon";
	    }
	}
	else {
	    print STDERR "NOPE";
	}

	print GFFOUT $gff_parts[0]."\t".
	    $program.$d.$src_suffix."\t".
	    $type."\t".
	    $start."\t".
	    $end."\t".
	    $gff_parts[5]."\t".
	    $gff_parts[6]."\t".
	    $gff_parts[7]."\t".
	    $attribute."\n";

    } # End of while INFILE

    close INFILE;
    close GFFOUT;

}

sub print_help {
    my ($help_msg, $podfile) =  @_;
    # help_msg is the type of help msg to use (ie. help vs. usage)
    
    print "\n";
    
    #-----------------------------+
    # PIPE WITHIN PERL            |
    #-----------------------------+
    # This code made possible by:
    # http://www.perlmonks.org/index.pl?node_id=76409
    # Tie info developed on:
    # http://www.perlmonks.org/index.pl?node=perltie 
    #
    #my $podfile = $0;
    my $scalar = '';
    tie *STDOUT, 'IO::Scalar', \$scalar;
    
    if ($help_msg =~ "usage") {
	podselect({-sections => ["SYNOPSIS|MORE"]}, $0);
    }
    else {
	podselect({-sections => ["SYNOPSIS|ARGUMENTS|OPTIONS|MORE"]}, $0);
    }

    untie *STDOUT;
    # now $scalar contains the pod from $podfile you can see this below
    #print $scalar;

    my $pipe = IO::Pipe->new()
	or die "failed to create pipe: $!";
    
    my ($pid,$fd);

    if ( $pid = fork() ) { #parent
	open(TMPSTDIN, "<&STDIN")
	    or die "failed to dup stdin to tmp: $!";
	$pipe->reader();
	$fd = $pipe->fileno;
	open(STDIN, "<&=$fd")
	    or die "failed to dup \$fd to STDIN: $!";
	my $pod_txt = Pod::Text->new (sentence => 0, width => 78);
	$pod_txt->parse_from_filehandle;
	# END AT WORK HERE
	open(STDIN, "<&TMPSTDIN")
	    or die "failed to restore dup'ed stdin: $!";
    }
    else { #child
	$pipe->writer();
	$pipe->print($scalar);
	$pipe->close();	
	exit 0;
    }
    
    $pipe->close();
    close TMPSTDIN;

    print "\n";

    exit 0;
   
}

1;
__END__

=head1 NAME

batch_augustus.pl - Run Augustus gene prediction in batches

=head1 VERSION

This documentation refers to program version $Rev$

=head1 SYNOPSIS

  USAGE:
    batch_augustus.pl -i indor -o outdir -c config

    -i,--infile        # Path to the input directory
    -o,--outfie        # Path to the output directory
    -c,--config        # Configuraiton file

=head1 DESCRIPTION

Run the Augustus gene prediction program in batch mode for
a directory of sequence contigs. The configuration file uses the param
name to tag sets of parameters to allow for multiple
parameter sets for each contig.

=head1 REQUIRED ARGUMENTS

=over 2

=item -i,--indir

Path of the directory containing the sequences to process.

=item -o,--outdir

Path of the directory to place the program output.

=item -c, --config

Path to a config file. This is a tab delimited text file
indicating the options for running the Augustus gene prediction
program. Lines beginning with # are ignored.

=back

=head1 OPTIONS

=over 2

=item -p,--program

The program name to use. The default name is AUGUSTUS.

=item --augustus-bin

The path to the augustus binary. The default behavior is to assume that
the binary is in your $PATH such that typing augusuts from the command
line will invoke the program.

=item --gff-ver

The version of GFF to use for output. By default GFF3 output is produced, but
GFF2 can also be produced. Valid values are therefore 'GFF2' or 'GFF3'.

=item -d,--delim

The delimiting character to use between the program and the parameter set name.
By default the colon (:) character is used. Howerver since Gbrowse does not 
allow the use of the colon character for the source field, another value
can be set with this option.

=item --apollo

Produce a GFF file result that can be properly rendered in the Apollo genome
annotation curation program. Since Apollo can not properly render CDS features
or genes as computational evidences, this option modifies the use of CDS
in the output file with the word exon. HOWEVER, THESE FEATURES MAY NOT
ALWAYS BE TRUE BIOLOGICAL EXONS. In some cases the 5' and 3' UTRs are not
represented in the exon structure produced in this conversion.

=item --usage

Short overview of how to use program from command line.

=item --help

Show program usage with summary of options.

=item --version

Show program version.

=item --man

Show the full program manual. This uses the perldoc command to print the 
POD documentation for the program.

=item --verbose

Run the program with maximum output.

=item -q,--quiet

Run the program with minimal output.

=item --test

Run the program without doing the system commands. This will
test for the existence of input files.

=back

=head1 DIAGNOSTICS

Error messages generated by this program and possible solutions are listed
below.

=over 2

=item ERROR: No fasta files were found in the input directory

The input directory does not contain fasta files in the expected format.
This could happen because you gave an incorrect path or because your sequence 
files do not have the expected *.fasta extension in the file name.

=item ERROR: Could not create the output directory

The output directory could not be created at the path you specified. 
This could be do to the fact that the directory that you are trying
to place your base directory in does not exist, or because you do not
have write permission to the directory you want to place your file in.

=back

=head1 CONFIGURATION AND ENVIRONMENT

=head2 Configuration File

The location of the configuration file is indicated by the --config option
at the command line.
This is a tab delimited text file
indicating options for running the Augustus gene prediction program.
Lines beginning with # are ignored, and data are in three
tab delimited columns as shown below:

=over 2

=item Col 1. Parameter set name (ie maize_default)

This is the name for the parameter set that will be used in the GFF
output file.

=item Col 2. Species name (ie. arabidopsis, maize)

This is the suffix which will be added to the end of your blast
output file. You can use this option to set different extensions
for different types of blast. For example *.bln for blastn
output and *.blx for blastx output.

=item Col 2. Additional Augustus command line options.

These are defined in the Augustus documentation.

=back

An example config file:

 #-----------------------------------------------------------+
 # DAWGPAWS AUGUSTUS CONFIG FILE
 #-----------------------------------------------------------+
 # PARAM	        SPECIES    ADDITIONAL
 # NAME          NAME       OPTIONS
 #-----------------------------------------------------------+
 # Maize, default settings without extrinsic information
 maize_def	maize
 # Arab
 arab_def	arabidopsis

=head2 Environment

The following variables may be set in the user environment:

=over

=item $AUGUSTUS_BIN_PATH

The path to the augustus binary

=item $AUGUSTUS_CONFIG_PATH

The path to the directory containing the augustus configuration files.

=back

=head1 DEPENDENCIES

=head2 Required Software

The following software is required for batch_augustus.pl:

=over 2

=item * Augustus

This program is dependent on the Augustus gene prediction software. Augustus
is available from:
L<http://augustus.gobics.de/>

=back

=head2 Required Perl Modules

=over

=item * Getopt::Long

This module is required to accept options at the command line.

=back

=head1 BUGS AND LIMITATIONS

Any known bugs and limitations will be listed here.

=head2 Bugs

If you find a bug with this software, file a bug report on the DAWGPAWS
Sourceforge website: http://sourceforge.net/tracker/?group_id=204962

=head2 Limitations

If this program has known limitations they will be listed here.

=head1 SEE ALSO

The program is part of the DAWGPAWS package of genome
annotation programs. See the DAWGPAWS web page 
( http://dawgpaws.sourceforge.net/ )
or the Sourceforge project page 
( http://sourceforge.net/projects/dawgpaws ) 
for additional information about this package.


=head1 REFERENCE

Please refer to the DAWGPAWS manuscript in Plant Methods when describing
your use of this program:

JC Estill and JL Bennetzen. 2009. 
"The DAWGPAWS Pipeline for the Annotation of Genes and Transposable 
Elements in Plant Genomes." Plant Methods. 5:8.

Use of the Augustus program should also acknowledge the following:

M Stanke, M Diekhans, R Baertsch, D Haussler. 2008.
"Using native and syntenically mapped cDNA alignments to improve de novo 
gene finding." Bioinformatics 24: 637-644. 
doi 10.1093/bioinformatics/btn013

=head1 LICENSE

GNU General Public License, Version 3

L<http://www.gnu.org/licenses/gpl.html>

THIS SOFTWARE COMES AS IS, WITHOUT ANY EXPRESS OR IMPLIED
WARRANTY. USE AT YOUR OWN RISK.

=head1 AUTHOR

James C. Estill E<lt>JamesEstill at gmail.comE<gt>

=head1 HISTORY

STARTED: 06/01/2011

UPDATED: 07/13/2011

VERSION: $Rev$

=cut

