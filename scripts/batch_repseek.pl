#!/usr/bin/perl -w
#-----------------------------------------------------------+
#                                                           |
# batch_repseek.pl - Run repseek in batch mode              |
#                                                           |
#-----------------------------------------------------------+
#                                                           |
#  AUTHOR: James C. Estill                                  |
# CONTACT: JamesEstill_@_gmail.com                          |
# STARTED: 09/06/2012                                       |
# UPDATED: 09/06/2012                                       |
#                                                           |
# DESCRIPTION:                                              |
#  Run the repseek program in batch mode for a given        |
#  configuration file.                                      |
#                                                           |
# USAGE:                                                    |
#  batch_repseek.pl -i indir -o outdir -c config.txt        |
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
my $config_file;

my @repseek_params = ();         # Repminer config options

# BOOLEANS
my $quiet = 0;
my $verbose = 0;
my $show_help = 0;
my $show_usage = 0;
my $show_man = 0;
my $show_version = 0;
my $do_test = 0;                  # Run the program in test mode

#-----------------------------+
# COMMAND LINE OPTIONS        |
#-----------------------------+
my $ok = GetOptions(# REQUIRED OPTIONS
		    "i|indir=s"   => \$indir,
                    "o|outdir=s"  => \$outdir,
                    "c|config=s"  => \$config_file,
		    # ADDITIONAL OPTIONS
		    "q|quiet"     => \$quiet,
		    "verbose"     => \$verbose,
		    # ADDITIONAL INFORMATION
		    "usage"       => \$show_usage,
		    "test"        => \$do_test,
		    "version"     => \$show_version,
		    "man"         => \$show_man,
		    "h|help"      => \$show_help,);

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
    print "\nbatch_mask.pl:\n".
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


print STDERR "NUMBER OF FASTA FILES TO PROCESS: $count_files\n" 
    if $verbose;

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
    next if m/^\#/;
    
    my @in_line = split (/\t/, $_); # Split input by tab
    my $num_in_line = @in_line; 
    
    # Can have just a name to run LTR_Finder with default settings
    # or can have two columns with additional parameter options
    # parameter options in the second columns.
    # I will currently stick with the two column config file
    # since there are so many options availabe with LTR_FINDER
    # that a multiple column config file would get messy.
    if ($num_in_line == 2 ||
	$num_in_line == 3) { 
	$repseek_params[$i][0] = $in_line[0];            # Config name
	$repseek_params[$i][1] = $in_line[1];            # Options
	$i++;
    } # End of if $num_in_line is 10
    else {
	print "\a";
	print STDERR "WARNING: Unexpected number of line in config".
	    " file line $config_line_num\n$config_file\n";
	print STDERR "Found ".$num_in_line." \n";
	print STDERR $in_line[0]."\n";
    }
    
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
print STDERR "$num_proc_total find_ltr runs to process\n\n" if $verbose;

#-----------------------------------------------------------+
# MAIN PROGRAM BODY                                         |
#-----------------------------------------------------------+
my $file_num = 0;
my $proc_num = 0;
for my $ind_file (@fasta_files) {

    $file_num++;
    my $name_root;
    # Get root file name
    if ($ind_file =~ m/(.*)\.hard\.fasta$/) {
	# file ends in .hard.fasta
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
    my $name_root_dir = $outdir.$name_root."/";
    unless (-e $name_root_dir) {
	mkdir $name_root_dir ||
	    die "Could not create dir:\n$name_root_dir\n";
	print STDERR "Making root name dir\n";
	print STDERR "\t".$name_root_dir."\n";
    }

    #-----------------------------+
    # CREATE REPSEEK OUTDIR       |
    #-----------------------------+
    # Dir to hold gene prediction output from local software
    my $repseek_dir = $name_root_dir."repseek/";
    unless (-e $repseek_dir) {
	mkdir $repseek_dir ||
	    die "Could not create ltr_finder out dir:\n$repseek_dir\n";
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
	
	my $rs_param_name = $repseek_params[$i][0];
	my $rs_options = $repseek_params[$i][1];
	my $rs_outfile = $repseek_dir.$name_root.".".$rs_param_name.".txt";
	
	my $rs_cmd = "repseek".
	    " ".$rs_options.
	    " -r ".$rs_outfile.
	    " ".$indir.$ind_file;
	
	system ($rs_cmd) unless $do_test;
	
	print STDERR $rs_cmd."\n";
	
    }
    

}

exit 0;

#-----------------------------------------------------------+ 
# SUBFUNCTIONS                                              |
#-----------------------------------------------------------+

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

Name.pl - Short program description. 

=head1 VERSION

This documentation refers to program version 0.1

=head1 SYNOPSIS

  USAGE:
    Name.pl -i InDir -o OutDir

    --infile        # Path to the input file
    --outfie        # Path to the output file

=head1 DESCRIPTION

This is what the program does

=head1 REQUIRED ARGUMENTS

=over 2

=item -i,--indir

Path of the directory containing the sequences to process.

=item -o,--outdir

Path of the directory to place the program output.

=item -c, --config

Path to a config file. This is a tab delimited text file
indicating the required information for each of the databases to blast
against. Lines beginning with # are ignored.

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

=item --verbose

Run the program with maximum output.

=item -q,--quiet

Run the program with minimal output.

=item --test

Run the program without doing the system commands. This will
test for the existence of input files.

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
indicating required information for each of the databases to blast
against. Lines beginning with # are ignored, and data are in six 
columns as shown below:

=over 2

=item Col 1. Blast program to use [ tblastx | blastn | blastx ]

The blastall program to use. DAWG-PAWS will support blastn,
tblastx, and blastx format.

=item Col 2. Extension to add to blast output file. (ie. bln )

This is the suffix which will be added to the end of your blast
output file. You can use this option to set different extensions
for different types of blast. For example *.bln for blastn
output and *.blx for blastx output.

=back

An example config file:

 #-----------------------------+
 # BLASTN: TIGR GIs            |
 #-----------------------------+
 blastn	bln	8	1e-5	TaGI_10	-a 2 -U
 blastn	bln	8	1e-5	AtGI_13	-a 2 -U
 blastn	bln	8	1e-5	ZmGI_17	-a 2 -U
 #-----------------------------+
 # TBLASTX: TIGR GIs           |
 #-----------------------------+
 tblastx	blx	8	1e-5	TaGI_10	-a 2 -U
 tblastx	blx	8	1e-5	AtGI_13	-a 2 -U
 tblastx	blx	8	1e-5	ZmGI_17	-a 2 -U

=head2 Environment

This program does not make use of variables in the user environment.

=head1 DEPENDENCIES

=head2 Required Software

=over

=item * Software Name

Any required software will be listed here.

=back

=head2 Required Perl Modules

=over

=item * Getopt::Long

This module is required to accept options at the command line.

=back

=head1 BUGS AND LIMITATIONS

Any known bugs and limitations will be listed here.

=head2 Bugs

=over 2

=item * No bugs currently known 

If you find a bug with this software, file a bug report on the DAWG-PAWS
Sourceforge website: http://sourceforge.net/tracker/?group_id=204962

=back

=head2 Limitations

=over

=item * Known Limitation

If this program has known limitations they will be listed here.

=back

=head1 SEE ALSO

The program is part of the DAWG-PAWS package of genome
annotation programs. See the DAWG-PAWS web page 
( http://dawgpaws.sourceforge.net/ )
or the Sourceforge project page 
( http://sourceforge.net/projects/dawgpaws ) 
for additional information about this package.

=head1 LICENSE

GNU General Public License, Version 3

L<http://www.gnu.org/licenses/gpl.html>

THIS SOFTWARE COMES AS IS, WITHOUT ANY EXPRESS OR IMPLIED
WARRANTY. USE AT YOUR OWN RISK.

=head1 AUTHOR

James C. Estill E<lt>JamesEstill at gmail.comE<gt>

=head1 HISTORY

STARTED: 09/06/2012

UPDATED: 09/06/2012

VERSION: $Rev$

=cut


#-----------------------------------------------------------+
# REPSEEK USAGE
#-----------------------------------------------------------+
#
#General usage is 'repseek [-v|-h|-H] [ opts ] { core_value } file.fst [file2.fst]'
#
#   Core values (select either at least one of the four -l/-p  -L/-P option)
#      seed minimum length: lmin
#	-l lmin      : set a minimum length for the seeds detection (no pvalue).
#	-p prob      : set a p-value that gives a minimum seed length (Karlin & Ost equation)
#      repeat minimum score: smin
#	-L smin      : set a minimum score for extended repeat (def is 0: no filter).
#	-P prob      : set a p-value that gives a minimum repeat score (Waterman & Vingron regression)
#
#   Optionnal values
#	-r file      : 'r'edirect output into file (otherwise stdout is used).
#	-T           : print out R-'T'able
#	               the R-Table (or R-Table2) shows all non-unique position and its degree of redundancy
#	-S           : output 'S'eeds and exit (do not perform extention).
#	-s seed_file : use seeds given in seed_file instead of detecting them
#	               File format is 'd|i begin end len [optionnal other fields]'.
#	-D mask_file : mask sequence regions during seed detection only (cannot be used with -s seed_file).
#	               File format is 'begin end [seq#]' (seq# is 1 or 2; nothing is treated as 1).
#	-R 0.##      : merge 'R'epeats when they share 0.## of their both copies (default = 0.90).
#	               When set to 1.0, merge repeats having exactly the same positions
#	-d           : detect only 'd'irect repeats.
#	-i           : detect only 'i'nverted repeats.
#	               (-d and -i are mutually exclusive).
#	-B           : if some seeds overlaps (i.e. in low-complexity sequence), just keep the 'B'iggest one.
#	-m #.#       : keep only seeds occurring at least the specified minimum number of times (min is a float).
#	-M #.#       : keep only seeds occurring less than the specified maximum number of times (Max is a float).
#	-O 0|1       : for 1 sequence, direct repeat can have their 2 copies 'O'verlapping (0: no, 1: yes -default-)
#	-c           : set the chromosome as 'c'ircular -default is linear- (unused for 2 sequences)
#	-o #.#       : set gap_open penalty -default 4.0- (cannot be change when -P is used)
#	-e #.#       : set gap_ext penalty -default 1.0- (cannot be change when -P is used)
#	-X #.#       : set 'X'g, the 'exploration value' -default 30.0-
#
#   Information
#	-v           : Print 'v'ersion and exit
#	-h           : Print a short 'h'elp and exit
#	-H           : Print a larger 'H'elp on a specific mode and exit
#
#   Input
#     Input file(s) should be in a fasta format.
#     All non-ACGT characters are considered as N, except for X in which no repeats are detected (mask)
#
#   Contact
#     if you need help, want to report a bug or suggest an option,
#     please contact G. Achaz by email: achaz(at)abi.snv.jussieu.fr

