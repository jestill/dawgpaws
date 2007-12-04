#!/usr/bin/perl -w
#-----------------------------------------------------------+
#                                                           |
# cnv_blast2gff.pl Convert BLAST output to gff 
#                                                           |
#-----------------------------------------------------------+
#  AUTHOR: James C. Estill                                  |
# CONTACT: JamesEstill_at_gmail.com                         |
# STARTED: 07/06/2006                                       |
# UPDATED: 12/04/2007                                       |
#                                                           |  
# DESCRIPTION:                                              | 
# Convert blast output to a Apollo compatible gff file.     |
#                                                           |
# LICENSE:                                                  |
#  GNU General Public License, Version 3                    |
#  http://www.gnu.org/licenses/gpl.html                     |  
#                                                           |
#-----------------------------------------------------------+

#-----------------------------+
# INCLUDES                    |
#-----------------------------+
use strict;                    # Follow the rules
use Getopt::Long;              # Get options from the command line
use Bio::SearchIO;             # Parse BLAST output
use File::Copy;                # Copy the gff output to the gff dir
use Pod::Select;               # Print subsections of POD documentation
#use Pod::Text;                 # Print POD doc as formatted text file
use IO::String;
use IO::Scalar;
use IO::File;
use IO::Pipe;
use Pod::Html qw( &pod2html );
use Pod::Text qw( &pod2text );

#-----------------------------+
# PROGRAM VARIABLES           |
#-----------------------------+
my ($VERSION) = q$Rev$ =~ /(\d+)/;

#-----------------------------+
# VARIABLE SCOPE              |
#-----------------------------+

# Set variable scope
my $indir;
my $outdir;
my $logfile;
my $name_root;
my $out_gff_path;
my $out_gff_dir;
my $out_gff_copy;
my $msg;
my $ind_blast_dir;
my $blast_file_num;

# Booleans
my $verbose = 0;
my $help = 0;
my $quiet = 0;
my $show_help = 0;
my $show_usage = 0;
my $show_man = 0;
my $show_version = 0;
my $do_append = 0;

# Arrays
my @blast_files;


#//////////////////////
# file_num_max is the number of seqs to process in test run
my $file_num_max = 1;
my $file_num = 0;
#\\\\\\\\\\\\\\\\\\\\\\


#-----------------------------+
# COMMAND LINE OPTIONS        |
#-----------------------------+
my $ok = GetOptions(# REQUIRED
		    "i|indir=s"   => \$indir,
                    "o|outdir=s"  => \$outdir,
		    # OPTIONS
		    "logfile=s"   => \$logfile,
		    # BOOLEANS
		    "verbose"     => \$verbose,
		    "append"      => \$do_append,
		    #"maxe"        => \$max_e,
		    #"minlen"      => \$min_len,
		    "usage"       => \$show_usage,
		    "version"     => \$show_version,
		    "man"         => \$show_man,
		    "h|help"      => \$show_help,
		    "q|quiet"     => \$quiet,);

#-----------------------------+
# SHOW REQUESTED HELP         |
#-----------------------------+
#if ($show_usage) {
#    print_help("");
#}


# Printing Help with POD::Select
if ($show_usage) {
    
    print "\n";
    
    #-----------------------------+
    # PIPE WITHIN PERL            |
    #-----------------------------+
    # This code made possible by:
    # http://www.perlmonks.org/index.pl?node_id=76409
    # Tie info developed on:
    # http://www.perlmonks.org/index.pl?node=perltie 
    #
    my $podfile = $0;
    my $scalar = '';
    tie *STDOUT, 'IO::Scalar', \$scalar;
    podselect({-sections => ["NAME|SYNOPSIS|MORE"]}, $0);
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
	# BEGIN AT WORK HERE
	#pod2html();	
	pod2text();
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

    #-----------------------------+
    # PRINT POD SECTION           |
    #-----------------------------+
    # The following prints NAME and SYNOPSIS
    #podselect({-sections => ["NAME|SYNOPSIS|MORE"]}, $0);

    #-----------------------------+
    # PRINT ENTIRE POD
    #-----------------------------+
    #podselect ($0);

    #-----------------------------+
    # PRINT RAW POD               |
    # OBJ ORIENTED                |
    #-----------------------------+
    #my $parser = new Pod::Select(-output => ">&STDERR");
    #my $pod_txt = Pod::Text->new (sentence => 0, width => 78);
    #$parser->select("NAME|SYNOPSIS|BLANK");
    #$parser->parse_from_file($0);

    #-----------------------------+
    # USE SYSTEM COMMANDS         |
    #-----------------------------+
    #my $usage_cmd = "podselect -section 'SYNOPSIS|MORE' $0 | pod2text";
    #system ($usage_cmd);

    print "\n";

    exit 0;

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
# CHECK REQUIRED OPTIONS      |
#-----------------------------+
if ( (!$indir) || (!$outdir) ) {
    print_help("full");
}

#-----------------------------+
# OPEN THE LOG FILE           |
#-----------------------------+
if ($logfile) {
    # Open file for appending
    open ( LOG, ">>$logfile" ) ||
	die "Can not open logfile:\n$logfile\n";
    my $time_now = time;
    print LOG "==================================\n";
    print LOG "  batch_mask.pl\n";
    print LOG "  JOB: $time_now\n";
    print LOG "==================================\n";
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

my $num_files = @fasta_files;

#-----------------------------+
# SHOW ERROR IF NO FILES      |
# WERE FOUND IN THE INPUT DIR |
#-----------------------------+
if ($num_files == 0) {
    print "\a";
    print "\nERROR: No fasta files were found in the input directory\n".
	"$indir\n".
	"Fasta files must have the fasta or fa extension.\n\n";
    exit;
}

#-----------------------------+
# CREATE THE OUT DIR          |
# IF IT DOES NOT EXIST        |
#-----------------------------+
unless (-e $outdir) {
    print "Creating output dir ...\n" unless $quiet;
    mkdir $outdir
	|| die "Could not create the output directory:\n$outdir";
}


#-----------------------------+
# CONVERT BLAST TO GFF FOR    |
# EACH FILE IN THE DIR        |
#-----------------------------+
for my $ind_file (@fasta_files)
{

    $file_num++;
    $blast_file_num = 0;

    #-----------------------------+
    # GET BASE FILE NAME          |
    #-----------------------------+
    if ($ind_file =~ m/(.*)\.masked\.fasta$/ ) {	    
	$name_root = "$1";
    }  
    elsif ($ind_file =~ m/(.*)\.fasta$/ ) {	    
	$name_root = "$1";
    }
    if ($ind_file =~ m/(.*)\.masked\.fa$/ ) {	    
	$name_root = "$1";
    }
    elsif ($ind_file =~ m/(.*)\.fa$/ ) {	    
	$name_root = "$1";
    } 
    else {
	$name_root = "UNDEFINED";
    }

    # PRINT PROCESS STATUS TO LOG
    print LOG "\n\nProcessing File $file_num of $num_files.\n" if $logfile;

    # PRINT PROCESS STATUS TO TERMINAL
    print STDERR 
	"\n\n+-----------------------------------------------------------+\n"
	if $verbose;
    print STDERR 
	"| Processing File $file_num of $num_files.\n" if $verbose;
    print STDERR 
	"+-----------------------------------------------------------+\n"
	if $verbose;
    
    #-----------------------------+
    # LOAD LIST OF BLAST FILES    |
    # TO PROCESS                  |
    #-----------------------------+ 
    $ind_blast_dir = $outdir.$name_root."/blast/";
    if (-e $ind_blast_dir) {
	opendir( BLASTDIR, $ind_blast_dir ) || 
	    die "Can't open directory:\n$ind_blast_dir"; 
	my @blast_files = grep /\.blo$|\.bln$|\.blx$/, readdir BLASTDIR ;
	closedir( BLASTDIR );


	#-----------------------------+
	# SET PATH TO OUTPUT GFF FILES|
	#-----------------------------+
	$out_gff_path = $ind_blast_dir.$name_root."_all_blast.gff";
	
	$out_gff_dir = $outdir.$name_root."/gff/";
	unless (-e $out_gff_dir) {
	    print "Creating output dir ...\n$out_gff_dir" if $verbose;
	    mkdir $out_gff_dir
		|| die "Could not create the output directory:\n$out_gff_dir";
	}
	
	$out_gff_copy = $out_gff_dir.$name_root."_all_blast.gff";

	
	print STDERR "IN:\t$ind_blast_dir\n";
	print STDERR "OUT:\t$out_gff_path\n";
	print STDERR "COPY:\t$out_gff_copy\n\n";

	#-----------------------------+
	# CONVERT EACH BLAST OUTFILE  |
	# TO GFF                      |
	#-----------------------------+
	$blast_file_num = 0;

	for my $ind_blast_file (@blast_files) {
	    $blast_file_num++; 

	    # For first blast output file overwrite any existing data
	    # This is done by setting the do_append boolean to 0
	    if ($blast_file_num == 1) {
		$do_append = 0;
	    }
	    else {
		$do_append = 1;
	    }
	    
	    my $blast_file_path = $ind_blast_dir.$ind_blast_file;
	    print "Converting: $ind_blast_file\n";

	    blast2gff ( $blast_file_path, $out_gff_path, 
			$do_append, $name_root);
	}

	
	#-----------------------------+
	# MAKE COPY OF GFF FILE IN    |
	# GFF DIR                     |
	#-----------------------------+
	$msg = "ERROR: Could not copy file from:\n".
	    "\t$out_gff_path\n".
	    "\t$out_gff_copy\n";
	copy ($out_gff_path, $out_gff_copy)
	    || print LOG "$msg";

    }
    else {
	print LOG "ERROR: Could not find BLAST dir:\n$ind_blast_dir\n";
    }

    # If max file number then exit
    # This is for debug tests
    #if ($file_num = $file_num_max) {exit;}

}

close LOG if $logfile;

exit;



#-----------------------------------------------------------+
# SUBFUNCTIONS                                              |
#-----------------------------------------------------------+

sub blast2gff {
# CONVERT BLAST TO GFF 
    
    # seqname - ID to use in the seqname field
    # blastin - path to the blast input file
    # gffout  - path to the gff output file
    # append  - boolean append data to existing file at gff out
    my ($blastin, $gffout, $append, $seqname) = @_;
    my $blastprog;        # Name of the blast program used (blastn, blastx)
    my $dbname;           # Name of the database blasted
    my $hitname;          # Name of the hit
    my $start;            # Start of the feature
    my $end;              # End of the feature
    my $strand;           # Strand of the hit

    # Open the BLAST report object
    my $blast_report = new Bio::SearchIO ( '-format' => 'blast',
					   '-file'   => $blastin) 
	|| die "Could not open BLAST input file:\n$blastin.\n";
    
    # Open file handle to the gff outfile    
    if ($append) {
	open (GFFOUT, ">>$gffout") 
	    || die "Can not open file:\n $gffout\n";
    }
    else {
	open (GFFOUT, ">$gffout") 
	    || die "Can not open file:\n $gffout\n";
    }
    
    while (my $blast_result = $blast_report->next_result())
    {
	$blastprog = $blast_result->algorithm;
	$dbname = $blast_result->database_name();

    	while (my $blast_hit = $blast_result->next_hit())
	{

	    while (my $blast_hsp = $blast_hit->next_hsp())
	    {

		my $hitname = $blast_hit->name();

		$strand = $blast_hsp->strand('query');
		
		if ($strand =~ "-1") {
		    $strand = "-";
		}
		elsif ($strand =~ "1") {
		    $strand = "+";
		}
		else {
		    die "Error parsing strand\n";
		}
		
		# Make certain that start coordinate is
		# less then the end coordinate
		if ( $blast_hsp->start() < $blast_hsp->end() ) {
		    $start = $blast_hsp->start();         # Start
		    $end = $blast_hsp->end();
		    #$strand = "+";
		    
		}
		else {
		    #
		    $start = $blast_hsp->end();         # Start
		    $end = $blast_hsp->start();
		    #$strand = "-";
		}


		#-----------------------------+
		# PRINT OUTPUT TO GFF FILE    |
		# WITH BAC DATA               |
		#-----------------------------+
		# Changing BLASTN to the Bac Name appears to allow 
		# these to be drawn on different levels.
		print GFFOUT 
		    #$blast_hit->name()."\t".
		    "$seqname\t".                            # Seqname
		    "$blastprog:$dbname\t".                  # Source
		    "exon\t".                                # Feature type name
		    "$start\t".                              # Start
		    "$end\t".                                # End
		    $blast_hsp->score()."\t".                # Score
		    "$strand\t".                             # Strand
		    ".\t".                                   # Frame
		    "$hitname\n";                            # Feature name

	    } # End of while next hsp
	} # End of while next hit
    } # End of while next result
    
    close GFFOUT;
    
}

sub print_help {

    # Print requested help or exit.
    # Options are to just print the full 
    my ($opt) = @_;

    my $usage = "USAGE:\n". 
	"cnv_blast2gff.pl -i InFile -o OutFile\n";
    my $args = "REQUIRED ARGUMENTS:\n".
	"  --indir       # Path to the input directory\n".
	"  --outdir      # Path to the output direcotry\n".
	"\n".
	"OPTIONS::\n".
	"  --append       # Append results to an existing file\n".
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

1;
__END__

=head1 NAME

batch_cnv_blast2gff.pl - Convert blast output to gff format

=head1 VERSION

This documentation refers to program version $Rev$

=head1 SYNOPSIS

=head2 Usage

    batch_cnv_blast2gff.pl -i InDir -o OutDir

=head2 Required Arguments

    -i, --indir    # Directory of fasta files to process
    -o, --outdir   # Path to the base output directory

=head1 DESCRIPTION

Convert NCBI BLAST output to an Apollo compatible GFF file. This can produce
a separate gff file for each BLAST report, or can merge all BLAST results
into a single GFF file.

=head1 REQUIRED ARGUMENTS

=over 2

=item -i,--indir

Path of the directory containing the sequences to process.

=item -o,--outdir

Path of the directory to place the program output.

=back

=head1 OPTIONS

=over 2

=item --logfile

Path to a file that will be used to log program status.
If the file already exists, additional information will be concatenated
to the existing file

=item --append

Append results to existing gff file.

=item --usage

Print a short overview of how to use program from the command line.

=item --help

Print a short program usage state with a summary of options.

=item --version

Show program version. This will print the SVN version of the script.

=item --man

Show the full program manual. This uses the perldoc command to print the 
POD documentation for the program.

=item -q,--quiet

Run the program with minimal output.

=item --verbose

Run the program with maximum output to the screen. This will provided a 
detailed status of the progress of the program.

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

This program does not depend on any configuration files or environmental
settings.

=head1 DEPENDENCIES

=head2 Required Perl Modules

=over

=item * Bio::SearchIO;

This module is part of the BioPerl package of programs. It is used to
parse the BLAST output. For information on downloading and installing
BioPerl see http://www.bioperl.org/wiki/Main_Page

=item * File::Copy

This module is required to copy the BLAST results.

=item * Getopt::Long

This module is required to accept options at the command line.

=back

=head1 BUGS AND LIMITATIONS

=head2 Bugs

=over 2

=item * No bugs currently known 

If you find a bug with this software, file a bug report on the DAWG-PAWS
Sourceforge website: http://sourceforge.net/tracker/?group_id=204962

=back

=head2 Limitations

=over

=item * Limited to NCBI BLAST

The current version is limited to using the NCBI version of BLAST.

=item * Limited file extensions are supported

BLAST output file must currently end with blo, bln, or blx. For example
a BLASTx output may be named BlastOut.blx while a BLASTN output
may be names BlastOut.bln. FASTA files must end with a fasta or fa extension.
For examples must have names like my_seq.fasta or my_seq.fa.

=back

=head1 SEE ALSO

The batch_blast.pl program is part of the DAWG-PAWS package of genome
annotation programs. See the DAWG-PAWS web page 
( http://dawgpaws.sourceforge.net/ )
or the Sourceforge project page 
( http://sourceforge.net/projects/dawgpaws ) 
for additional information about this package.

=head1 LICENSE

GNU GENERAL PUBLIC LICENSE, VERSION 3

http://www.gnu.org/licenses/gpl.html   

=head1 AUTHOR

James C. Estill E<lt>JamesEstill at gmail.comE<gt>

=head1 HISTORY

STARTED: 08/06/2007

UPDATED: 08/06/2007

VERSION: $Rev$

=cut

#-----------------------------------------------------------+
# HISTORY                                                   |
#-----------------------------------------------------------+
#
# 08/06/2007
# -Program started, all base subfunctions added
#
# 12/04/2007
# - POD Documentation updated
# - Version number changed to SVN Revision Number
#
# 12/04/2007
# - Changed terminal output to STDERR
# - Moved POD documentation to end of the code
# - Trying to add POD select to print help and usage
#   message from the POD documentation
# - Addd an end statement
