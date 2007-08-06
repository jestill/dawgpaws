#!/usr/bin/perl -w
#-----------------------------------------------------------+
#                                                           |
# cnv_blast2gff.pl Convert BLAST output to gff 
#                                                           |
#-----------------------------------------------------------+
#  AUTHOR: James C. Estill                                  |
# CONTACT: JamesEstill_at_gmail.com                         |
# STARTED: 07/06/2006                                       |
# UPDATED: 08/06/2007                                       |
#                                                           |  
# DESCRIPTION:                                              | 
# Convert blast output to a Apollo compatible gff file.     |
#                                                           |
# USAGE:                                                    |
#                                                           |
#-----------------------------------------------------------+

=head1 NAME

batch_cnv_blast2gff.pl - Short program description. 

=head1 VERSION

This documentation refers to program version 0.1

=head1 SYNOPSIS

 Usage:
  Name.pl -i InDir -o OutDir

=head1 DESCRIPTION

This is what the program does

=head1 REQUIRED ARGUMENTS

=over 2

=item -i,--infile

Path of the input file.

=item -o,--outfile

Path of the output file.

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

Other modules or software that the program is dependent on.

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
use strict;                    # Follow the rules
use Getopt::Long;              # Get options from the command line
use Bio::SearchIO;             # Parse BLAST output
use File::Copy;                # Copy the gff output to the gff dir

#-----------------------------+
# PROGRAM VARIABLES           |
#-----------------------------+
my $VERSION = "0.1";

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
                    "o|outdir=s"   => \$outdir,
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
    else {
	$name_root = "UNDEFINED";
    }

    # PRINT PROCESS STATUS TO LOG
    print LOG "\n\nProcessing File $file_num of $num_files.\n" if $logfile;

    # PRINT PROCESS STATUS TO TERMINAL
    print "\n\n+-----------------------------------------------------------+\n"
	unless $quiet;
    print "| Processing File $file_num of $num_files.\n" unless $quiet;
    print "+-----------------------------------------------------------+\n"
	unless $quiet;
    
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

	
	print "IN:\t$ind_blast_dir\n";
	print "OUT:\t$out_gff_path\n";
	print "COPY:\t$out_gff_copy\n\n";

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

	    blast2gff ( $blast_file_path, $out_gff_path, $do_append, $name_root);
	}

	
	#-----------------------------+
	# MAKE COPY OF GFF FILE IN    |
	# GFF DIR                     |
	#-----------------------------+
	$msg = "ERROR: Could not copy file from:\n".
	    "\t$out_gff_path\n".
	    "\t$out_gff_copy\n";
	copy ($out_gff_path, $out_gff_copy)
	    || print LOG "msg";

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
	"  --infile       # Path to the input file\n".
	"  --outfile      # Path to the output file\n".
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


=head1 HISTORY

STARTED: 08/06/2007

UPDATED: 08/06/2007

VERSION: $Id: $

=cut

#-----------------------------------------------------------+
# HISTORY                                                   |
#-----------------------------------------------------------+
#
# 08/06/2007
# -Program started
