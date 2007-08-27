#!/usr/bin/perl -w
#-----------------------------------------------------------+
#                                                           |
# batch_audit.pl - Audit DAWG PAWS analysis directory       |
#                                                           |
#-----------------------------------------------------------+
#                                                           |
#  AUTHOR: James C. Estill                                  |
# CONTACT: jestill_at_sourceforge.net                       |
# STARTED: 07/30/2007                                       |
# UPDATED: 08/02/2007                                       |
#                                                           |
# DESCRIPTION:                                              |
#  Audit the dawg paws analysis output.                     |
#                                                           |
#-----------------------------------------------------------+

=head1 NAME

batch_audit.pl - Audit the DAWG PAWS analysis directory

=head1 VERSION

This documentation refers to batch_audit version 1.0

=head1 SYNOPSIS

 Usage:
 batch_audit.pl -i DirToProcess -o OutDir

=head1 DESCRIPTION

Audit the dawg paws analysis output.

=head1 REQUIRED ARGUMENTS

=over 2

=item -i,--indir

Path of the directory containing the sequences to process.

=item -o,--outdir

Path of the directory to place the program output.

=back

=head1 OPTIONS

=over 2

=item --color

Print the terminal output in fancy ANSI color.

=item --deep

Do a deep audit. 
This process is slower and includes the following additional steps: 

=over 2

=item *

BLAST output files that are empty are moved to a subdir

=back

=item --blast

Audit the BLAST output

=item --rm

Audit the repeatmasker output

=item --ta

Audit the TriAnnotation output

=item --full

Run a full audit. This will currently audit blast output, 
repeatmasker output and TriAnnotation output.

=item --logfile

Path to a file that will be used to log program status.
If the file already exists, additional information will be concatenated
to the existing file.

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

=item --test

Run the program without doing the system commands.

=back

=head1 DIAGNOSTICS

The error messages that can be generated will be listed here.

=head1 CONFIGURATION AND ENVIRONMENT

Names and locations of config files
environmental variables
or properties that can be set.

=head1 DEPENDENCIES

=head2 Required Software

=over

=item *

RepeatMasker
(http://www.repeatmasker.org/)

=item *

Apollo (Genome Annotation Curation Tool)
http://www.fruitfly.org/annot/apollo/

=back

=head2 Required Perl Modules

=over

=item *

File::Copy

=item *

Getopt::Long

=back

=head1 BUGS AND LIMITATIONS

=head2 TO DO

=over 2

=item *

Make the results compatable for an upload to a chado
database.

=item *

Make it a variable to possible to put the gff output (1) all in positive 
strand, (2) all in negative strand, (3) alignment to positive or
negative strand, (4) cumulative in both positive and negative strand.
Current behavior will be to do number 4 above.

=back

=head2 Limitations

=over

=item *

Currently must use short names in the FASTA file.

=item *

This program has been tested with RepeatMasker v  3.1.6

=back

=head1 LICENSE

GNU LESSER GENERAL PUBLIC LICENSE

http://www.gnu.org/licenses/lgpl.html

=head1 AUTHOR

James C. Estill E<lt>JamesEstill at gmail.comE<gt>

=cut

print "\n";

#-----------------------------+
# INCLUDES                    |
#-----------------------------+
use File::Copy;                # Used to move files
use Getopt::Long;              # Get cmd line options
use Bio::SearchIO;             # Parse BLAST output

#-----------------------------+
# PROGRAM VARIABLES           |
#-----------------------------+
my $ver = "1.0";

#-----------------------------+
# VARIABLE SCOPE              |
#-----------------------------+
my $logfile;                   # Path to a logfile to log error info
my $indir;                     # Directory containing the seq files to process
my $outdir;                    # Directory to hold the output
my $msg;                       # Message printed to the log file
my $name_root;                 # Root name to be used for output etc

# ARRAYS
my @pass_files;              # Array to hold the list of passed seqs

# BOOLEANS
my $show_help = 0;             # Show program help
my $show_version = 0;          # Show program version
my $show_man = 0;              # Show program manual page using peldoc
my $show_usage = 0;            # Show program usage command             
my $quiet = 0;                 # Boolean for reduced output to STOUT
my $apollo = 0;                # Path to apollo and apollo variables
my $test = 0;                  # Run the program in test mode
my $verbose = 0;               # Run the program in verbose mode
my $do_audit_blast = 0;        # Audit the blast output
my $do_audit_rm = 0;           # Audit the repeatmasker output
my $do_audit_ta = 0;           # Audit the TriAnnotation output
my $do_full_audit = 0;         # Audit everything
my $blast_ok = 0;              # Blast output is okay
my $ta_ok = 0;                 # TriAnnotation output is okay
my $rm_ok = 0;                 # RepeatMask outoupt is okay
my $print_color = 0;           # Print output in fancy ANSI color
my $pass_audit = 1;            # Does the Seq pass the audit 
my $do_deep_audit = 0;         # Do a deep audit, move empty BLAST etc

# COUNTERS
my $num_pass = 0;              # Number of seqs that pass the audit
my $num_proc = 1;              # Number of processors to use

#-----------------------------+
# COMMAND LINE OPTIONS        |
#-----------------------------+
my $ok = GetOptions(
		    # Required
		    "i|indir=s"    => \$indir,
                    "o|outdir=s"   => \$outdir,
		    # Optional strings
		    "logfile=s"    => \$logfile,
		    # Booleans
		    "deep"         => \$do_deep_audit,
		    "color"        => \$print_color,
		    "copy-gff"     => \$copy_gff,
		    "blast"        => \$do_audit_blast,
		    "rm"           => \$do_audit_rm,
		    "ta"           => \$do_audit_ta,
		    "full"         => \$do_full_audit,
		    "verbose"      => \$verbose,
		    "test"         => \$test,
		    "usage"        => \$show_usage,
		    "version"      => \$show_version,
		    "man"          => \$show_man,
		    "h|help"       => \$show_help,
		    "q|quiet"      => \$quiet,);


my $ProcNum = 0;

#//////////////////////
my $file_num_max = 2;
my $file_num = 0;
#\\\\\\\\\\\\\\\\\\\\\\

#-----------------------------+
# SHOW REQUESTED HELP         |
#-----------------------------+
if ($show_usage) {
    print_help("");
}


if ($show_man) {
    # User perldoc to generate the man documentation.
    system ("perldoc $0");
    exit($ok ? 0 : 2);
}

if ($show_help || (!$ok) ) {
    print_help("full");
}

if ($show_version) {
    print "\nbatch_audit.pl:\n".
	"Version: $ver\n\n";
    exit;
}

# Show full help when required options
# are not present
if ( (!$indir) || (!$outdir) ) {
    print_help("full");
}

if ($print_color) {
    use Term::ANSIColor;    
}

#-----------------------------+
# DETERMINE AUDIT STATUS      |
#-----------------------------+
if ($do_full_audit) {
    $do_audit_blast = 1;
    $do_audit_rm = 1;
    $do_audit_ta = 1;
}

#-----------------------------+
# OPEN THE LOG FILE           |
#-----------------------------+
if ($logfile) {
    # Open file for appending
    open ( LOG, ">$logfile" ) ||
	die "Can not open logfile:\n$logfile\n";
    my $time_now = time;
    print LOG "==================================\n";
    print LOG "  batch_audit.pl\n";
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
print "Identifying fasta files\n" if $verbose;
opendir( DIR, $indir ) || 
    die "Can't open directory:\n$indir"; 
my @fasta_files = grep /\.fasta$|\.fa$/, readdir DIR ;
closedir( DIR );

my $num_files = @fasta_files;
@fasta_files = sort(@fasta_files);

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
# EXIT IF THE OUTPUT DIR      |
# DOES NOT EXIST              |
#-----------------------------+
unless (-e $outdir) {
    print "\a";
    print "ERROR: The output dir does not exist at:\n$outdir\n";
    exit;
}

#-----------------------------+
# RUN REQUESTED AUDIT FOR     |
# EACH FASTA FILE             |
#-----------------------------+

for my $ind_file (@fasta_files)
{

    $file_num++;

    # Temp exit while debuging
    #if ($file_num > $file_num_max) {exit;}
    
    if ($ind_file =~ m/(.*)\.hard\.fasta$/ ) {	    
	$name_root = "$1";
    } 
    elsif ($ind_file =~ m/(.*)\.masked\.fasta$/ ) {	    
	$name_root = "$1";
    }  
    elsif ($ind_file =~ m/(.*)\.fasta$/ ) {	    
	$name_root = "$1";
    } 
    else {
	$name_root = "UNDEFINED";
    }

    # Print header to logfile
    if ($logfile) {
	print LOG "=========================================\n";
	print LOG " $name_root\n";
	print LOG "=========================================\n";
    }


    print "=========================================\n";
    print "Processing $name_root\n";
    print "=========================================\n\n";

    #-----------------------------+
    # MAKE SURE THE GFF AND GAME  |
    # DIRS ARE PRESENT            |
    #-----------------------------+
    $dir_game_out = $outdir.$name_root."/game/";
    unless (-e $dir_game_out) {
	print "Creating dir:\n$dir_game_out\n" if $verbose;
	mkdir $dir_game_out ||
	    die "Could not create dir:\n$dir_game_out\n";
    }
    
    $dir_gff_out = $outdir.$name_root."/gff/";
    unless (-e $dir_gff_out) {
	print "Creating dir:\n$dir_gff_out\n" if $verbose;
	mkdir $dir_gff_out ||
	    die "Could not create dir:\n$dir_gff_out\n";
    }

    #-----------------------------+
    # AUDIT REPEAT MASKER OUTPUT  |
    #-----------------------------+
    if ($do_audit_rm) {
	#print "Auditing the repeatmasker output\n" if $verbose;
	$rm_ok = audit_rm($outdir, $name_root);
	
	if ($rm_ok) {
	    print LOG "$name_root RepeatMasker complete\n" if $logfile;
	    print color 'green' if $print_color;
	    print "$name_root RepeatMasker complete\n";
	    print color 'reset' if $print_color;
	} 
	else {
	    $pass_audit = 0;
	    print LOG "$name_root RepeatMasker incomplete\n" if $logfile;
	    print color 'red' if $print_color;
	    print "$name_root RepeatMasker incomplete\n";
	    print color 'reset' if $print_color;
	}
    }


    #-----------------------------+
    # AUDIT BLAST OUTPUT          |
    #-----------------------------+
    if ($do_audit_blast) {
	#print "Auditing BLAST for $name_root\n" if $verbose;
	$blast_ok = audit_blast($outdir, $name_root, $do_deep_audit);

	if ($blast_ok) {
	    print color 'green' if $print_color;
	    print "$name_root BLAST complete\n";
	    print color 'reset' if $print_color;
	}
	else {
	    $pass_audit = 0;
	    print color 'red' if $print_color;
	    print "$name_root BLAST incomplete\n";
	    print color 'reset' if $print_color;
	}

    }

    #-----------------------------+
    # AUDIT TRIANNOTATION OUTPUT  |
    #-----------------------------+
    if ($do_audit_ta) {
	#print "Auditing TriAnnotation output\n" if $verbose;
	$ta_ok = audit_ta($outdir, $name_root);

	if ($ta_ok) {
	    print color 'green' if $print_color;
	    print "$name_root TriAnnotation complete\n";
	    print color 'reset' if $print_color;
	} 
	else {
	    $pass_audit = 0;
	    print color 'red' if $print_color;
	    print "$name_root TriAnnotation incomplete\n";
	    print color 'reset' if $print_color;
	}

    }

    # Print spacer
    print LOG "\n" if $logfile;
    print "\n";


    # Reset the pass_audit
    if ($pass_audit) {
	push @pass_files, $name_root;
	$num_pass++;
    }
    $pass_audit = 1;

} # End of for each file in the input folder


#-----------------------------+
# PRINT SUMMARY OF AUDIT      |
#-----------------------------+
# Try to sort the output array
#sort {$b cmp $a}(@pass_files);

print "=========================================\n";
print "  SUMMARY\n";
print "=========================================\n\n";

print "Total seqs: $num_files\n";
print "Passed audit: $num_pass\n";
for my $ind_pass_file (@pass_files) {
    print "\t$ind_pass_file\n";
}

close LOG if $logfile;

exit;

#-----------------------------------------------------------+
# SUBFUNCTIONS                                              |
#-----------------------------------------------------------+

sub audit_gene {
# Audit locally produced gene annotation results

}

sub audit_rm {

    #-----------------------------+
    # AUDIT REPEAT MASKER OUTPUT  |
    #-----------------------------+
    my ($outdir, $name_root) = @_;
    my $msg;
    my $rm_fail = 0;

    my $rm_dir = $outdir.$name_root."/rm/";

    # Expected rm files
    my @rm_files = ( $name_root."_TREP9_EL.gff",
		     $name_root."_ALLDB.rm.gff", 
		     );
    

    if (-e $rm_dir) {

	for my $ind_rm_file (@rm_files) {
	    my $rm_file_path = $rm_dir.$ind_rm_file;
	    if (-e $rm_file_path) {
		#print "\tokay\n";
	    }
	    else {
		print LOG "MISSING: $name_root: $ind_rm_file\n" if $logfile;
		print "MISSING: $name_root: $ind_rm_file\n" if $verbose;
		$rm_fail = 1;
	    }

	}

    }
    else {
	$rm_fail = 1;
	print "The RepeatMasker output dir does not exist at:\n";
	print "    $rm_dir\n";
    }

    #-----------------------------+
    # RETURN TRUE IF PASSED AUDIT |
    #-----------------------------+
    if ($rm_fail) {
	return 0;
    }
    else {
	return 1;
    }


}

sub audit_ta {

    #-----------------------------+
    # AUDIT TRIANNOTATION OUTPUT  |
    #-----------------------------+

    my ($outdir, $name_root) = @_;
    my $msg;
    my $ta_fail = 0;
    
    my $ta_dir = $outdir.$name_root."/ta/";
    my $gff_dir = $outdir.$name_root."/gff/";
    my $game_dir = $outdir.$name_root."/game/";

    # Expected ta files
    my @ta_files = (# Original output from TriAnnotation
		    "1trf.gff",
		    "2eugOs.gff",
		    "2fGh.gff",
		    "2gID.gff",
		    "2gmHv.gff",
		    "2gmOs.gff",
		    "2gmTa.gff",
		    "2gmZm.gff",
		    # Apollo formatted output
		    "1trf.ap.gff",
		    "2eugOs.ap.gff",
		    "2fGh.ap.gff",
		    "2gID.ap.gff",
		    "2gmHv.ap.gff",
		    "2gmOs.ap.gff",
		    "2gmTa.ap.gff",
		    "2gmZm.ap.gff",
		    );
    
    if (-e $ta_dir) {
	# Check that expected files are present
	for my $ind_ta_file (@ta_files) {
	    my $ta_file_path = $ta_dir.$ind_ta_file;
	    if (-e $ta_file_path) {
		# If this is a ap.gff file then copy to the gff dir
		if ($ind_ta_file =~ m/.*\.ap\.gff$/ ) {
		    print "\tMaking copy of $ind_ta_file\n" if $verbose;
		    my $gff_file_path = $gff_dir.$ind_ta_file;
		    copy ($ta_file_path, $gff_file_path)
		}
		#print LOG "";
	    }
	    else {
		$ta_fail = 1;
		print LOG "$name_root MISSING: $ind_ta_file\n" if $logfile;
		print "$name_root: MISSING $ind_ta_file\n" if $verbose;
	    }
	}



    }
    else {
	$ta_fail = 1;
	print "The TriAnnotation output dir does not exist at:\n" if $verbose;
	print "    $ta_dir\n" if $verbose;
	
	print LOG "The TA output dir does not exist at:\n" if $logfile;
	print LOG "    $ta_dir\n" if $logfile;

    }

    #-----------------------------+
    # RETURN TRUE IF PASSED AUDIT |
    #-----------------------------+
    if ($ta_fail) {
	return 0;
    }
    else {
	return 1;
    }

}

sub audit_blast {
    
    #-----------------------------+
    # AUDIT BLAST OUTPUT          |
    #-----------------------------+
    
    my ($outdir, $name_root, $do_deep_audit) = @_;
    my $msg;
    my $blast_fail = 0;

    # Expected location of blast output
    my $blast_dir = $outdir.$name_root."/blast/";
    my $no_hit_dir = $blast_dir."/no_hits/";
    
    if ( (-e $blast_dir) & ($do_deep_audit) ) {
	# If the output directory exists and this is a deep 
	# audit request
	#////////////////////////////
	# BEGIN CHECK EACH FILE
	#\\\\\\\\\\\\\\\\\\\\\\\\\\\\
	opendir( BLASTDIR, $blast_dir ) || 
	    die "Can't open directory:\n$blast_dir"; 
	my @blast_files = grep /\.blx$|\.bln$/, readdir BLASTDIR;
	close BLASTDIR;

	#-----------------------------+
	# MAKE DIR FOR NO HITS        |
	#-----------------------------+
	my $no_hit_dir = $blast_dir."/no_hits/";
	unless (-e $no_hit_dir) {
	    mkdir $no_hit_dir;
	}
	
	#-----------------------------+
	# CHECK EACH BLAST FILE       |
	#-----------------------------+
	for my $ind_blast (@blast_files) {
	    
	    my $blast_file_path = $blast_dir.$ind_blast;
	    
	    #print "Processing: $ind_blast" if $verbose;
	    #print "\tProcessing: $blast_file_path\n";
	    
	    # The blast object
	    $msg = "Could not open BLAST FILE:\n$blast_file_path\n";
	    my $blobj = new Bio::SearchIO ( '-format' => 'blast',
					    '-file'   => $blast_file_path) 
		|| die "$msg\n";
	    
	    # The blast result
	    while ( my $blres = $blobj->next_result() ) {		

		if ($verbose) {
		    print "$ind_blast\tHITS:".$blres->num_hits()."\n";
		}
		
		if ($logfile) {
		    print LOG "$ind_blast\tHITS:".$blres->num_hits()."\n";
		}
		
		# If not hits, move the output to the no_hits dir
		if ($blres->num_hits() == 0) {
		    my $new_file_path = $blast_dir."/no_hits/".$ind_blast;
		    $msg = "Could not move from\n".
			"\t$blast_file_path\n".
			"To:".
			"\t$new_file_path\n";
		    move ($blast_file_path, $new_file_path) ||
			die "$msg\n";
		}
		    
	    }
	} # End of for each blast output file

	#////////////////////////////
	# END OF CHECK EACH FILE
	#\\\\\\\\\\\\\\\\\\\\\\\\\\\\

    } 
    elsif (-e $blast_dir) {
	# The output dir exists but this is not a deep audit
	# Maybe just count the number of blast output files and 
	# pass fail depending on number of files or file list
    }	
    else {
	$blast_fail = 1;
	print "The Blast output dir does not exist at:\n" if $verbose;
	print "    $blast_dir\n" if $verbose;
    } # End of if the blast_dir exists


    #-----------------------------+
    # RETURN TRUE IF PASSED AUDIT |
    #-----------------------------+
    if ($ta_fail) {
	return 0;
    }
    else {
	return 1;
    }

}

sub print_help {

    # Print requested help or exit.
    # Options are to just print the full 
    my ($opt) = @_;
    
    my $usage = "USAGE:\n".
	"  batch_audit.pl -i DirToProcess -o OutDir";
    my $args = "REQUIRED ARGUMENTS:\n".
	"  --indir        # Path to the directory containing the sequences\n".
	"                 # to process. The files must have one of the\n".
	"                 # following file extensions:\n".
	"                 # [fasta|fa]\n".
	"  --outdir       # Path to the output directory\n".
	"\n".
	"OPTIONS:\n".
	"  --color        # Print output in fancy ANSI color\n".
	"  --deep         # Does a deep audit\n".
	"                 #   -Moves empty BLAST output to subdir\n".
	"  --logfile      # Path to file to use for logfile\n".
	"  --version      # Show the program version\n".     
	"  --usage        # Show program usage\n".
	"  --help         # Show this help message\n".
	"  --man          # Open full program manual\n".
	"  --test         # Run the program in test mode\n".
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

STARTED: 07/30/2007

UPDATED: 08/02/2007

=cut

#-------------------------------------------------------+
# HISTORY                                               |
#-------------------------------------------------------+
#
# 07/30/2007
# - Program started
#
# 07/31/2007
# - Changed BLAST audit to a subfunction
# - added audit_rm subfunction
# - added audit_ta subfunction
# - added booleans to keep track of pass/fail of 
#   individual audit jobs
#
# 08/02/2007
# - Added ability to print output in color
# - Added counter for seqs that pass the audit
# - Added array to hold the list of passed seqs
#   and these get printed at the end
# - Added sort of the fasta_files array
