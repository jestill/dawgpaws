#!/usr/bin/perl -w
#-----------------------------------------------------------+
#                                                           |
# seq_oligocount.pl - Count oligos of size k
#                                                           |
#-----------------------------------------------------------+
#                                                           |
#  AUTHOR: James C. Estill                                  |
# CONTACT: JamesEstill_@_gmail.com                          |
# STARTED: 10/11/2007
# UPDATED: 10/12/2007
#                                                           |
# DESCRIPTION:                                              |
#  Short Program Description                                |
#                                                           |
# USAGE:                                                    |
#  ShortFasta Infile.fasta Outfile.fasta                    |
#                                                           |
# VERSION: $Rev: 348 $                                            |
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
use Bio::SeqIO;                # Seq IO used to to work with the input file
use strict;
use Getopt::Long;

#-----------------------------+
# PROGRAM VARIABLES           |
#-----------------------------+
my ($VERSION) = q$Rev$ =~ /(\d+)/;

#-----------------------------+
# VARIABLE SCOPE              |
#-----------------------------+
my $indir;                     # Directory containing fasta files to process
my $outdir;                    # Base directory for output
my $indexdir;                  # Directory containing the index files
my $file_config;               # Path to the configuration file
my $index;                     # Path to the sequence index file
#my $kmer_len = 20;             # Oligomer length, default is 20
#my $seq_name;                  # Sequence name used in the output

# Booleans
my $quiet = 0;
my $verbose = 0;
my $show_help = 0;
my $show_usage = 0;
my $show_man = 0;
my $show_version = 0;
my $do_gff = 0;               # Boolean to create the gff file
my $thresh = '50';

# ARRAYS
my @params;                    # Parameters array for running oligocounts
                               # Index starts at zero
#-----------------------------+
# COMMAND LINE OPTIONS        |
#-----------------------------+
my $ok = GetOptions(# REQUIRED OPTIONS
		    #"i|infile=s"  => \$infile,
		    "i|indir=s"   => \$indir,
		    "d|db=s"      => \$indexdir,
                    "o|outdir=s"  => \$outdir,
		    "c|config=s"  => \$file_config,
		    # ADDITIONAL OPTIONS
		    "t|thresh=s"  => \$thresh,
		    # BOOLEAN OPTIONS
		    "q|quiet"     => \$quiet,
		    "verbose"     => \$verbose,
		    "gff"         => \$do_gff,
		    # ADDITIONAL INFORMATION
		    "usage"       => \$show_usage,
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

if ($show_version) {
    print "\n$0:\nVersion: $VERSION\n\n";
    exit;
}

if ($show_man) {
    # User perldoc to generate the man documentation.
    system("perldoc $0");
    exit($ok ? 0 : 2);
}

# Throw error if required options not present
if ( (!$indir) || (!$outdir) || (!$indexdir) ) {

    print "\a";
    print "\n";
    print "ERROR: Input file path required\n" if (!$indir);
    print "ERROR: Output directory required\n" if (!$outdir);
    print "ERROR: Index dir path required\n" if (!$indexdir);
    print_help("");

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

unless ($indexdir =~ /\/$/ ) {
    $outdir = $indexdir."/";
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

print STDERR "NUMBER OF FILES TO PROCESS: $num_files\n" if $verbose;

#-----------------------------+
# CREATE THE OUT DIR          |
# IF IT DOES NOT EXIST        |
#-----------------------------+
unless (-e $outdir) {
    print "Creating output dir ...\n" unless $quiet;
    mkdir $outdir ||
	die "ERROR: Could not create the output directory:\n$outdir";
}

#-----------------------------+
# GET INFO FROM CONFIG FILE   |
#-----------------------------+
# Load to the 2d dbs array
open (CONFIG, $file_config) ||
    die "Can't open the config file:\n$file_config";
open (CONFIG, $file_config) ||
    die "Can't open the config file:\n$file_config";
$i = 0;
my $line_num = 0;
while (<CONFIG>) {
    $line_num ++;
    unless (m/^\#/) {
	chomp;
	my @tmpary = split (/\t/);
	my $count_tmp = @tmpary;
	
	if ($count_tmp == 2) {
	 
	    $params[$i][0] = $tmpary[0];
	    $params[$i][0] = $tmpary[0];
	    $i++;
	    
	}
	else {
	    print "ERROR: Config file line number $line_num\n";
	    print "       Only $line_num variables were found\n" 
	}
	
    }
}
close CONFIG;

# Number of parameter sets specified in the config file
my $num_par_sets = $i;

my $num_proc_total = $num_files * $num_par_sets;

print STDERR "$num_proc_total find_ltr runs to process\n";

#-----------------------------+
# OLIGO COUNT FOR EACH FILE   |
# FOR EACH CONFIG SET         |
#-----------------------------+ 
for my $ind_file (@fasta_files) {
    
    $file_num++;

    #-----------------------------+
    # Get the root name of the    |
    # file to mask                |
    #-----------------------------+
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
    # Create parent dir if it     |
    # does not already exist      |
    #-----------------------------+
    my $dir_parent = $outdir.$name_root."/";
    unless (-e $dir_parent) {
	print "creating dir: $dir_parent\n";
	mkdir $dir_parent ||
	    die "Could not creat the output dir:\n$dir_parent\n";
    }

    #-----------------------------+
    # Create the dir to hold all  |
    # the oligocount outupt       |
    #-----------------------------+
    my $dir_oligocount = $dir_parent."oligocount/";
    unless (-e $dir_oligocount ) {
	print "Creating output dir\n: $dir_oligocount\n" if $verbose;
	mkdir $dir_oligocount ||
	    die "Could not create the output directory:\n$dir_oligocount";
    }

    #-----------------------------+
    # Create the dir to hold all  |
    # the oligocount outupt       |
    #-----------------------------+
    my  $dir_gff = $dir_parent."gff/";
    unless (-e $dir_gff ) {
	print "Creating output dir\n: $dir_gff\n" if $verbose;
	mkdir $dir_gff ||
	    die "Could not create the output directory:\n$dir_gff";
    }


    for ($i=0; $i<$num_par_sets; $i++) {

	my $kmer_len = $params[$i][0];
	my $db_name  = $params[$i][0];
	
	#-----------------------------+
	# RUN KMER COUNT FOR INFILE   |
	#-----------------------------+
	seq_kmer_count ($infile,$dir_parent,$index, $kmer_len, $name_root,
			$do_gff);

    }

}

exit 1;

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

sub seq_kmer_count {

    my ($fasta_in, $outdir, $vmatch_index, $k, $seq_name,
	$do_gff) = @_;
    
    # Array of threshold values
    #my @thresh = ("200");
    my $thresh = 50;

    my $in_seq_num = 0;

    # TODO:
    # Add option to do something here when infile not encountered
    my $inseq = Bio::SeqIO->new( -file => "<$fasta_in",
				 -format => 'fasta');

    
    #-----------------------------+
    # CREATE OUTPUT DIR IF NEEDED |
    #-----------------------------+
    


    # Counts of the number of occurences of each oligo
    my @counts = ();    
    my @vdat;   # temp array to hold the vmatch data split by tab
    my $pos = 0;

    my $start_pos;
    my $end_pos;
    my $i;         # Array index value

    #-----------------------------+
    # FILE PATHS                  |
    #-----------------------------+
    # These will need to be modified later to more useful names
    my $temp_fasta = $outdir."split_seq_".$k."_mer.fasta";
    my $vmatch_out = $outdir."vmatch_out.txt";
    my $gff_count_out = $outdir."vmatch_out.gff";
    my $ocount_out = $outdir."fast_count.txt";

    #-----------------------------+
    # WRITE FIRST LINE OF OCOUNT  |
    # FILE                        |
    #-----------------------------+
    open (OCOUNT,">$gff_count_out") ||
	die "Can not open oligo count file for output:\n$ocount_out\n";
    print OCOUNT ">".$seq_name."_".$k."mer\n";
    my $ocount_col = 1;             # Set ocount col to one

    while (my $seq = $inseq->next_seq) {

	$in_seq_num++;
	if ($in_seq_num == 2) {
	    print "\a";
	    die "Input file should be a single sequence record\n";
	}

	# Calculate base cooridate data
	my $seq_len = $seq->length();
	my $max_start = $seq->length() - $k;
	
	# Print some summary data if running in verbose mode
	print STDERR "\n==============================\n" if $verbose;
	print STDERR "SEQ LEN: $seq_len\n" if $verbose;
	print STDERR "MAX START: $max_start\n" if $verbose;
	print STDERR "==============================\n" if $verbose;
	
	#-----------------------------------------------------------+
	# CREATE FASTA FILE OF ALL K LENGTH OLIGOS                  |
	# IN THE INPUT SEQUENCE                                     |
	#-----------------------------------------------------------+
	print STDERR "Creating oligo fasta file\n" if $verbose;
	open (FASTAOUT, ">$temp_fasta") ||
	    die "Can not open temp fasta file:\n $temp_fasta\n";

	for ($i=0; $i<=$max_start; $i++) {

	    $start_pos = $i + 1;
	    $end_pos = $start_pos + $k - 1;

	    my $oligo = $seq->subseq($start_pos, $end_pos);

	    # Set counts array to zero
	    $counts[$i] = 0;
	    
	    print FASTAOUT ">$start_pos\n";
	    print FASTAOUT "$oligo\n";
	    
	}

	close (FASTAOUT);

	# QUERY OLIGO FASTA FILE AGAINST VMATCH INDEX
	my $vmatch_cmd = "vmatch -q $temp_fasta -complete".
	    " $vmatch_index > $vmatch_out";
#	# OUTPUT TO STDOUT FOR DEBUG
#	my $vmatch_cmd = "vmatch -q $temp_fasta -complete".
#	    " $vmatch_index";
	print STDERR "\nVmatch cmd:\n$vmatch_cmd\n" if $verbose;
	system ($vmatch_cmd);

	# PARSE VMATCH OUTPUT FILE
	# increment
	unless (-e $vmatch_out) {
	    print STDERR "Can not file the expected vmatch output file\n";
	    die;
	}

	# PARSE THE VMATCH OUTPUT FILE AND INCREMENT
	# THE COUNTS ARRAY
	open (VMATCH, "<$vmatch_out") ||
	    die "Can not open vmatch output file:\n$vmatch_out\n";

	# Count oligos hits in index file
	print STDERR "\nCounting oligos ...\n" if $verbose;
	while (<VMATCH>) {
	    # ignore comment lines
	    unless (m/^\#/) {
		chomp;
		@vdat = split;
		
		my $len_vdat = @vdat;
		# print the following for debu
		#print "VDAT LEN: $len_vdat\n";
		#print $vdat[5]."\n";
		
		# Get the seq file in the new sub se
		# Counts index starts at zero
		# It may be possible to increment without needing
		# to add one
		#$counts[$vdat[5]]++;
		$counts[$vdat[5]] = $counts[$vdat[5]] + 1;

	    }
	}
	close (VMATCH);

	#///////////////////////////////////
	# NEED SEGMENTATION STEP HERE
	# THIS WILL JOIN INDIVIDUAL PIPS
	# FOR A GIVEN THRESHOLD COVERAGE
	# This could also be done with a 
	# separate script operating on the
	# gff output from this program.
	#\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

	# Output oligo counts to gff file
	
	# This is the pip segment file, a single pip with depth coverage data
	# for eack oligo segment
	if ($do_gff) {
	    open (GFFCOUNT, ">$gff_count_out") ||
		die "Can not open gff out file:\n$gff_count_out\n";
	}


	print STDERR "\nCreating output files ...\n" if $verbose;

	for ($i=0; $i<=$max_start; $i++) {

	    $start_pos = $i + 1;
	    $end_pos = $start_pos + $k - 1;
	    
	    #-----------------------------+
	    # PRINT OCOUNT OUT            |
	    #-----------------------------+
	    print OCOUNT $counts[$i];
	    $ocount_col++;
	    if ($ocount_col > 80) {
		print OCOUNT "\n";
		$ocount_col = 1;
	    }
	    else {
		print OCOUNT " ";
	    }
	    
	    

	    #-----------------------------+
	    # GFF PIP OUTPUT FILE         |
	    #-----------------------------+
	    # The count will be placed in the score position
	    my $seq_str = $seq->subseq($start_pos, $end_pos);

	    if ($do_gff) {
		print GFFCOUNT "$seq_name\t".  # Ref sequence
		    "vmatch\t".                # Source
		    "wheat_count\t".           # Type
		    "$start_pos\t".            # Start
		    "$end_pos\t".              # End
		    $counts[$i]."\t".          # Score
		    ".\t".                     # Strand
		    ".\t".                     # Phase
		    "Vmatch ".$k."mer\n";    # Group
	    }

	    
	    #-----------------------------+
	    # PRINT OLIGOS OVER THRESHOLD |
	    # TO STDOUT                   |
	    #-----------------------------+
	    # PRINT OUT SEQS EXCEEDING THRESHOLD VALUES TO STDOUT
	    if ($counts[$i] > $thresh) {
		my $thresh_seq = $seq->subseq($start_pos, $end_pos);
		print STDOUT "$start_pos\t".$counts[$i]."\t".
		    "$thresh_seq\n";
	    }

	}
	
	close (GFFCOUNT) if ($do_gff);
	close (OCOUNT);

    } # End of while seq object

    # May want to make a single fasta file of all oligos for the
    # qry fasta_sequence and parse the results from vmatch from that
    
}

1;
__END__

# OLD print_help subfunction
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

seq_oligocount.pl - Count oligos from an input sequence

=head1 VERSION

This documentation refers to program version $Rev: 348 $

=head1 SYNOPSIS

=head2 Usage

    seq_oligocount.pl -i InFile -o OutDir -db index.fasta
                      -n SeqName -k 20

=head2 Required Argumenta

    -i,--infile   # Path to the input fasta file
    -d,--db       # Path to the mkvtree index file
    -o,--outdir   # Path to the base output directory
    -n,--name     # Name to assign to the sequence file
    -k,--kmer     # Oligomer query length

=head1 DESCRIPTION

The seq_oligocount program will take a query sequence, break it
into subsequences of size k and query it against an persistent index
created by the mkvtree program. It produces a GFF output file describing
the number of copies of every oligomer in the query sequence in the 
subject index database.

The variables I am considering incorporating

    --thresh     # Comma delimited array of integers to describe
                 # Currently using a single value, default is 50
                 # Currently just determins if this value is printed
                 # to STDOUT
    --center     # Location to center point value of summary info as
                 # a point, 
                 #  start of oligomer
                 #  center of oligomer
                 #  end of oligomer
                 # [left|right|center]
    --gff        # BOOLEAN. Don't always create a gff output


=head1 REQUIRED ARGUMENTS

=over 2

=item -i,--infile

Path of the input file.

=item -o,--outdir

Path of the output file.

=item -n,--name

Name to assign to the sequence file

=item -d,--db

Path to the fasta file that was indexed with the mkvtree program.

=back

=head1 OPTIONS

=over 2

=item -k,--kmer

Length of the kmer to index. The default value of this variable is 20.

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

=item ERROR: Could not create the output directory

The output directory could not be created at the path you specified. 
This could be do to the fact that the directory that you are trying
to place your base directory in does not exist, or because you do not
have write permission to the directory you want to place your file in.

=back

=head1 CONFIGURATION AND ENVIRONMENT

An external configuration file is not required for this program, and
it does not make use of any variables set in the users environment.

=head1 DEPENDENCIES

=head2 Required Software

=over 2

=item Vmatch

This program requires the Vmatch package of programs.
http://www.vmatch.de . This software is availabe at no cost for
noncommercial academic use. See program web page for details.

=back

=head2 Required Perl Modules

=over

=item * File::Copy

This module is required to copy the BLAST results.

=item * Getopt::Long

This module is required to accept options at the command line.

=item * Bio:SeqIO

The Bio:SeqIO module is a component of the BioPerl package

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

=item * Limit to query sequence file

Since this program will generate a large file of your original sequence
broken into oligomers, you are somewhat limited to the size of query
sequence that is feasible to use this program with. 

=item * Output limited to GFF file

The output for this program is currently limited to GFF file of 
points with the copy number of each oligomer designeated.

=back

=head1 SEE ALSO

The seq_oligocount.pl program is part of the DAWG-PAWS package of genome
annotation programs. See the DAWG-PAWS web page 
( http://dawgpaws.sourceforge.net/ )
or the Sourceforge project page 
( http://sourceforge.net/projects/dawgpaws ) 
for additional information about this package.

=head1 LICENSE

GNU General Public License, Version 3

L<http://www.gnu.org/licenses/gpl.html>

=head1 AUTHOR

James C. Estill E<lt>JamesEstill at gmail.comE<gt>

=head1 HISTORY

STARTED: 12/21/2007

UPDATED: 12/21/2007

VERSION: $Rev$

=cut

#-----------------------------------------------------------+
# HISTORY                                                   |
#-----------------------------------------------------------+
#
# 12/21/2007
# - Started batch oligo count from seq_oligocount.pl
# - Modified to accept a directory of fasta files as input
#   instead of a single fasta file
# - Takess directory paths as input:
#    --indir      # Path to indir with fasta files
#    --outdir     # Path to base outdir
#    --indexdir   # Path to dir with index databases
# - Now use a config file to set k, etc
