#!/usr/bin/perl -w
#-----------------------------------------------------------+
#                                                           |
# cnv_mgescan2gff.pl - Converte MGEScan output to GFF       |
#                                                           |
#-----------------------------------------------------------+
#                                                           |
#  AUTHOR: James C. Estill                                  |
# CONTACT: JamesEstill_@_gmail.com                          |
# STARTED: 01/24/2012                                       |
# UPDATED: 01/24/2012                                       |
#                                                           |
# DESCRIPTION:                                              |
#  Convert MGEScan fasta output to a GFF format output file | 
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
# Term in LINE_element

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
use Bio::SeqIO;                # Read and write seq files in different formats

#-----------------------------+
# PROGRAM VARIABLES           |
#-----------------------------+
my ($VERSION) = q$Rev$ =~ /(\d+)/;
# Get GFF version from environment, GFF2 is DEFAULT
my $gff_ver = uc($ENV{DP_GFF}) || "GFF2";

#-----------------------------+
# VARIABLE SCOPE              |
#-----------------------------+
my $indir;
my $outfile;
my $delim = "_";
my $param;
my $program = "MGEScan-nonLTR";

# http://www.sequenceontology.org/browser/current_cvs/term/SO:0000186
my $parent_feature = "LINE_element";     # The SO complient term
# Currently SO does not support additional resolution annotation of LINE elements

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
		    "i|infile=s"  => \$indir,
                    "o|outfile=s" => \$outfile,
		    # ADDITIONAL OPTIONS
		    "p|param=s"   => \$param,
		    "program=s"   => \$program,
		    "delim=s"     => \$delim,
		    "gff-ver=s"   => \$gff_ver,
		    "q|quiet"     => \$quiet,
		    "verbose"     => \$verbose,
		    # ADDITIONAL INFORMATION
		    "usage"       => \$show_usage,
		    "test"        => \$do_test,
		    "version"     => \$show_version,
		    "man"         => \$show_man,
		    "h|help"      => \$show_help,);


#-----------------------------+
# PRINT REQUESTED HELP        |
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
    print "\ncnv_mgescan2gff.pl:\n".
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
# Get the FASTA files from the|
# directory provided by the   |
# var $indir                  |
#-----------------------------+
# Expect input as dir but will accept single files and
# load this to the fasta file array
my @fasta_files;
if (-d $indir) {
    opendir( DIR, $indir ) || 
	die "Can't open directory:\n$indir"; 
    @fasta_files = grep /\.fasta$|\.fa$/, readdir DIR ;
    closedir( DIR );
    
    # add slash
    unless ($indir =~ /\/$/ ) {
	$indir = $indir."/";
    }

}
elsif (-f $indir) {
    print STDERR "Expecting this to be a file.\n" if $verbose;
    push (@fasta_files, $indir);
}


#-----------------------------+
# SHOW ERROR IF NO FILES      |
# WERE FOUND IN THE INPUT DIR |
#-----------------------------+
my $count_files = @fasta_files;
if ($count_files == 0) {
    print STDERR "\a";
    print STDERR "\nERROR: No fasta files were found in the input directory\n".
	"$indir\n".
	"Fasta files must have the fasta or fa extension.\n\n";
    exit;
}

print STDERR "NUMBER OF FILES TO PROCESS: $count_files\n" if $verbose;


#-----------------------------+
# MAIN PROGRAM BODY           |
#-----------------------------+



# OUTPUT FILE
if ($outfile) {
    open(GFFOUT, $outfile) ||
	die "Can not open output file $outfile\n";
}
else {
    open (GFFOUT, ">&STDOUT") ||
	die "Can not print to STDOUT\n";
}
# Print GFF 3 header
if ($gff_ver =~ "GFF3") {
    print GFFOUT "##gff-version 3\n";
}


if ($param) {
    $program = $program.$delim.$param;
}

my $model_num = 0;
for my $ind_seq_file (@fasta_files) {

    # Get input file path by adding dir, with exception
    # for situations where indir was a file
    my $infile;
    if (-d $indir) {
	$infile = $indir.$ind_seq_file;
    }
    elsif (-f $indir) {
	$infile = $indir;
    }


    
    my $seq_in;
    # OPEN THE INPUT FILE
    my $infile_format = "fasta";
    if ($infile) {
	$seq_in = Bio::SeqIO->new(-file   => "<$infile",
				  -format => $infile_format) 
	    || "Can not connect to input file $infile\n";
    }

  
    while (my $inseq = $seq_in->next_seq) {

	# Increment model number, this should increment across
	# all the fasta files in the input dir and profide a 
	# unique identifier.
	$model_num++;
	
	my $seq_primary_id = $inseq->primary_id();
	my $seq_length = $inseq->length();


	#-----------------------------+
	# GET FEATURE LOCATIONS       |
	#-----------------------------+
	# Split the id in the input field to get the start point
	# in the genome
	my @header_parts = split (/\_/, $seq_primary_id);
	my $num_parts = @header_parts;
	# Get the start point from the split array
	# this should be the last item in the array
	my $feat_start = $header_parts[$num_parts - 1];
	my $feat_end = $feat_start + $seq_length;

	#-----------------------------+
	# GET SCAFFOLD ID             |
	#-----------------------------+
	# This will need to be extracted from the fasta header
	# If the input contains extra header information
	if ($seq_primary_id =~ m/(.*)\|(.*)/) {
	    $seq_primary_id = $2;
	}
	if ($seq_primary_id =~ m/(.*)\.fasta(.*)/) {
	    $seq_primary_id = $1;
	}
	if ($seq_primary_id =~ m/(.*)\.fa(.*)/) {
	    $seq_primary_id = $1;
	}

	#-----------------------------+
        # SET GFF ATTRIBUTE           |
	#-----------------------------+
	my $attribute;
	if ($gff_ver =~ "GFF3") {
	    $attribute = "ID=".$seq_primary_id."_".
		$program."_".$model_num.";".
		"Name=".$program."_".$model_num;
	}
	elsif ($gff_ver =~ "GFF2") {
	    $attribute = $program."_".$model_num;
	}
	else {
	    print STDERR "ERROR: GFF Version not recognized".
		$gff_ver."\n";
	    die;
	}

	# PRINT OUT TO STDERR FILE FOR TEST
	if ($verbose) {
	    print STDERR $seq_primary_id."\t";
	    print STDERR $seq_length."\t";
	    print STDERR $feat_start."\t";
	    print STDERR $feat_end."\t";
	    print STDERR "\n";
	}


	# Results are currently not stranded in the fasta output
	# but could perhaps be made stranded by mapping the 
	# MGEScan fasta file to the sequence contig used in 
	# the annotation.

	#-----------------------------+
	# PRINT GFF OUTPUT            |
	#-----------------------------+
	print GFFOUT $seq_primary_id."\t".  # 1.
	    $program."\t".                  # 2. program
	    $parent_feature."\t".           # 3. feature type
	    $feat_start."\t".               # 4. feature start
	    $feat_end."\t".                 # 5. feature end
	    $seq_length."\t".               # 6. score (length)
	    ".\t".                          # 7. strand
            ".\t".                          # 8. frame
            $attribute."\n"                 # 9. attribute


    } # End of next seq in sequence record

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

cnv_mgescan2gff.pl - Converte MGEScan output to GFF foramt.

=head1 VERSION

This documentation refers to program version 0.1

=head1 SYNOPSIS

=head2 Usage

    cnv_mgescan2gff.pl -i infile.fasta -o outfile.gff

=head2 Required Arguments

    --infile        # Path to the input file or dir
    --outfie        # Path to the output file

=head1 DESCRIPTION

Takes a fasta file dervied from MGEScan or a dir containing the fasta files
and converts these to a GFF3 file

=head1 REQUIRED ARGUMENTS

=over 2

=item -i,--indir

Path of the input file of MGEScan results or the dir of fasta files containing
all MGEScan results of interest.

=item -o,--outfile

Path of the output file. This will be a GFF format file.

=back

=head1 OPTIONS

=over 2

=item --program

The program name to use to identify the program in the GFF3 output. By default
the program name is set to MGEScan-nonLTR. This will be reportd in
column 2 of the GFF output file.

=item --param

The parameter name to use to tag the GFF output. This will be appended
to column 2 of the GFF output file. This allows for segregating different
runs of the MGEScan program that used different parameters for running
the program.

=item --delim

The delimiting character to use to separate the program name from the
parameter name in column 2 of the GFF output file. Be default this is
the '_' character. For example given the program MGEScan and the 
param name default, the output will be MGEScan_default.

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

=head1 EXAMPLES

The following are examples of how to use this script

=head2 Typical Use

  cnv_mgescan2gff.pl -i MGEScanResultFile.fasta -o MGEScanresult.gff 
                     --gff-ver GFF3 --param default

=head1 DIAGNOSTICS

=over 2

=item * Expecting input from STDIN

If you see this message, it may indicate that you did not properly specify
the input sequence with -i or --infile flag. 

=back

=head1 CONFIGURATION AND ENVIRONMENT

Names and locations of config files
environmental variables
or properties that can be set.

=head1 DEPENDENCIES

Other modules or software that the program is dependent on.

=head2 MGEScan-nonLTR

This program requires output from the MGEScan-nonLTR program to process.
This program is available from the authors 
http://darwin.informatics.indiana.edu/cgi-bin/evolution/nonltr/nonltr.pl

=head1 BUGS AND LIMITATIONS

Any known bugs and limitations will be listed here.

=head1 REFERENCE

Use of the DAWGPAWS suite of tools should cite:

JC Estill and JL Bennetzen. 2009. "The DAWGPAWS pipeline for the annotation 
of genes and transposable elements in plant genomes." Plant Methods. 5:8.

=head1 LICENSE

GNU General Public License, Version 3

L<http://www.gnu.org/licenses/gpl.html>

=head1 AUTHOR

James C. Estill E<lt>JamesEstill at gmail.comE<gt>

=head1 HISTORY

STARTED: 01/24/2012

UPDATED: 01/24/2012

VERSION: $Rev$

=cut

#-----------------------------------------------------------+
# HISTORY                                                   |
#-----------------------------------------------------------+
#
