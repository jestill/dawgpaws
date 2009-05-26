#!/usr/bin/perl -w
#-----------------------------------------------------------+
#                                                           |
# tsd_check.pl - Target Site Duplication Checker            |
#                                                           |
#-----------------------------------------------------------+
#  AUTHOR: James C. Estill                                  |
# CONTACT: jestill_at_sourceforge.net                       |
# STARTED: 06/08/2006                                       |
# UPDATED: 06/08/2006                                       |
#                                                           |
# DESCRIPTION:                                              |
#  Extract target site duplication sequences from a genomic |
#  sequence given the start and stop positions of an LTR    |
#  retrotransposon.                                         |
#                                                           |
# INPUT VARIABLES:                                          |
#   [f] FASTA file containing ALL contigs of interest       |
#       This can be a single chrom, multiple chrom,         |
#       multiple BACs etc.                                  | 
#   [a] Annotation file in GFF Format                       |
#   [o] Output file path                                    |
#   [l] Length of LTR seq to return                         |
#       Default value is 3                                  |
#   [t] Length of putative TDS to return                    |
#       Default value is 5                                  |
#   [q] Quiet mode. Do not print program status to screen   | 
#                                                           |
# REQUIRES:                                                 |
#   - bioperl                                               |
#     Makes use of the BioPERL seqIO function               |
#                                                           |
# USAGE:                                                    |
#  tsd_sel.pl -f SeqFile.fasta -a LTRAnn.gff -o OutFile.txt |

# Can use multiple tsds by using multiple t in input
# ie. -t 4 -t 5 -t 6
#-----------------------------------------------------------+

#-----------------------------+
# INCLUDES                    |
#-----------------------------+
use Bio::SeqIO;                # Allows for fasta sequence input/output
use Getopt::Long;              # Get input from the command line

#-----------------------------+
# GLOBAL VARIABLES            |
#-----------------------------+
# These variables need to be local to the program
# but need to be globally avaiable to the 
# subfunctions. The user will not need to 
# modify these variables.
my @StartGeneExons;            # Array to hold information for start exon
my @Annot;                     # Array to hold all annotation info for BAC
my @StartGeneList;             # Array to hold the list of gene names
                               # to use as the start genes.
my $StartPos;                  # Random start position on the BAC for
                               # one end of a simulated clone
my $SearchPos;                 # Search position to look for gene on other
                               # end of the fake clone.
my $CumSum = "0";              # Var to hold cumulative sum
my $ans;                       # Var to hold answers to questions
my $NumExonStartGene;          # Number of exons in the gene model
my $UsableLen;                 # Usable length of the exon
my $SeqId;                     # Sequence ID of contig sequence record to use
my $RandCumSum;                # Random cumulative sum
my $SeqLength;                 # Length of the sequence file
my $NumTSDMatch;               # Number of matches in TSD
my $MatchString;               # String to visualize match
my $MaxNum;                    # Max position of LTR elements when
                               # stored in the array

# Vars with default values
#my @TSDLens = ("4");           # Array of TSD Lengths to try
#my $TSDLen = "5";              # Length of the target site duplication
my $LTRLen = "3";              # Length of LTR shown in verbose output
my $gff_out;                   # Output in the gff2 format
my $gff3_out;                  # Output in the gff3 format
                               # For now this assumes input in gff2 format

# Booleans
my $quiet = 0;
my $verbose = 0;
my $test= 0;
my $show_version = 0;
my $show_man = 0;
my $show_help = 0;
my $show_usage = 0;
my $print_color = 0;


#-----------------------------+
# COMMAND LINE OPTIONS        |
#-----------------------------+
my $ok = GetOptions(# REQUIRED ARGUMENTS
		    "f|fasta=s"      => \$FastaFile,
		    "a|annot|gff=s"  => \$AnnotFile,
		    "o|outfile=s"    => \$OutFile,    # Text file output
		    # OPTIONS
		    "l|te-len=s"     => \$LTRLen,
                    "t|tsd-len=s"    => \@TSDLens,
		    "gff-out=s"      => \$gff_out,    # Output in gff2 
		    "gff3-out=s"     => \$gff3_out,   # Output in gff3 format
		    "q|quiet"        => \$quiet,
		    "verbose"        => \$verbose,
		    "color"          => \$print_color,
		    # ADDITIONAL INFORMATION
		    "usage"          => \$show_usage,
		    "test"           => \$test,
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
    print "\nbatch_mask.pl:\n".
	"Version: $VERSION\n\n";
    exit;
}

#-----------------------------+
# LOAD DEFAULT TSD-LENS       |
#-----------------------------+
# If TSDLens were not provided, try 4 5 and 6
unless (@TSDLens) {
    @TSDLens = (4,5,6);
}


# Load ANSICOLOR if printing in color
if ($print_color) {
    use Term::ANSIColor;           # Allows for pretty print of ANSI output
}

#-----------------------------+
# FILEHANDLES I/O             |
#-----------------------------+
if ($OutFile) {
    open (OUTFILE, ">$OutFile") ||
	die ("Can not open outfile:\n$OutFile\n");
}
else {
    open (OUTFILE, ">&STDOUT") ||
	die "Can not print to STDOUT\n";
}

# GFF2 Format output
if ($gff_out) {
    open (GFFOUT, ">$gff_out") ||
	die ("Could not open gff outfile:\n $gff_out\n");
}

# GFF 3 Format output if desired
if ($gff3_out) {
    open (GFF3OUT, ">$gff3_out") ||
	die ("Could not open gff3 outfile:\n $gff3_out\n");
}

# Get the input sequence file in the bioperl SeqIO format
my $seqio_object = Bio::SeqIO->new( '-file' => $FastaFile,
				    '-format' => "fasta" ) ||
    die ("Can not open fasta file:\n$FastaFile");


#-----------------------------+
# LOAD ANNOTATION DATA        |
#-----------------------------+
&LoadAnnot($AnnotFile);        # Load full annotation file to 2d Array
my $LenAnnot = @Annot; 
print STDERR "$LenAnnot Annotation Objects\n" if $verbose;


# PRINT METADATA TO HEADER
print OUTFILE "# FASTA-IN:\t$FastaFile\n";
print OUTFILE "# GFF-IN:\t$AnnotFile\n";
if ($OutFile) {
    print OUTFILE "# OUT-FILE:\t$OutFile\n";
}
else {
    print OUTFILE "# OUT-FILE:\tSTDOUT\n";
}
if ($gff_out) {
    print OUTFILE "# GFF-OUT:\t$gff_out\n";
}
print OUTFILE "# LTR LENGTH:\t$LTRLen\n";

print OUTFILE "# TSD LENGTH:\t";
my $tmp_num_tsd = 1;
foreach $TSDLen (@TSDLens) {
    print OUTFILE "|" unless $tmp_num_tsd == 1;    
    print OUTFILE "$TSDLen";
    $tmp_num_tsd++;
}
print OUTFILE "\n";

#-----------------------------+
# FOR EACH SEQUENCE RECORD    |
# IN THE FASTA FILE SEARCH FOR|
# LTRS TO EXTRACT             |
#-----------------------------+
while (my $inseq = $seqio_object->next_seq) {

    $SeqId = $inseq->primary_id;
    
    if ($verbose) {
	print STDERR color 'bold';
	print STDERR "SEQUENCE: ".$SeqId."\n";
	print STDERR color 'reset';
    }

    #-----------------------------+
    # FOR EACH ANNOTATED LTR      |
    # IN THE GFF FILE             |
    #-----------------------------+
    $MaxNum = $LenAnnot - 1;
    for ($i=0; $i<=$MaxNum; $i++) {

	# The following checks that the seq_id from the gff matches the
	# 
	if ($SeqId =~ $Annot[$i][0]) {
	    
	    # Vars that do not change for TSD_len
	    my $LTR_Seq = $Annot[$i][0];
	    my $LTR_ID = $Annot[$i][8];
	    my $LTR_Start = int( $Annot[$i][3] );
	    my $LTR_End = int( $Annot[$i][4] );

	    foreach $TSDLen (@TSDLens) {
		
		if ($verbose) {
		    print STDERR "\tFeature:".$LTR_ID." : ".
			$LTR_Start."-".$LTR_End."\n";
		}		
		
		# This first pulls out TSD and portion of the LTR
		# and then just selects the TSD for the comparison
		my $LeftEnd = $inseq->subseq( $LTR_Start - $TSDLen, 
					      $LTR_Start + $LTRLen-1);
		
		my $RightEnd = $inseq->subseq( $LTR_End - $LTRLen+1,
					       $LTR_End + $TSDLen );
		
		
		my $LeftMatch = substr($LeftEnd, 0, $TSDLen);
		my $RightMatch = substr($RightEnd, $LTRLen);
		my $MatchLen = length($LeftMatch);
		
		#-----------------------------+
		# GENERATE A STRING SHOWING   |
		# MATCHING RESIDUES           |
		#-----------------------------+
		# DOING THIS MATCH ON A BASE BY BASE COMPARISION
		# SO THAT N charcaters can be avoided
		$NumTSDMatch = "0";
		$MatchString = "";
		
		for ($j=0; $j<=$MatchLen; $j++) {
		    
		    my $L = substr($LeftMatch,$j,1);
		    my $R = substr($RightMatch,$j,1);
		    
		    # DO NOT COUNT N BELOW
		    
		    if ($L =~ $R && 
			uc($L) ne "N" && 
			uc ($R) ne "N") {
			$NumTSDMatch++;
			$MatchString = $MatchString."|";  
		    }
		    else {
			$MatchString = $MatchString." "; 
		    }
		    
		}

		
		#-----------------------------+
		# PRINT GFF OUTPUT            |
		# IF POSITIVE FOR TSD MATCH   |
		#-----------------------------+
		if ($NumTSDMatch == $TSDLen) {
		    
		    $CumSum++;
		    
		    # PRINT OUTPUT IN GFF2 FORMAT
		    # ASSUMES INPUT IN GFF2
		    if ($gff_out) {
			print GFFOUT $Annot[$i][0]."\t".
			    $Annot[$i][1]."\t".
			    $Annot[$i][2]."\t".
			    $Annot[$i][3]."\t".
			    $Annot[$i][4]."\t".
			    $Annot[$i][5]."\t".
			    $Annot[$i][6]."\t".
			    $Annot[$i][7]."\t".
			    $Annot[$i][8]."\n";
		    }
		    
		    
		    # PRINT OUTPUT IF GFF3 FORMAT
		    if ($gff3_out) {
			
			# Feature ID
			my $gff3_id = $Annot[$i][8];
			$gff3_id = "ID=".$gff3_id."_".$CumSum;
			
			# Feature Name
			my $gff3_name = "Name=".$Annot[$i][8];
			
			# Features Note, separated by pipes
			my $gff3_note = "Note=tsd_length|".$TSDLen."|".
			    "tsd_seq|".$LeftMatch;
			
			# Feature Attribute
			my $gff3_attribute = $gff3_id.";".
			    $gff3_name.";".
			    $gff3_note;
			
			print GFF3OUT $Annot[$i][0]."\t".
			    $Annot[$i][1]."\t".
			    $Annot[$i][2]."\t".
			    $Annot[$i][3]."\t".
			    $Annot[$i][4]."\t".
			    $Annot[$i][5]."\t".
			    $Annot[$i][6]."\t".
			    $Annot[$i][7]."\t".
			    $gff3_attribute."\n";
		    }
		    
		    
		}  # End of output for TSD matches
		
		print OUTFILE $Annot[$i][0]."\t".
		    $Annot[$i][1]."\t".
		    $Annot[$i][2]."\t".
		    $Annot[$i][3]."\t".
		    $Annot[$i][4]."\t".
		    $Annot[$i][5]."\t".
		    $Annot[$i][6]."\t".
		    $Annot[$i][7]."\t".
		    $Annot[$i][8]."\t".
		    $NumTSDMatch."\t".
		    #$LeftEnd."\t".
		    #$RightEnd."\t".
		    #$MatchString."\t".
		    $LeftMatch."\t".
		    $RightMatch."\n";
		
		#-----------------------------+
		# PRETTY PRINT OUTPUT TO THE  |
		# SCREEN                      |
		# THIS WILL BE SORTED BY THE  |
		# BAC OCCURRENCE IN THE FASTA |
		# INPUT FILE                  |
		#-----------------------------+
		if ($verbose) {
		    if ($NumTSDMatch == $TSDLen) {
			print STDERR color 'green';
			print STDERR "\t\tMATCH:".$NumTSDMatch."\n";
			print STDERR "\t\t";
			print STDERR " " x $LTRLen;
			print STDERR $LeftEnd."\n";
			print STDERR "\t\t";
			print STDERR " " x $LTRLen;
			print STDERR $MatchString."\n";
			print STDERR "\t\t".$RightEnd."\n";
			print STDERR color 'reset';
		    }
		    else {
			print STDERR color 'red';
			print STDERR "\t\tMATCH:".$NumTSDMatch."\n";
			print STDERR "\t\t";
			print STDERR " " x $LTRLen;
			print STDERR $LeftEnd."\n";
			print STDERR "\t\t";
			print STDERR " " x $LTRLen;
			print STDERR $MatchString."\n";
			print STDERR "\t\t".$RightEnd."\n";
			print STDERR color 'reset';
		    }
		}
		
		#-----------------------------+
		# SAVE RESULTS TO ANNOTATION  |
		# ARRAY                       |
		#-----------------------------+
		$Annot[$i][9] = $NumTSDMatch;
		$Annot[$i][10] = $LeftEnd;
		$Annot[$i][11] = $RightEnd;
		$Annot[$i][12] = $MatchString;
		$Annot[$i][13] = $LeftMatch;
		$Annot[$i][14] = $RightMatch;
		
	    } # End of for each tsd_length in the array
	    
	} # End of found seq/annotation match
	
    } # End of for each annotation in the gff file

} # End of for each sequence in the output file


if ($verbose) {

    print STDERR "\n\n";
    
    print STDERR color 'bold';    
    print STDERR color 'black';
    print STDERR "CONTIG\tELEMENT\tNUM_MATCH\tLEFT\tRIGHT\t\n";
    print STDERR color 'reset';
    
    for ($i=0; $i<=$MaxNum; $i++) {
	print STDERR $Annot[$i][0]."\t";
	print STDERR $Annot[$i][8]."\t";
	print STDERR $Annot[$i][9]."\t";
	print STDERR $Annot[$i][10]."\t";
	print STDERR $Annot[$i][11]."\t";
	print STDERR $Annot[$i][13]."\t";
	print STDERR $Annot[$i][14]."\t";
	print STDERR "\n";
    }

}

#-----------------------------+
# Close output files
#-----------------------------+

if ($gff3_out) {
    close (GFF3OUT);
}


if ($gff_out){
    close (GFFOUT);
}

close (OUTFILE);

exit;

#-----------------------------------------------------------+
# SUBFUNCTIONS                                              |
#-----------------------------------------------------------+

#-----------------------------+
# LOAD BAC ANNOTATION DATA    |
#-----------------------------+
# This should be in tab delim gff format.
sub LoadAnnot {
    my $GffFile = $_[0];
    my (@line, $SeqName, $Source, $Feature, $Start, $End);
    my ($Score, $Strand, $Frame);

    open (INFILE, $GffFile) || 
	die ("Can not load file $GffFile \n");
    while (<INFILE>) {

	# Get in from tab delim text file 
	chomp;         # Get rid of newline character
	my @line = split(/\t/, $_); 
	my $Seqname = $line[0];
	my $Source = $line[1];
	my $Feature = $line[2];
	my $Start = $line[3];
	my $End = $line[4];
	my $Score = $line[5];
	my $Strand = $line[6];
	my $Frame = $line[7];
	my $Attribute = $line[8]; # Gene name column


	# Can get rid of comment below to show the loading
	# of the individual features
	#print "Load: $Attribute $Start\-$End\n";

	# Load the information into the @Annot array
	push (@Annot, [$Seqname, $Source, $Feature, $Start, $End,
		       $Score, $Strand, $Frame, $Attribute]);

    }
    close (INFILE);

    my $TestAryLen = @Annot;

    if ($TestAryLen < 1) {
	print color 'red';
	print "ERROR: I did not find any annotation data\n";
	print color 'reset';
	exit;
    }

}

#-----------------------------+
# USER VARIABLE CHECK         |
#-----------------------------+
# Let's the user check the input variables

sub UserVarCheck {
    print STDERR "\nYOU HAVE SELECTED THE FOLLOWING VARIBLES:\n";
    #-----------------------------+
    # FASTA FILE                  |
    #-----------------------------+
    print STDERR "FASTA File:\n";
    if (-e $FastaFile) {
	print STDERR "\t$FastaFile\n";
    }
    else {
	print STDERR "\tWARNING: $FastaFile \n\tDoes Not Exist\n";
    }
    
    #-----------------------------+
    # ANNOTATION FILE             |
    #-----------------------------+
    print "ANNOTATION FILE:\n";
    if (-e $AnnotFile) {
	print STDERR "\t$AnnotFile\n";
    }
    else {
	print STDERR "\tWARNING: $AnnotFile \n\tDoes Not Exist\n";
    }
    
    #-----------------------------+
    # OUTPUT FILE                 |
    #-----------------------------+
    print "OUTPUT FILE:\n";
    if (-e $OutFile) {
	print "\t$OutFile already exists\n";
	print "\tThe existing file will be overwritten\n";
    }
    else {
	print "\t$OutFile\n";
    }
    
    #-----------------------------+
    # LTR LENGTH                  |
    #-----------------------------+
    print "LTR LENGTH:\n\t$LTRLen\n";

    #-----------------------------+
    # TSD LENGTH                  |
    #-----------------------------+
    print "TSD LENGTH\n\t$TSDLen\n";
    
    #-----------------------------+
    # QUIET MODE                  |
    #-----------------------------+
    # If this is in quiet mode the answer to
    # the question is Y and the user will not
    # need to supply feedback at the command line.
    if (! $quiet){
	$ans = &UserFeedback("\nDo you want to continue (y/n)?");
    }else{
	$ans = "Y";
    }
    
    if ($ans =~"N" )
    { 
	print "Goodbye\n";
	exit;
    }elsif ($ans =~"Y"){
	print "Starting the process...\n";	
    }    
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


__END__

=head1 NAME

tsd_check.pl - Check sequence features for TSDs

=head1 VERSION

This documentation refers to program version $Rev$

=head1 SYNOPSIS

=head2 Usage

    tsd_check.pl -f seq_file.fasta -a annot_file.gff
                 -o outfile.txt
                 -gff3-out gff_outfile_gff3.gff

=head2 Required Arguments

    -f,--fasta    # Path to the input fasta file
    -a,-annot     # Path to the gff format annotation file
                  # This is assumed to be in gff2 format
    -o,--outfile  # Path to output tab delimited txt file
                  # This will print to STDOUT otherwise

=head1 DESCRIPTION

This program takes as its input a fasta file and a gff annotation file
describing a sequence that has been annotated for transposable elements
that should contain a putative target site duplication (tsd). The 
tsd_check.pl program will then check for target site duplications in the
fasta file just outside of the locations of the element indicated
by the start and stop locations in the gff intput file. The program will
search for TSDs that are 4,5, and 6 bp in length, and the list of
acceptable TSD lengths can be modified with the --tsd-len option.

The output is a text file describing the sequence of the area around
the element as well the number of matches between these sequences. It 
is also possible to return the positions of the elements that have
flanking TSDs as a gff file in the gff2 or gff3 format. The output
in the gff3 format will include the sequence and length of the 
flanking TSD elements. Using the --verbose flag will also show
the alignment in the STDERR output stream.

=head1 REQUIRED ARGUMENTS

=over 2

=item -f, --fasta

Path to the input fasta file containing the annotated elements.

=item -a, --anot

Path to the gff2 format annotation file describing the location of
the transposable elements.

=item -o,--outfile

Path of the output file providing all the data of the putative TSDs
around the annotated elements. If an outfile is not provided, the
program will print its output to STDOUT

=back

=head1 OPTIONS

=over 2

=item -l, --te-len

The length of the transposable element to include in the verbose output.

=item -t, --tsd-len

The length of the target site duplication to use. Multiple lengths of 
tsds can be supplied by using this argument more then once.

=item --gff-out

Annotation file uutput of the elements that are positive for a TSD. This
file is in the gff version 2 format.

=item --gff3-out

Annotation file output of the elements that are positive for a TSD.
This file in the the gff version 3 format. The ID field uses the 
intput id in gff2 format. The note field contains the sequence of 
the TSD and well as the length of the TSD. This data is providided 
in a format that is delimited by the pipe '|' character.

=item --verbose

Running the program in verbose format will show the test alignments
for all of the TSDs.

=item --color

This optin will print the --verbose output in color. In this format the
sequence alignments that do not have TSDs are printed in red and those
sequence alignments that do have TSDs are printed in green. This option
makes use of the Term::ANSIColor perl module.

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

The typical use of this program will be to search elements annotated 
in a genome and see if they contain target site dupliations. This output
would then typically be saved in gff3 output format.

    tsd_check.pl -f chrom01.fasta -a chrom01_ltrs.gff
                 -o chrom01_ltr_tsd.txt
                 -gff3-out chrom01_ltrs_wtsd.gff

=head 2 Specifying List of TSD Lengths

It is also possible to specify a list of Target Site Duplication
lenghts to use by using the t option multiple times at the 
command line. For example, if you wanted to check for TSDs that
were 3,4,5, or 6 bp in length:

    tsd_check.pl -f chrom01.fasta -a chrom01_ltrs.gff
                 -t 3 -t 4 -t 5 -t 6
                 -o chrom01_ltr_tsd.txt
                 -gff3-out chrom01_ltrs_wtsd.gff

=head 2 Print alignment output in color

The program will print alignment of the TSDs if the --verbose option
is used. It is also possible to print alignments without matching
TSDs in red and TSDs with matching TSDs in green. This can be invoked
with the --color option. For example:

    tsd_check.pl -f chrom01.fasta -a chrom01_ltrs.gff
                 -o chrom01_ltr_tsd.txt
                 -gff3-out chrom01_ltrs_wtsd.gff
                 --verbose --color

=head1 DIAGNOSTICS

=over 2

=item * Can not print to STDOUT

The program is attempting to print outout to the STDOUT screen. If
you are having this error, try specifying an output path with the
-o or --output option.

=back

=head1 CONFIGURATION AND ENVIRONMENT

This program does not make use of configuration files or options
set in the user environment.

=head1 DEPENDENCIES

=head2 Software

This program is not dependent on external software.

=head2 Perl Modules

This program makes use of the following PERL modules.

=over2

=item Term::ANSIColor

This module is used to print alignment output in color. This module
will only be used if the --color option is specified at the 
command line.

=item Bio::SeqIO

The SeqIO module of Bioperl is used to parse subregions of the
input sequence.

=back

=head1 BUGS AND LIMITATIONS

This program does not currently contain any known bugs. A current
limitation is that the -t option can not accept a list as input;
the values must be specified multiple times usin the -t option
for each value to check.

=head1 REFERENCE

A manuscript is being submitted describing the DAWGPAWS program. 
Until this manuscript is published, please refer to the DAWGPAWS 
SourceForge website when describing your use of this program:

JC Estill and JL Bennetzen. 2009. 
The DAWGPAWS Pipeline for the Annotation of Genes and Transposable 
Elements in Plant Genomes.
http://dawgpaws.sourceforge.net/

=head1 LICENSE

GNU General Public License, Version 3

L<http://www.gnu.org/licenses/gpl.html>

=head1 AUTHOR

James C. Estill E<lt>JamesEstill at gmail.comE<gt>

=head1 HISTORY

STARTED: 06/08/2006

UPDATED: 05/21/2009

VERSION: $Rev$

=cut

#-----------------------------------------------------------+
# HISTORY                                                   |
#-----------------------------------------------------------+
#
# 06/08/2006
# - Begin modification of DblHit to do a dts_sel for clem
#
# 05/15/2009
# - Modifiying to long intput format
# - Removing user feedback to facilitate use on r cluster
# - Added output of fetures with tsds to a separate file
#
# 05/18/2009
# - Added count of cumulative sum
# - Added option to provide output in GFF format
#     - This assumes input is GFF2 format
#     - This will append cum_sum to input ID to insure unique id
#     - This will give information on the TSD in the Note field 
#       of the GFF3 formatted attributes file. This note 
#       sep by pipe characters
# 
# 05/21/2009
# - Adding POD documentation
#
#-----------------------------------------------------------+
# MODEL                                                     |
#-----------------------------------------------------------+
#
# Given a genomic sequence and Start/End Positions of the
# LTR Retrotransposon for a given genomic sequence:
#
# GENOME
# |=========================================================================|
#                    
#                  START                           END
#                   ||                             ||
# LTR               \/                             \/
#                    [=5'LTR=]-------------[3'=LTR=]
#                                      
# RETURNED          
# SEQS              =:                              :=
#                  5'Seq                             3'Seq
#
# GIVEN: 
# START is position of start of 5' LTR
# END is position of end of 3' LTR
