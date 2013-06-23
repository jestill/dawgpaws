#!/usr/bin/perl -w
#-----------------------------------------------------------+
#                                                           |
# batch_blast2seq_.pl
#                                                           |
#-----------------------------------------------------------+
#                                                           |
#  AUTHOR: James C. Estill                                  |
# CONTACT: JamesEstill_@_gmail.com                          |
# STARTED: 09/07/2012                                       |
# UPDATED: 09/11/2012                                       |
#                                                           |
# DESCRIPTION:                                              |
#  Do needle alignment and report values that exceed        |
#  threshold.                                               |
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

package DAWGPAWS;

#-----------------------------+
# INCLUDES                    |
#-----------------------------+
use GD;                        # The gd program will be used to draw the output matrix
use strict;
use Getopt::Long;
# The following needed for printing help
use Pod::Select;               # Print subsections of POD documentation
use Pod::Text;                 # Print POD doc as formatted text file
use IO::Scalar;                # For print_help subfunction
use IO::Pipe;                  # Pipe for STDIN, STDOUT for POD docs
use File::Spec;                # Convert a relative path to an abosolute path

use Bio::SearchIO;             # Parse BLAST output
use Bio::SeqIO;                # use to get the sequence length
#-----------------------------+
# PROGRAM VARIABLES           |
#-----------------------------+
my ($VERSION) = q$Rev$ =~ /(\d+)/ || "pre-release";
# Get GFF version from environment, GFF2 is DEFAULT
my $gff_ver = uc($ENV{DP_GFF}) || "GFF2";

#-----------------------------+
# VARIABLE SCOPE              |
#-----------------------------+
my $indir;
my $outdir;
my $summary_out;
my $img_out;
my $img_thickness = 2;

# BOOLEANS
my $quiet = 0;
my $verbose = 0;
my $show_help = 0;
my $show_usage = 0;
my $show_man = 0;
my $show_version = 0;
my $do_test = 0;                  # Run the program in test mode
my $do_keep_temp = 0;             # Keep the temp file
my $do_image = 0;                 # Boolean for drawing the image of the graph
my $img_scale = 1;
# The threoshold value to report in the matrix

my $threshold = 80;

#-----------------------------+
# COMMAND LINE OPTIONS        |
#-----------------------------+
my $ok = GetOptions(# REQUIRED OPTIONS
		    "i|indir=s"     => \$indir,
                    "o|outdir=s"    => \$outdir,
		    # outfile for summary output
		    "s|summary=s"    => \$summary_out,
		    # ADDITIONAL OPTIONS
		    "image"         => \$do_image,
		    "thickness=s"   => \$img_thickness,
		    "scale=s"       => \$img_scale,
		    "k|keep-temp"   => \$do_keep_temp,
		    "t|threshold=s" => \$threshold,
		    "gff-ver=s"     => \$gff_ver,
		    "q|quiet"       => \$quiet,
		    "verbose"       => \$verbose,
		    # ADDITIONAL INFORMATION
		    "usage"         => \$show_usage,
		    "test"          => \$do_test,
		    "version"       => \$show_version,
		    "man"           => \$show_man,
		    "h|help"        => \$show_help,);

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
    print "\nbatch_mask.pl:\n".
	"Version: $VERSION\n\n";
    exit;
}

if ($summary_out) {
    open (SUMOUT, ">$summary_out") ||
	die "Can not open matrix output file\n";
}
else {
    open (SUMOUT, ">&STDOUT") ||
	die "Can not open STDOUT for writing\n";
}

#-----------------------------+
# MAIN PROGRAM BODY           |
#-----------------------------+


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
my @fasta_files = grep /\.fasta$|\.fa$|\.fna$/, readdir DIR ;
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
# CREATE OUTPUT DIR IF NEEDED
#-----------------------------+
unless (-e $outdir) {
    mkdir $outdir ||
	die "Could not create dir:\n$outdir\n"
}


#-----------------------------+
# DO AN ALIGNMENT             |
# FOR EACH SEQUENCE IN THE    |
# INPUT DIR                   |
#-----------------------------+
my $file_num = 0;
my $proc_num = 0;

for my $ind_a_file (@fasta_files) {

    $file_num++;
    print STDERR "Processing ".$file_num." of ".$count_files."\n";
    my $a_name_root;
#    my $a_seq_length;

    # Get root file name
    if ($ind_a_file =~ m/(.*)\.hard\.fasta$/) {
	# file ends in .hard.fasta
	$a_name_root = "$1";
    }
    elsif ($ind_a_file =~ m/(.*)\.masked\.fasta$/) {
	# file ends in .masked.fasta
	$a_name_root = "$1";
    }
    elsif ($ind_a_file =~ m/(.*)\.fasta$/ ) {	    
	# file ends in .fasta
	$a_name_root = "$1";
    }  
    elsif ($ind_a_file =~ m/(.*)\.fa$/ ) {	    
	# file ends in .fa
	$a_name_root = "$1";
    } 
    elsif ($ind_a_file =~ m/(.*)\.fna$/ ) {	    
	# file ends in .fa
	$a_name_root = "$1";
    } 
    else {
	$a_name_root = $ind_a_file;
    }
    
    
    # Get sequece information
    my $a_file_path = $indir.$ind_a_file;
    my $a_seq_in  = Bio::SeqIO->new(-file => $a_file_path , 
				    '-format' => 'Fasta');

    my $a_seq_obj = $a_seq_in->next_seq();
    my $a_seq_length = $a_seq_obj->length;
    my $a_seq_id = $a_seq_obj->id();
    
    # Need the seq length to compare the sum of the tiled alignment to the 
    # original sequence length

    # For now will use a as both query and subject
    my $seq_path = $indir.$ind_a_file;
    my $blast_outfile = $outdir.$a_name_root.".xml";
    # The following uses the options for blast2seq
    # as set by default in the NCBI GUI for blast2seq
    # using the blastn option
    my $blast_cmd = "blastn".
	" -query ".$seq_path.
	" -subject ".$seq_path.
	" -out ".$blast_outfile.
	# output in xml
	" -outfmt 5".
	" -evalue 10".
	" -word_size 11".
	" -gapopen 5".
	" -gapextend 2".
	" -penalty -3".
	" -reward 2";

    # Can write this to an output file
    # or optionally could write to STDOUT
    # and pipe results 
    
    print STDERR "\n".$blast_cmd."\n"
	if $verbose;

    unless ($do_test) {
	system ($blast_cmd);
    }


    # process the blast result
#    my $blast_report = new Bio::SearchIO ( '-format' => 'blast',
    my $blast_report = new Bio::SearchIO ( '-format' => 'blastxml',
					   '-file'   => $blast_outfile) ||
					       die "Cant open $blast_outfile\n";
    
    my $sum_hsp_len = 0;
    my $qry_seq_id;

    # Set the image canvas for drawing the image
    my $img;
    my $white;
    my $black;
    my $red;
    my $green;
    my $blue;
    my $gray;

    if ( $do_image ) {
	my $max_x = $a_seq_length * $img_scale;
	my $max_y = $a_seq_length * $img_scale;

	$img = new GD::Image($max_x, $max_y);
	$img->setThickness($img_thickness);

	
	# allocate colors for later use
	$white = $img->colorAllocate   ( 255,  255, 255 );
	$gray = $img->colorAllocate    ( 100,  100, 100 );
	$black = $img->colorAllocate   (   0,    0,   0 );
	$red = $img->colorAllocate     ( 255,    0,   0 );
	$green = $img->colorAllocate   (   0,  255,   0 );
	$blue  = $img->colorAllocate   (   0,    0, 255 );

	# Set file path for ouput
	# Will do png for now
	$img_out = $outdir.$a_name_root.".".$a_name_root.".png";
	print STDERR "Image to:".$img_out."\n";
	
    }

    while (my $blast_result = $blast_report->next_result() ) {
	$qry_seq_id = $blast_result->query_name;
	while (my $blast_hit = $blast_result->next_hit() ) {
	    
	    while (my $blast_hsp = $blast_hit->next_hsp()) {

		#-----------------------------+
		# GET QRY START AND END       |
		#-----------------------------+
		# Make certain that start coordinate is
		# less then the end coordinate
		my $start;
		my $end;
		if ( $blast_hsp->start() < $blast_hsp->end() ) {
		    $start = $blast_hsp->start();
		    $end = $blast_hsp->end();
		}
		else {
		    $start = $blast_hsp->end();
		    $end = $blast_hsp->start();
		}
		my $qry_strand = $blast_hsp->strand('query');
		my $qry_len = $end - $start + 1;
		$sum_hsp_len = $sum_hsp_len + $qry_len;

		# The following are reported as 1 and -1
		my $qry_strand = $blast_hsp->strand('query');
		my $hit_strand = $blast_hsp->strand('hit');
		print STDERR "QRY STRAND: ".$qry_strand.":".
		    " HIT STRAND:".$hit_strand."\n";
		print STDERR "\tQRY:".$blast_hsp->start('query').":".
		    $blast_hsp->end('query');
		print STDERR "\tHIT:".$blast_hsp->start('hit').":".
		    $blast_hsp->end('hit')."\n";
		#-----------------------------+
		# GET SUBJECT START AND END   |
		#-----------------------------+
		#my $hit_start;
		#my $hit_end;
		#if ( $blast_hsp->start('hit') < $blast_hsp->end('hit') ) {
		#    $hit_start = $blast_hsp->start('hit');
		#    $hit_end = $blast_hsp->end('hit');
		#}
		#else {
		#    $hit_start = $blast_hsp->end('hit');
		#    $hit_end = $blast_hsp->start('hit');
		#}


		# Draw a LINE for the HSP if requested
		if ( $do_image ) {
		    # THE FOLLOWING ASSUMES THAT Q
		    # Need strand information for this
		    # image is as $image->line($x1,$y1,$x2,$y2,$color)
		    # since the test is a self image can do this:
		    # Could mulitiply to determine if same or different
		    if ( $blast_hsp->strand('query') == 1) {
			if ( $blast_hsp->strand('hit') == 1) {
			    $img->line( $blast_hsp->start(),
					$blast_hsp->start('hit'), 
					$blast_hsp->end, 
					$blast_hsp->end('hit'), 
					$black );
			    
			}
			else {
			    $img->line( $blast_hsp->start(),
					$blast_hsp->end('hit'), 
					$blast_hsp->end, 
					$blast_hsp->start('hit'), 
					$black );
			    
			}
		    }
		    # query on negative strand
		    else {
			if ( $blast_hsp->strand('hit') == 1) {
			    $img->line( $blast_hsp->start(),
					$blast_hsp->end('hit'), 
					$blast_hsp->end, 
					$blast_hsp->start('hit'), 
					$black );
			}
			else {
			    $img->line( $blast_hsp->start(),
					$blast_hsp->start('hit'), 
					$blast_hsp->end, 
					$blast_hsp->end('hit'), 
					$black );
			}
			
		    }
		}
		
	    } # End of next hsp 
	} # End of next hit
    } # End of next result


    #-----------------------------+
    # WRITE IMAGE OUTPUT TO FILE  |
    #-----------------------------+

#    $image = $sourceImage->copyFlipHorizontal();
#    $image = $sourceImage->copyFlipVertical();
    # Test of a transpose of the image to make it looke
    # like what fols expect from BLAST2Seq visualization
    # on NCBI BLAST2Seq
#    my $img_transpose = $img->copyTranspose();
    my $img_transpose = $img->copyFlipHorizontal();

	# 
    if ( $do_image ) {
	open (IMGOUT, ">$img_out") ||
	    die "Can not open $img_out";
	binmode IMGOUT;
#	print IMGOUT $img->png;
	print IMGOUT $img_transpose->png;
	close IMGOUT;
    }

    # Get summary information
    my $hsp_ratio = $sum_hsp_len/$a_seq_length;
    $hsp_ratio = sprintf("%.3f", $hsp_ratio);
    
    # Write summary output
    if ($verbose) {
	print STDERR "\n\n";
	print STDERR "Sequence ID:".$qry_seq_id."\n";
	print STDERR "Sum HSP length: ".$sum_hsp_len."\n";
	print STDERR "HSP Length ratio:".$hsp_ratio."\n";
    }
    # Write output to summary file
#    $a_seq_id
#    print SUMOUT $qry_seq_id."\t".$hsp_ratio."\n";
    print SUMOUT $a_seq_id."\t".$hsp_ratio."\n";

}


close (SUMOUT);
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

batch_blast2seq_self.pl - blast2seq for all seqs in dir

=head1 VERSION

This documentation refers to program version 0.1

=head1 SYNOPSIS

=head2 Usage

    batch_pairwise_global_align.pl -i indir -o outdir -s summary.txt

=head2 Required Arguments

    --indir,i          # Path to the output dir
    --outdir,-o        # Path to the input dir
    --summary,-s       # Summary file output

=head1 DESCRIPTION

For a directory of FASTA files, run blast2seq against the indiviudal
seuqneces and report the length of the HSPs against the length of the 
sequences. For completely unique sequences, this will be 1.0. As the
sequences have lots of internal repeats, this number will be larger
than one. The basic idea is that this will be useful information when
screening for false positives of LTR retro annotation.

=head1 REQUIRED ARGUMENTS

=over 2

=item -i,--indir

Path to the input dir containing the fasta files to align.

=item -o,--outdir

Path to the dir where the blast reports will be stored.

=item -s, --summary

Path to the summary file that gives the HSP ratio for each sequence.

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

=head1 EXAMPLES

The following are examples of how to use this script

=head2 Typical Use

To run the program with default settings

 batch_pairwise_global_align.pl -i input_dir -o output_dir -m matrix_out.txt

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

=head2 Software

=over

=item EMBOSS::needle

This program requires the needle global alignment program from EMBOSS.
http://emboss.sourceforge.net/

=back

=head1 BUGS AND LIMITATIONS

Any known bugs and limitations will be listed here.

=head1 REFERENCE

Please refer to the DAWGPAWS manuscript in Plant Methods when describing
your use of this program:

JC Estill and JL Bennetzen. 2009. 
"The DAWGPAWS Pipeline for the Annotation of Genes and Transposable 
Elements in Plant Genomes." Plant Methods. 5:8.

=head1 LICENSE

GNU General Public License, Version 3

L<http://www.gnu.org/licenses/gpl.html>

=head1 AUTHOR

James C. Estill E<lt>JamesEstill at gmail.comE<gt>

=head1 HISTORY

STARTED: 09/07/2012

UPDATED: 09/07/2012

VERSION: $Rev$

=cut

#-----------------------------------------------------------+
# HISTORY                                                   |
#-----------------------------------------------------------+
#
