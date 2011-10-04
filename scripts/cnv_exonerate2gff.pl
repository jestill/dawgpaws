#!/usr/bin/perl -w
#-----------------------------------------------------------+
#                                                           |
# cnv_exonerate2gff.pl - Convert Exonerate output GFF       |
#                                                           |
#-----------------------------------------------------------+
#                                                           |
#  AUTHOR: James C. Estill                                  |
# CONTACT: JamesEstill_@_gmail.com                          |
# STARTED: 09/30/2011                                       |
# UPDATED: 10/4/2011                                        |
#                                                           |
# DESCRIPTION:                                              |
#  Convert exonerate output to GFF format output. Exonerate |
#  is useful for generating annotation evidences for genes  |
#  since it is an intro-aware alignment algorithm.          |
#                                                           |
# USAGE:                                                    |
#  cnv_exonerate2gff.pl -i infile.txt -o outfile.gff        |
#                                                           |
# VERSION: $Rev$                                     |
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
# Get GFF version from environment, GFF2 is DEFAULT
my $gff_ver = uc($ENV{DP_GFF}) || "GFF2";

#-----------------------------+
# VARIABLE SCOPE              |
#-----------------------------+
my $infile;
my $outfile;

# BOOLEANS
my $quiet = 0;
my $verbose = 0;
my $show_help = 0;
my $show_usage = 0;
my $show_man = 0;
my $show_version = 0;
my $do_test = 0;                  # Run the program in test mode
my $do_append = 0;

my $qry_name;
my $feature_type;
my $program_name;
my $param_name;
my $delim_char = ":";

#-----------------------------+
# COMMAND LINE OPTIONS        |
#-----------------------------+
my $ok = GetOptions(# REQUIRED OPTIONS
		    "i|infile=s"     => \$infile,
                    "o|outfile=s"    => \$outfile,
		    # ADDITIONAL OPTIONS
		    # Not recommended for exonerate
		    "n|s|seqname|name=s"  => \$qry_name,
		    "delim=s"        => \$delim_char,
		    "gff-ver=s"      => \$gff_ver,
		    "q|quiet"        => \$quiet,
		    "verbose"        => \$verbose,
		    "f|feature=s"    => \$feature_type,
		    "p|program=s"    => \$program_name,
		    "param=s"        => \$param_name,
		    "append"         => \$do_append,
		    # ADDITIONAL INFORMATION
		    "usage"          => \$show_usage,
		    "test"           => \$do_test,
		    "version"        => \$show_version,
		    "man"            => \$show_man,
		    "h|help"         => \$show_help,);



 
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

#-----------------------------+
# MAIN PROGRAM BODY           |
#-----------------------------+
exonerate2gff ( $infile, $outfile, $do_append, $qry_name,
    $program_name, $param_name, $delim_char, $feature_type);

exit 0;

#-----------------------------------------------------------+ 
# SUBFUNCTIONS                                              |
#-----------------------------------------------------------+

# Basic subfunction to convert exonerate output to GFF3 format
# for use with Apollo or EVM
sub exonerate2gff {

    my ($exin, $gffout, $append, $seqname, $prog, $suffix, $d, $feature) = @_;

    #-----------------------------+
    # INPUT FROM FILE PATH OR     |
    # STDINT OTHERWISE            |
    #-----------------------------+
    if ($exin) {
	open (EXIN, $exin) ||
	    die "Con not open input file $exin\n";
    }
    else {
	print STDERR "Expecting input from STDIN\n";
	open (EXIN, "<&STDIN") ||
	    die "Can not access STDIN for input";
    }


    #-----------------------------+
    # OUTPUT TO FILE PATH OR      |
    # STDOUT OTHERWISE            |
    #-----------------------------+
    if ($gffout) {
	open (GFFOUT, $gffout) ||
	    die "Can not open $gffout for output\n";
    }
    else {
    	open (GFFOUT, ">&STDOUT") ||
	    die "Can not open STDOUT for writing\n";
    }
    # PRINT GFF3 HEADER
    if ($gff_ver =~ "GFF3") {
	print GFFOUT "##gff-version 3\n";
    }    


    # Vars used across lines in the input GFF file
    my $feat_match_id;
    my $feat_gene_score;
    my $gene_num = 0;      # Appending gene_number keeps unique_id
                           # Used for cases with more than one match
                           # per scaffold or contig
    
    # PROCESS THE INPUT FILE
    while (<EXIN>) {
	
	# Need to pay attantion of when we are startin 
	# and ending GFF information, and will also
	# may need to fetch some identifying information from
	# the first occurence of gene in the input file
	
	
	my $feat_seqid;
	my $feat_source;
	my $feat_type;
	my $feat_start;
	my $feat_end;
	my $feat_score;
	my $feat_strand;
	my $feat_frame;
	my $feat_attributes;
	my @attribute_parts;
	
	
	# Skip non-gff lines
	next if m/^\#/; 
	next if m/^\>/; 
	next if m/^vulgar\:/; 
	chomp;
	
	
#	print STDERR $_."\n" 
#	    if $verbose;
	
	# Get the poarts from the GFF string
	my @gff2_parts = split (/\t/, $_);
	my $num_gff2_parts = @gff2_parts;
#	print STDERR $num_gff2_parts."\n" 
#	    if $verbose;
	
	# Do test for expected number of results
	# then load featur parts to usable variabl names
	# PROBLEM: utr5 feature does not have 9 parts, just 8
	if ($num_gff2_parts == 9) {
	    $feat_seqid = $gff2_parts[0];
	    $feat_source = $gff2_parts[1];
	    $feat_type  = $gff2_parts[2];
	    $feat_start  = $gff2_parts[3];
	    $feat_end  = $gff2_parts[4];
	    $feat_score  = $gff2_parts[5];
	    $feat_strand  = $gff2_parts[6];
	    $feat_frame  = $gff2_parts[7];
	    $feat_attributes  = $gff2_parts[8];
	}
	else {
	    next;
	}
	
	# override file vars with vars from cmd line
	# this allows for using vars as output by program, or
	# replacing with custom variables.
	if ($prog) {
	    $feat_source = $prog;
	}
	
	if ($suffix) {
	    $feat_source = $feat_source.$d.$suffix;
	}
	
	# For the first go will just use exon and gene features
	# since that is what the EVM parser uses
	if ( $feat_type =~ "exon" ||
	     $feat_type =~ "gene") {

	    if ($feat_type =~ "gene") {
		$gene_num++;
		$feat_gene_score = $feat_score;
	    }
	    if ($feat_attributes =~ m/sequence (.*);/) {
		$feat_match_id = $1;
		# Remove trailing whitespaces
		$feat_match_id =~ s/\s+$//;
	    }
	    
	    #-----------------------------------------------------------+
	    # Currently just print gff string for 'exons'
	    #-----------------------------------------------------------+
	    if ($feat_type =~ "exon") {
		
		$feat_score = $feat_gene_score;

		$feat_attributes = "Id=".$feat_seqid."_".
		    "ex_".$feat_match_id."_".$gene_num.
		    ";Name=".$feat_match_id.
		    ";Target=".$feat_match_id;

		
		$feat_type = "match";
		# Allow override of feature type at cmd line
		if ($feature) {
		    $feat_type = $feature;
		}
		
		my $gff_string = $feat_seqid."\t".
		    $feat_source."\t".
		    $feat_type."\t".
		    $feat_start."\t".
		    $feat_end."\t".
		    $feat_score."\t".
		    $feat_strand."\t".
		    $feat_frame."\t".
		    $feat_attributes.
		    "\n";
		
		print GFFOUT $gff_string;
	    }
	    
	} # End of if gene or exon
    } # End of processing each line of input file
    
} # subfunction end

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

cnv_exonerate2gff.pl -  Convert Exonerate output to GFF

=head1 VERSION

This documentation refers to program version $Rev$

=head1 SYNOPSIS

=head2 Usage

    cnv_exonerate2gff.pl -i InFile -o OutFile

=head2 Required Arguments

    --infile        # Path to exonerate output file to convert
    --outfie        # Path to the output GFF3 format ile

=head1 DESCRIPTION

This take exonerate output and converts the exonerate style GFF2 format
to a GFF3 format. Currently only the 'exon' matches are returned,
and these are converted to match features.

=head1 REQUIRED ARGUMENTS

=over 2

=item -i,--infile

Path of the input file. If no input file path is given, the program will
expect input from STDIN.

=item -o,--outfile

Path of the output file. If not output file is given, the program will
write output to STDOUT.

=back

=head1 OPTIONS

=over 2

=item -n,--seqname

Name of the sequence contig. This will overwrite the names as used
in the exonerate file.

=item -p,--program

The program name to use in the second colum of the GFF output file.

=item --param

The parameter name to use in the second column of the GFF output file.
This can be used to identify use of the same program using different
parameter settings.

=item --delim

The delimiting character to use to separate program and parameter name
in the second column of the GFF output file. By default the colon ':'
character is used for delimiting program and parameter settings.

=item -f,--feature

The feature name to use for the third column of the GFF output file.
By default, these will be set to 'match' but this command line options
allows for an override of this default.

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

The following are example uses of the program.

=head2 Using all the default options:

 cnv_exonerate2gff.pl -i test_protein2genome.exonerate 

This will produce output like the following:

 scaffold00095	exonerate:protein2genome:local	match	1515255	1515425	361	-	.	Id=scaffold00095_ex_AT1G09500.2|PACid:19651388_37; Name=AT1G09500.2|PACid:19651388
 scaffold00095	exonerate:protein2genome:local	match	1495360	1495545	361	-	.	Id=scaffold00095_ex_AT1G09500.2|PACid:19651388_37; Name=AT1G09500.2|PACid:19651388
 scaffold00095	exonerate:protein2genome:local	match	1494229	1494391	361	-	.	Id=scaffold00095_ex_AT1G09500.2|PACid:19651388_37; Name=AT1G09500.2|PACid:19651388
 scaffold00095	exonerate:protein2genome:local	match	1493525	1493697	361	-	.	Id=scaffold00095_ex_AT1G09500.2|PACid:19651388_37; Name=AT1G09500.2|PACid:19651388
 scaffold00015	exonerate:protein2genome:local	match	2840379	2840552	342	-	.	Id=scaffold00015_ex_AT1G09500.2|PACid:19651388_38; Name=AT1G09500.2|PACid:19651388
 scaffold00015	exonerate:protein2genome:local	match	2838517	2838705	342	-	.	Id=scaffold00015_ex_AT1G09500.2|PACid:19651388_38; Name=AT1G09500.2|PACid:19651388
 scaffold00015	exonerate:protein2genome:local	match	2836318	2836488	342	-	.	Id=scaffold00015_ex_AT1G09500.2|PACid:19651388_38; Name=AT1G09500.2|PACid:19651388

=head2 Override the program name

It is possible to set the program name at the command line:

 cnv_exonerate2gff.pl -i test_protein2genome.exonerate --program exonerate

This will produce output like the following: 

 scaffold00095	exonerate	match	1515255	1515425	361	-	.	Id=scaffold00095_ex_AT1G09500.2|PACid:19651388_37; Name=AT1G09500.2|PACid:19651388
 scaffold00095	exonerate	match	1495360	1495545	361	-	.	Id=scaffold00095_ex_AT1G09500.2|PACid:19651388_37; Name=AT1G09500.2|PACid:19651388
 scaffold00095	exonerate	match	1494229	1494391	361	-	.	Id=scaffold00095_ex_AT1G09500.2|PACid:19651388_37; Name=AT1G09500.2|PACid:19651388
 scaffold00095	exonerate	match	1493525	1493697	361	-	.	Id=scaffold00095_ex_AT1G09500.2|PACid:19651388_37; Name=AT1G09500.2|PACid:19651388
 scaffold00015	exonerate	match	2840379	2840552	342	-	.	Id=scaffold00015_ex_AT1G09500.2|PACid:19651388_38; Name=AT1G09500.2|PACid:19651388
 scaffold00015	exonerate	match	2838517	2838705	342	-	.	Id=scaffold00015_ex_AT1G09500.2|PACid:19651388_38; Name=AT1G09500.2|PACid:19651388
 scaffold00015	exonerate	match	2836318	2836488	342	-	.	Id=scaffold00015_ex_AT1G09500.2|PACid:19651388_38; Name=AT1G09500.2|PACid:19651388

=head2 Identify the parameter used 

It is also possible to identify the paramters used or the database searched
against by appending a parameter name to the program column.

 cnv_exonerate2gff.pl -i test_protein2genome.exonerate --program exonerate --param arabidopsis

This will produce results like:

 scaffold00095	exonerate:arabidopsis	match	1515255	1515425	361	-	.	Id=scaffold00095_ex_AT1G09500.2|PACid:19651388_37; Name=AT1G09500.2|PACid:19651388
 scaffold00095	exonerate:arabidopsis	match	1495360	1495545	361	-	.	Id=scaffold00095_ex_AT1G09500.2|PACid:19651388_37; Name=AT1G09500.2|PACid:19651388
 scaffold00095	exonerate:arabidopsis	match	1494229	1494391	361	-	.	Id=scaffold00095_ex_AT1G09500.2|PACid:19651388_37; Name=AT1G09500.2|PACid:19651388
 scaffold00095	exonerate:arabidopsis	match	1493525	1493697	361	-	.	Id=scaffold00095_ex_AT1G09500.2|PACid:19651388_37; Name=AT1G09500.2|PACid:19651388
 scaffold00015	exonerate:arabidopsis	match	2840379	2840552	342	-	.	Id=scaffold00015_ex_AT1G09500.2|PACid:19651388_38; Name=AT1G09500.2|PACid:19651388
 scaffold00015	exonerate:arabidopsis	match	2838517	2838705	342	-	.	Id=scaffold00015_ex_AT1G09500.2|PACid:19651388_38; Name=AT1G09500.2|PACid:19651388
 scaffold00015	exonerate:arabidopsis	match	2836318	2836488	342	-	.	Id=scaffold00015_ex_AT1G09500.2|PACid:19651388_38; Name=AT1G09500.2|PACid:19651388

=head2 Change the delimiting character

You can change the delmiting charater with the --delim option:

 cnv_exonerate2gff.pl -i test_protein2genome.exonerate --program exonerate --param arabidopsis --delim _

To prodouce results like:

 scaffold00095	exonerate_arabidopsis	match	1515255	1515425	361	-	.	Id=scaffold00095_ex_AT1G09500.2|PACid:19651388_37; Name=AT1G09500.2|PACid:19651388
 scaffold00095	exonerate_arabidopsis	match	1495360	1495545	361	-	.	Id=scaffold00095_ex_AT1G09500.2|PACid:19651388_37; Name=AT1G09500.2|PACid:19651388
 scaffold00095	exonerate_arabidopsis	match	1494229	1494391	361	-	.	Id=scaffold00095_ex_AT1G09500.2|PACid:19651388_37; Name=AT1G09500.2|PACid:19651388
 scaffold00095	exonerate_arabidopsis	match	1493525	1493697	361	-	.	Id=scaffold00095_ex_AT1G09500.2|PACid:19651388_37; Name=AT1G09500.2|PACid:19651388
 scaffold00015	exonerate_arabidopsis	match	2840379	2840552	342	-	.	Id=scaffold00015_ex_AT1G09500.2|PACid:19651388_38; Name=AT1G09500.2|PACid:19651388
 scaffold00015	exonerate_arabidopsis	match	2838517	2838705	342	-	.	Id=scaffold00015_ex_AT1G09500.2|PACid:19651388_38; Name=AT1G09500.2|PACid:19651388
 scaffold00015	exonerate_arabidopsis	match	2836318	2836488	342	-	.	Id=scaffold00015_ex_AT1G09500.2|PACid:19651388_38; Name=AT1G09500.2|PACid:19651388

=head2 Change the feature type

You can change the feature type with the --feat option

 cnv_exonerate2gff.pl -i test_protein2genome.exonerate --program exonerate --feat match_part


To produce results like:

 scaffold00095	exonerate	match_part	1515255	1515425	361	-	.	Id=scaffold00095_ex_AT1G09500.2|PACid:19651388_37; Name=AT1G09500.2|PACid:19651388
 scaffold00095	exonerate	match_part	1495360	1495545	361	-	.	Id=scaffold00095_ex_AT1G09500.2|PACid:19651388_37; Name=AT1G09500.2|PACid:19651388
 scaffold00095	exonerate	match_part	1494229	1494391	361	-	.	Id=scaffold00095_ex_AT1G09500.2|PACid:19651388_37; Name=AT1G09500.2|PACid:19651388
 scaffold00095	exonerate	match_part	1493525	1493697	361	-	.	Id=scaffold00095_ex_AT1G09500.2|PACid:19651388_37; Name=AT1G09500.2|PACid:19651388
 scaffold00015	exonerate	match_part	2840379	2840552	342	-	.	Id=scaffold00015_ex_AT1G09500.2|PACid:19651388_38; Name=AT1G09500.2|PACid:19651388
 scaffold00015	exonerate	match_part	2838517	2838705	342	-	.	Id=scaffold00015_ex_AT1G09500.2|PACid:19651388_38; Name=AT1G09500.2|PACid:19651388
 scaffold00015	exonerate	match_part	2836318	2836488	342	-	.	Id=scaffold00015_ex_AT1G09500.2|PACid:19651388_38; Name=AT1G09500.2|PACid:19651388



=head1 DIAGNOSTICS

=over 2

=item * Expecting input from STDIN

If you see this message, it may indicate that you did not properly specify
the input sequence with -i or --infile flag. 

=back

=head1 CONFIGURATION AND ENVIRONMENT

This program does not make use of a configuration file or rely on 
options set in the user environment.

=head1 DEPENDENCIES

This program does not have any dependencies beyound the modules typically
installed with Perl.

=head1 BUGS AND LIMITATIONS

This program currently only produces GFF3 output.

=head1 REFERENCE

To cite your use of DAWGPAWS programs:

JC Estill and JL Bennetzen. 2009. "The DAWGPAWS pipeline for the annotation 
of genes and transposable elements in plant genomes." Plant Methods. 5:8.

=head1 LICENSE

GNU General Public License, Version 3

L<http://www.gnu.org/licenses/gpl.html>

=head1 AUTHOR

James C. Estill E<lt>JamesEstill at gmail.comE<gt>

=head1 HISTORY

STARTED: 09/30/2011

UPDATED: 10/04/2011

VERSION: $Rev$

=cut

#-----------------------------------------------------------+
# HISTORY                                                   |
#-----------------------------------------------------------+
#
