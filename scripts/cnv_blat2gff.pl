#!/usr/bin/perl -w
#-----------------------------------------------------------+
#                                                           |
# cnv_blat2gff.pl - Convert BLAST output to gff             |
#                                                           |
#-----------------------------------------------------------+
#  AUTHOR: James C. Estill                                  |
# CONTACT: JamesEstill_at_gmail.com                         |
# STARTED: 07/06/2006                                       |
# UPDATED: 09/13/2010
#                                                           |  
# DESCRIPTION:                                              | 
#  Simple modification of blast conversion program.
#                                                           |
# USAGE:                                                    |
#  cnv_blast2gff.pl -i infile.bln -o blast_out.gff          |
#                                                           |
# VERSION:                                                  |
#  $Rev: 923 $                                              |
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
use Bio::SearchIO;             # Parse BLAST output
# The following needed for printing help
use Pod::Select;               # Print subsections of POD documentation
use Pod::Text;                 # Print POD doc as formatted text file
use IO::Scalar;                # For print_help subfunction
use IO::Pipe;                  # Pipe for STDIN, STDOUT for POD docs
use File::Spec;                # Convert a relative path to an abosolute path

#-----------------------------+
# PROGRAM VARIABLES           |
#-----------------------------+
my ($VERSION) = q$Rev: 923 $ =~ /(\d+)/;
# Get GFF version from environment, GFF2 is DEFAULT
my $gff_ver = uc($ENV{DP_GFF}) || "GFF2";

#-----------------------------+
# LOCAL VARIABLES             |
#-----------------------------+
my $max_e;
my $min_len;

# Set variable scope
#my $seqname = "seq";
my $qry_name;
my $infile;                     # Default to STDIN if not given
my $outfile;                    # Default to STDOUT if not given
my $blast_program;
my $blast_alignment = 0;        # Alternatives include 8 and 9
my $param;
my $feature_type = "match";

# Booleans
my $simple = 1; # Doing something simple while SearchIO not working
my $verbose = 0;
my $help = 0;
my $quiet = 0;
my $show_help = 0;
my $show_usage = 0;
my $show_man = 0;
my $show_version = 0;
my $do_append = 0;
# Min t length from BLAT result
my $min_t_length;  
my $num_tgaps_allowed;
my $no_tgaps_allowed = 0;
my $min_tratio;             # Minimum transcript ratio
my $min_match_length;

#-----------------------------+
# COMMAND LINE OPTIONS        |
#-----------------------------+
my $ok = GetOptions(# REQUIRED
		    "i|infile=s"     => \$infile,
                    "o|outfile=s"    => \$outfile,
		    "n|s|seqname|name=s"  => \$qry_name,
		    # OPTIONS
		    "simple"         => \$simple,
		    "min-match-len=i"=> \$min_match_length,
		    "num-tgaps=f"    => \$num_tgaps_allowed,
		    "min-tratio=f"   => \$min_tratio,
		    "no-tgaps"       => \$no_tgaps_allowed,
		    "min-t-length=s" => \$min_t_length,
		    "gff-ver=s"      => \$gff_ver,
		    "f|feature=s"    => \$feature_type,
		    "p|program=s"    => \$blast_program,
		    "m|align=s"      => \$blast_alignment,
		    "d|database=s"   => \$param,
		    "e|maxe=s"       => \$max_e,
		    "l|minlen=s"     => \$min_len,
		    # Booleans
		    "verbose"        => \$verbose,
		    "append"         => \$do_append,
		    "usage"          => \$show_usage,
		    "version"        => \$show_version,
		    "man"            => \$show_man,
		    "h|help"         => \$show_help,
		    "q|quiet"        => \$quiet,);
 
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
# ALLOWING FOR ZERO GAPS      |
#-----------------------------+
if ($num_tgaps_allowed) {
    print STDERR "Evaluating gaps\n";
    print STDERR "$num_tgaps_allowed\n";
    $num_tgaps_allowed = int($num_tgaps_allowed);
    if ($num_tgaps_allowed < 1 ) {
	print STDERR "Setting zero man\n";
	$num_tgaps_allowed = "0 but true";
    }
    else {
	$num_tgaps_allowed = int($num_tgaps_allowed);
	print STDERR "Number transcript gaps is more than zero\n";
    }
}

# THE NO GAPS OPTION WILL OVERRIDE NUMBER OF GAPS VAR
#if ($no_tgaps_allowed) {
#    $num_tgaps_allowed = "0 but true";
#
#}

# CHECK HOW NUMBER OF GAPS WAS PASSED
#print STDERR "The number of gaps is $num_tgaps_allowed\n";
#if ($num_tgaps_allowed <1) {
#    print STDERR "This evaluates to less than one\n";
#}
#else {
#    print STDERR "This evalulates as greater than or equal to one\n";
#}
#exit;

#-----------------------------+
# SHOW REQUESTED HELP         |
#-----------------------------+
if ( $show_usage ) {
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
    print "\ncnv_blast2gff.pl:\n".
	"Version: $VERSION\n\n";
    exit;
}


#-----------------------------------------------------------+
# MAIN PROGRAM BODY                                         |
#-----------------------------------------------------------+
if ($simple) {
    simple_blat2gff ($infile, $outfile, $do_append, $qry_name, $blast_alignment,
	       $param, $blast_program, $feature_type);
} 
else {
    blast2gff ($infile, $outfile, $do_append, $qry_name, $blast_alignment,
	       $param, $blast_program, $feature_type);
}

exit;


#-----------------------------------------------------------+
# SUBFUNCTIONS                                              |
#-----------------------------------------------------------+

sub simple_blat2gff {
    print STDERR "Doing simple\n" if $verbose;

    my ($blastin, $gffout, $append, $seqname, $align, $suffix, $prog,
	$feature) = @_;

    # Assume BLAT report does not have header
    my $has_header = 0;
    my $sum_tblock_length;

    #-----------------------------+
    # FILE HANDLES                |
    #-----------------------------+
    open (BLATIN, $blastin) ||
	die "Con not open input file $blastin\n";

    if ($gffout) {
	open (GFFOUT, $gffout) ||
	    die "Can not open $gffout for output\n";
    }
    else {
    	open (GFFOUT, ">&STDOUT") ||
	    die "Can not open STDOUT for writing\n";
    }

    # TO FIX LATER
    my $source = "BLAT";
    my $line_count = 0;


    while (<BLATIN>) {

	$line_count ++;

	#///
#
#        CURRENT DEFAULT IS TO PRINT ALL LINES
#
	#///
	my $do_print = 1;

	chomp;

	my @blat_parts = split( /\t/, $_);
	my $num_blat_parts = @blat_parts;


	#-----------------------------+
	# CODE TO SKIP HEADER LINES   |
	#-----------------------------+
	if ($blat_parts[0]) {
	    if ($blat_parts[0] =~ "psLayout") {
		$has_header = 1;
		print STDERR "It has a header\n" if $verbose;
	    }
	}

	if ($has_header) {
	    if ($line_count < 6) {
		next;
	    }
	}

	#-----------------------------+
	# GETTING VARS FROM INPUT     |
	#-----------------------------+
	my $gff_seqname;
	if ($seqname) {
	    $gff_seqname = $seqname;
	}
	else {
	    $gff_seqname = $blat_parts[9];
	}

	my $match_length = $blat_parts[0];
	my $strand = $blat_parts[8];
	my $t_gap_count = int($blat_parts[6]);
	my $start = $blat_parts[11];
	my $end = $blat_parts[12];
	my $attribute = $blat_parts[13];
	my $t_name = $blat_parts[13];
	# Assume that BLAT results start at line six
	my $t_size = $blat_parts[14];
	my $t_start = $blat_parts[15];
	my $t_end = $blat_parts[16];
        my $t_match_len = $t_end - $t_start;
	my $block_count = $blat_parts[17];
	my $score = "1";                               # FAKE SCORE

	# GET BLOCK INFORMATION
	my $block_sizes = $blat_parts[18];
	my @block_size_parts = split (/\,/, $block_sizes);
	my $num_blocks = @block_size_parts;

	$sum_tblock_length = 0;
	if ($num_blocks ) {
	    $sum_tblock_length += $_ for @block_size_parts;
#	    $sum += $_ for @a;
	}

	my $tblock_length;
#	$tblock_length = $sum_tblock_length - $num_blocks + 1;
#	print STDERR $_."\n";
#	print STDERR $num_blocks.":".
#	    $tblock_length.":".
#	    $sum_tblock_length."\n\n";

	my $match_ratio = $match_length/$t_size;

	#-----------------------------------------------------------+
	# FILTERING RESULTS                                         |
	#-----------------------------------------------------------+

	#-----------------------------+
	# FILTER BY HIT LENGTH        |
	#-----------------------------+
	if ($min_len) {
	    print STDERR "You have a minimun length\n";
	}

	#-----------------------------+
	# FILTER BY NUMBER OF GAPS    |
	#-----------------------------+
	if ($num_tgaps_allowed) {
	    if ($t_gap_count > $num_tgaps_allowed) {
		$do_print = 0;
	    }
	}

	#-----------------------------+
	# FILTER BY TRANSCRIPT LENGTH |
	#-----------------------------+
	# The minimum length of the transcript that is an acceptable match
	if ($min_t_length) {
	    if ( $t_size <= $min_t_length) {
		$do_print = 0;
	    }
	}

	# FILTER BY MINIMUM TOTAL MATCH LENGTH
	if ($min_match_length) {
	    if ($match_length < $min_match_length) {
		$do_print = 0;
	    }

	}
	
	#-----------------------------+
	# FILTER BY MATCH RATIO       |
	#-----------------------------+
	if ($min_tratio) {
	    if ( $match_ratio >= $min_tratio) {
		#$min_tratio;
	    }
	}

	#-----------------------------+
	# SETTING ATTRIBUTE LINE      |
	#-----------------------------+
	if ($gff_ver =~ "GFF3") {
	    $feature = "match_part"; 
	    $attribute = "ID=".$source."_".$blat_parts[13].
		"; ".
		"Name=".$blat_parts[13]."; ".
		"Target=".$blat_parts[13]." ".
		$start." ".$end;
	}

	#-----------------------------+
	# PRINT OUTPUT TO GFF FILE    |
	#-----------------------------+
	if ($do_print) {
	    print GFFOUT "$gff_seqname\t".           # Seqname
		"$source\t".                     # Source
		"$feature\t".                    # Feature type name
		"$start\t".                      # Start
		"$end\t".                        # End
		"$score\t".                      # Score
		"$strand\t".                     # Strand
		"."."\t".                        # Frame
#		$line_count.                     # Attribute
		$attribute.                      # Attribute
		"\n";                            # newline

	    # TEST WORKING OUT LENGTH OF MATCHES
	    print STDERR $_."\n";
	    print STDERR "\t".$num_blocks.":".
		$sum_tblock_length."\n";
	    print STDERR "\tMatch Ratio: ".$match_ratio."\n";
	    print STDERR "\n";
	}

#	if ($line_count == 100) {
#	    exit;
#	}

    } # End of while BLATIN

}

sub blast2gff {
# CONVERT BLAST TO GFF 
    
    # seqname - ID to use in the seqname field
    # blastin - path to the blast input file
    # gffout  - path to the gff output file
    # append  - boolean append data to existing file at gff out
    # align   - the alignment type of the blast output
    # suffix  - the suffix to add to type of blast
    # prog    - blast program (blastn, blastx, wublast etc)
    # gff second column is built as progparam
    my ($blastin, $gffout, $append, $seqname, $align, $suffix, $prog,
	$feature) = @_;
    my $blastprog;        # Name of the blast program used (blastn, blastx)
    my $dbname;           # Name of the database blasted
    my $hitname;          # Name of the hit
    my $start;            # Start of the feature
    my $end;              # End of the feature
    my $strand;           # Strand of the hit

    # Recode from number for tab delimited blast
    if ($align == 8 || $align ==9) {
	$align = "tab";
    }
    
    #-----------------------------+
    # BLAST INTPUT OBJECT         |
    #-----------------------------+
    # This allows for input from NCBI-BLAST, and WUBLAST in either 
    # default or tab delimited output as well as 
    # INPUT FROM FILE PATH
    my $blast_report;
    if ($blastin) {
	
	# TAB ALIGNED BLAST OUTPUT
	if ($align =~ "tab" ) {
	    print STDERR "Parsing Tab BLAST\n" if $verbose;
	    $blast_report = new Bio::SearchIO (
#	    $blast_report = new Bio::SearchIO::blasttable ( 
		'-format' => 'blasttable',
		'-file'   => $blastin )
		|| die "Could not open BLAST input file:\n$blastin.\n";
	}
	else {
	    print STDERR "Parsing BLAT\n" if $verbose;
	    $blast_report = new Bio::SearchIO ( '-format' => 'blast',
						'-file'   => $blastin,
						'-signif' => $max_e,
						'-min_query_len' => $min_len) 
		|| die "Could not open BLAST input file:\n$blastin.\n";
	} # Default to default blast alignment
    }

    # INPUT FROM STDIN
    else {
	if ($align =~ "tab") {
	    $blast_report = new Bio::SearchIO::blasttable ( 
		'-format' => 'blasttable',
		'-fh'     => \*STDIN )
		|| die "Could not open STDIN for blast table input\n";
	}
	else {
	    print STDERR "Expecting input from STDIN\n";
	    $blast_report = new Bio::SearchIO ( '-format' => 'blast',
						'-fh'     => \*STDIN,
						'-signif' => $max_e,
						'-min_query_len' => $min_len) 
		|| die "Could not open STDIN for blast input\n";
	}
    }
    
    #-----------------------------+
    # FILEHANDLE FOR OUTPUT       |
    #-----------------------------+
    if ($gffout) {
	if ($append) {
	    open (GFFOUT, ">>$gffout") 
		|| die "Can not open file:\n $gffout\n";
	}
	else {
	    open (GFFOUT, ">$gffout") 
	    || die "Can not open file:\n $gffout\n";
    
	    if ($gff_ver =~ "GFF3") {
		print GFFOUT "##gff-version 3\n";
	    }

	}
    }
    else {
	open (GFFOUT, ">&STDOUT") ||
	    die "Can not print to STDOUT\n";
	if ($gff_ver =~ "GFF3") {
	    print GFFOUT "##gff-version 3\n";
	}
    }
    
    while (my $blast_result = $blast_report->next_result())
    {

	# Not using the following for the moment
	# 01/30/2009
	#$blastprog = $blast_result->algorithm;
	$blastprog = $blast_result->algorithm;
	#$dbname = $blast_result->database_name();
	
	#-----------------------------+
	# GET SEQNAME IF NOT PROVIDED |
	#-----------------------------+
	# First attempt to extract from blast report
	# otherwise use 'seq' as the indentifier.
	# For multiple query sequence, this will use
	# the same identifier for all
	unless ($seqname) {
	    if ($blast_result->query_name) {
		$seqname = $blast_result->query_name();
	    }
	    else {
		$seqname = "seq";
	    }
	}
	# Remove prohibited characters
	if ($gff_ver =~ "GFF3") {
	    $seqname = seqid_encode($seqname);
	}

	#-----------------------------+
	# GET ALGORITHM IF UNKNOWN    |
	#-----------------------------+
	# This attempts to get the algorithm from the
	# blast report, otherwise the generic 'blast' is used
	unless ($prog) {
	    if ($blast_result->algorithm) {
		$prog = $blast_result->algorithm;
	    }
	    else {
		$prog = "blast";
	    }
	}
	# Sanitize program name
	if ($gff_ver =~ "GFF3") {
	    $prog = gff3_encode($prog);		
	}
	
        #-----------------------------+
	# APPEND DB IDENTIFIER        |
	#-----------------------------+
	# Providing a suffix at the command line will override
	# the name provided in the blast report
	my $source = $prog;
	my $blast_db;
	if ($suffix) {
	    if ($gff_ver =~ "GFF3") {
		$suffix = gff3_encode($suffix);
	    }
	    $source = $source.":".$suffix;
	}
	elsif ($blast_result->database_name()) {
	    # Remove trailing white space from db name
	    $blast_db = $blast_result->database_name();
	    # Trim trailing white space
	    $blast_db =~ s/\s+$//;
	    if ($gff_ver =~ "GFF3") {
		$blast_db = gff3_encode($blast_db);
	    }
	    $source = $source.":".$blast_db;
	}

    	while (my $blast_hit = $blast_result->next_hit()) {
	    print STDERR "hit\n";
	    my $hitname = $blast_hit->name();	    
	    if ($gff_ver =~ "GFF3") {
		$hitname = gff3_encode($hitname);
	    }

	    my $hsp_num = 0;

	    while (my $blast_hsp = $blast_hit->next_hsp())
	    {

		$hsp_num++;
		my $hsp_id = sprintf("%04d", $hsp_num);

		#-----------------------------+
		# GET STRAND                  |
		#-----------------------------+
		# NOTE: This defaults to positive strand!!
		$strand = $blast_hsp->strand('query');
		if ($strand =~ "-1") {
		    $strand = "-";
		}
		elsif ($strand =~ "1") {
		    $strand = "+";
		}
		else {
		    $strand = "+";
		}
		
		#-----------------------------+
		# GET QRY START AND END       |
		#-----------------------------+
		# Make certain that start coordinate is
		# less then the end coordinate
		if ( $blast_hsp->start() < $blast_hsp->end() ) {
		    $start = $blast_hsp->start();
		    $end = $blast_hsp->end();
		}
		else {
		    $start = $blast_hsp->end();
		    $end = $blast_hsp->start();
		}

		#-----------------------------+
		# GET HIT START AND END       |
		#-----------------------------+
		my $hit_start;
		my $hit_end;
		if ( $blast_hsp->start('hit') < $blast_hsp->end('hit') ) {
		    $hit_start = $blast_hsp->start('hit');
		    $hit_end = $blast_hsp->end('hit');
		}
		else {
		    $hit_start = $blast_hsp->end('hit');
		    $hit_end = $blast_hsp->start('hit');
		}


		#-----------------------------+
		# Set attribute               |
		#-----------------------------+
		my $attribute;
		if ($gff_ver =~ "GFF3") {
		    $feature = "match_part"; 
		    $attribute = "ID=".$source."_".$hitname.
			";".
			"Name=".$hitname.";".
			"Target=".$hitname." ".
			$hit_start." ".$hit_end;
		}
		else {
		    $attribute = $hitname;
		}
		
		#-----------------------------+
		# PRINT OUTPUT TO GFF FILE    |
		#-----------------------------+

		print GFFOUT "$seqname\t".           # Seqname
		    "$source\t".                     # Source
		    "$feature\t".                    # Feature type name
		    "$start\t".                      # Start
		    "$end\t".                        # End
		    $blast_hsp->score()."\t".        # Score
		    "$strand\t".                     # Strand
		    ".\t".                           # Frame
		    $attribute.                      # Attribute
		    "\n";                            # newline

	    } # End of while next hsp
	} # End of while next hit
    } # End of while next result
    
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

sub seqid_encode {
    # Following conventions for GFF3 v given at http://gmod.org/wiki/GFF3
    # Modified from code for urlencode in the perl cookbook
    # Ids must not contain unescaped white space, so spaces are not allowed
    my ($value) = @_;
    $value =~ s/([^[a-zA-Z0-9.:^*$@!+_?-|])/"%" . uc(sprintf "%lx" , unpack("C", $1))/eg;
    return ($value);
}

sub gff3_encode {
    # spaces are allowed in attribute, but tabs must be escaped
    my ($value) = @_;
    $value =~ s/([^[a-zA-Z0-9.:^*$@!+_?-| ])/"%" . uc(sprintf "%lx" , unpack("C", $1))/eg;
    return ($value);
}

=head1 NAME

cnv_blast2gff.pl - Convert blast output to GFF

=head1 VERSION

This documentation refers to program version $Rev: 923 $

=head1 SYNOPSIS

=head2 Usage

    cnv_blast2gff.pl -i blast_report.bln -o blast_out.gff -n seq_id

=head2 Required Arguments

    -i      # Path to the input file
            # If not specified the program will expect input from STDIN
    -o      # Path to the output file
            # If not specified the program will write to STDOUT
    -s      # Name of the query sequence
            # If a name is not specified, the name will be
            # extracted from the blast report or use "seq"

=head1 DESCRIPTION

This program will translate a blast report for a single query sequence into
the GFF format. Since this uses the general bioperl BLAST parser code,
this should also be able to parse output from BLAT or wublast. This code
works best for converting a blast report of a single query sequence against
a single database.

=head1 REQUIRED ARGUMENTS

=over 2

=item -i,--infile

Path of the input file. If an input file is not specified, the program
will expect input from STDIN.

=item -o,--outfile

Path of the output file. If an output file is not specified, the program
will write output to STDOUT.

=back

=head1 OPTIONS

=over 2

=item --gff-ver

The GFF version for the output. This will accept either gff2 or gff3 as the
options. By default the GFF version will be GFF2 unless specified otherwise.
The default GFF version for output can also be set in the user environment
with the DP_GFF option. The command line option will always override the option
defined in the user environment. 

=item -s,--seqname

Identifier of the sequence that has been used as the query sequence
in the blast report. This will be used as the first column in the 
gff output file.

=item -p,--program

The blast program used. This will be used to identify the source program in the
second column of the GFF output file. Example of valid values include
blastn, blastx, or wublast.

=item --feature

The type of feature. Be default, this is set to match_part in order to
follow the Sequence Ontology guidelines and to facilitate visualizing
this blast report in Apollo and Gbrowse. It is also possible to set this to an
ontology complient name such as match or expressed_sequence_match.

=item -d,--database

The name for the database that was blasted against. If provided, this
will be appended to the program variable in the second colum of the 
GFF output file.

=item -m,--align

The alignment format use in the BLAST report to be parsed.
The program will assume that you are using the default alginment format for
blast. Otherwise, you can specify 'tab' or '8' or '9' for tab delimited blast.

=item -e,--maxe

The maximum e value threshold to accept.

=item -l,--min-len

The minimum length to accept.

=item --verbose

Run the program with maximum reporting of error and status messages.

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

The following are examples on how to use the cnv_blast2gff.pl program.

=head2 Typical Use

The typical use of this program will be to convert an existing blast output
to to the GFF file format:

  cnv_blast2gff.pl -i blast_result.bln -o parsed_result.gff

This will generate a GFF format file named parsed_result.gff.

=head2 Piping BLAST Result Directly to the Conversion Utility

It is also possible to directly send the blast result to the cnv_blast2gff.pl
program using the standard streams.

  blastall -p blastin .. | cnv_blast2gff.pl -o blast_result.gff

This will take the blast output from NCBI's blastall program and 
convert the output to the gff format.

=head2 Combining Blast Results in GFF format

Since the cnv_blast2gff.pl program will write the results to the
standard output stream if no file path is specified, it is possible to
use standard unix commands to combined results. Consider the following
set of commands:

  cnv_blast2gff.pl blast_result01.bln > combined_results.gff
  cnv_blast2gff.pl blast_result02.bln >> combined_results.gff
  cnv_blast2gff.pl blast_result03.bln >> combined_results.gff

This will combined the blast results from the 01, 02 and 03 search into
a single gff file named combined_results.gff. 
Alternatively if the blast results are in the default format 
(not tab delimited -m8/m9), this can also be done by taking advantage
of the STDIN stream:

  cat *.bln | cnv_blast2gff.pl > combined_results.gff

=head2 Specify the Sequence ID with --name

The first column in the GFF output file results indicates the id of the
sequence that is being annotated. By default, the cnv_blast2gff.pl program
will attempt to extract this ID from the blast result. It is also possible
to specify this from the command line using the --name option. For example
consider you had a blast report that gave the following result"

  cnv_blast2gff.pl -i bl_result.bln -o gff_result.gff

that generated a gff file like the following

 HEX3045G05   blast:mips   exon     8537    8667    39	 +    .	 rire1
 HEX3045G05   blast:mips   exon     9911    9996    38	 +    .	 rire1
 HEX3045G05   blast:mips   exon     10025   10191   36	 +    .	 rire1
 HEX3045G05   blast:mips   exon     76161   76235   35	 +    .	 rire1
 HEX3045G05   blast:mips   exon     81151   81200   34	 +    .	 rire1
 ...

where HEX304GO5 indicates the sequence id. This sequence identifier could
be modified using the --name option:

  cnv_blast2gff.pl -i bl_result.bln -o gff_result.gff --seqname HEX001

this would give the following result:

 HEX001   blast:mips   exon     8537    8667    39	 +    .	 rire1
 HEX001   blast:mips   exon     9911    9996    38	 +    .	 rire1
 HEX001   blast:mips   exon     10025   10191   36	 +    .	 rire1
 HEX001   blast:mips   exon     76161   76235   35	 +    .	 rire1
 HEX001   blast:mips   exon     81151   81200   34	 +    .	 rire1
 ...

=head2 Specify the Database name with --database

By default the cnv_blast2gff.pl program will identify the database in the
second column as a suffix to the blast program, separated by a colon.
The command:

  cnv_blast2gff.pl -i bl_result.bln -o gff_result.gff

that generated a gff file like the following

 HEX3045G05   blast:mips   exon     8537    8667    39	 +    .	 rire1
 HEX3045G05   blast:mips   exon     9911    9996    38	 +    .	 rire1
 HEX3045G05   blast:mips   exon     10025   10191   36	 +    .	 rire1
 HEX3045G05   blast:mips   exon     76161   76235   35	 +    .	 rire1
 HEX3045G05   blast:mips   exon     81151   81200   34	 +    .	 rire1
 ...

Could have the database suffix modified using the --database option as follows:

  cnv_blast2gff.pl -i bl_result.bln -o gff_result.gff --name tes

This would modify the gff output to the following:

 HEX3045G05   blast:tes   exon     8537    8667    39	 +    .	 rire1
 HEX3045G05   blast:tes   exon     9911    9996    38	 +    .	 rire1
 HEX3045G05   blast:tes   exon     10025   10191   36	 +    .	 rire1
 HEX3045G05   blast:tes   exon     76161   76235   35	 +    .	 rire1
 HEX3045G05   blast:tes   exon     81151   81200   34	 +    .	 rire1

=head1 DIAGNOSTICS

=over 2

=item * Expecting input from STDIN

If you see this message, it may indicate that you did not properly specify
the input sequence with -i or --infile flag. 

=back

=head1 CONFIGURATION AND ENVIRONMENT

=over 2

=item DP_GFF

The DP_GFF variable can be defined in the user environment to set
the default GFF version output. Valid settings are 'gff2' or
'gff3'.

=back

=head1 DEPENDENCIES

=head2 Required Software

This program requires output from a BLAST program. Since this program makes
use of the BioPerl blast parser, it should be possible to convert
local alignment results from any of the following programs:

=over 2

=item * BLAST

=item * PSIBLAST

=item * PSITBLASTN

=item * RPSBLAST

=item * WUBLAST

=item * bl2seq

=item * WU-BLAST

=item * BLASTZ

=item * BLAT

=item * Paracel 

=item * BTK

=back

=head2 Required Perl Modules

The following perl modules are required for this program:

=over 2

=item * Bio::SearchIO

The SearchIO module is part of the bioperl module.

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

=item * Limited BLAST testing

This program has only been tested with BLAST output from the NCBI-BLAST package
using standard and tab delimited output. If you find that there are limiations
to this program that limit your use, please email the author or file a 
bug report on the Sourceforge website:
http://sourceforge.net/tracker/?group_id=204962

=back

=head1 REFERENCE

Your use of the DAWGPAWS programs should reference the following manuscript:

JC Estill and JL Bennetzen. 2009. 
"The DAWGPAWS Pipeline for the Annotation of Genes and Transposable 
Elements in Plant Genomes." Plant Methods. 5:8.

=head1 LICENSE

GNU General Public License, Version 3

L<http://www.gnu.org/licenses/gpl.html>

THIS SOFTWARE COMES AS IS, WITHOUT ANY EXPRESS OR IMPLIED
WARRANTY. USE AT YOUR OWN RISK.

=head1 AUTHOR

James C. Estill E<lt>JamesEstill at gmail.comE<gt>

=head1 HISTORY

STARTED: 08/06/2007

UPDATED: 01/12/2010

VERSION: $Rev: 923 $

=cut

#-----------------------------------------------------------+
# HISTORY                                                   |
#-----------------------------------------------------------+
#
# 08/06/2007
# -Program started
# -Main body of program written
# 01/30/2009
# -Updated POD documentation
# -Moved POD documentation to end of the program
# -Updated to new print_help subfunction that extracts
#  help from POD
# -Added default to print to STDOUT if not outfile
#  path given at command line
# -Added default input from STDIN in no infile path
#  given at command line
# -Added package name DAWGPAWS
# -Added options to pass database name, query sequence id,
#  and program used. Otherwise these are attempted to
#  be extracted from the command line
# 01/31/2009
# -Fixed print_help to extract help from POD
# 03/09/2009
# -Added some examples
# 03/30/2009
# -added -s and --seqname
# -added support for --feature, kept exon as default
# 01/12/2010
# - Added support for GFF3 output format
# - Updated POD code to include changes
# - Still need to add code to sanitize hit names for
#   illegal characters
# 01/13/2010
# - Added encode subfunctions to sanitize names
#    *seqid_encode
#    *attribute_encode
# - Moved some variable definitions to within the 
#   'for each blast_result' loop to speed things up
# - For GFF3 set 3rd col to match instead of exon
# 01/14/2010
# - Modified match to match_part. This will allow parts
#   to be joined as single molecules in Apollo
# - Updated examples to show using cat to combine results
#   for conversion using the STDIN stream.
# - Fixed bug where $prog was not being set to the 
#   BLAST algorithm.

