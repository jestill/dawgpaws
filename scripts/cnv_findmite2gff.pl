#!/usr/bin/perl -w
#-----------------------------------------------------------+
#                                                           |
# cnv_findmite2gff.pl - Convert findmite results to gff     |
#                                                           |
#-----------------------------------------------------------+
#                                                           |
#  AUTHOR: James C. Estill                                  |
# CONTACT: JamesEstill_@_gmail.com                          |
# STARTED: 08/31/2007                                       |
# UPDATED: 08/31/2007                                       |
#                                                           |
# DESCRIPTION:                                              |
#  Parses findmite results to a gff file. Assigns unique
#  name to all of the results
#
#-----------------------------------------------------------+
#
# TO DO: 
#  -Create output file of all putative mites
#  -Add subfuncton to the batch_findmite program

=head1 NAME

cnv_findmite2gff.pl - Convert FINDMITE output to gff format

=head1 VERSION

This documentation refers to program $Rev$

=head1 SYNOPSIS

  USAGE:
    cnv_findmite2gff -i InFile -o OutFile
    
=head1 DESCRIPTION

Converts FINDMITE output to GFF format.

=head1 COMMAND LINE ARGUMENTS

=head2 Required Arguments

=over 2

=item -i,--infile

Path of the input file.

=item -o,--outfile

Path of the output file.

=back


=head2 Additional Options

=over 2

=item -q,--quiet

Run the program with minimal output.

=back

=head2 Additional Information

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

package DAWGPAWS;

#-----------------------------+
# INCLUDES                    |
#-----------------------------+
use strict;
use Getopt::Long;

#-----------------------------+
# PROGRAM VARIABLES           |
#-----------------------------+
my ($VERSION) = q$Rev$ =~ /(\d+)/;

#-----------------------------+
# VARIABLE SCOPE              |
#-----------------------------+
my $infile;                    # Findmite outupt file to convert
my $outfile;                   # The gff output file produced

# Booleans
my $help = 0;
my $quiet = 0;
my $show_help = 0;
my $show_usage = 0;
my $show_man = 0;
my $show_version = 0;
my $verbose = 0;


#-----------------------------+
# COMMAND LINE OPTIONS        |
#-----------------------------+
my $ok = GetOptions(# Required arguments
		    "i|infile=s"  => \$infile,
                    "o|outfile=s" => \$outfile,
		    # Additional options
		    "q|quiet"     => \$quiet,
		    "verbose"     => \$verbose,
		    # Additional information
		    "usage"       => \$show_usage,
		    "version"     => \$show_version,
		    "man"         => \$show_man,
		    "h|help"      => \$show_help,);

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

#-----------------------------------------------------------+
# MAIN PROGRAM BODY                                         |
#-----------------------------------------------------------+
print "test\n";

findmite2gff ($infile,$outfile, 0, "HEXTEST");

exit;

#-----------------------------------------------------------+ 
# SUBFUNCTIONS                                              |
#-----------------------------------------------------------+

sub findmite2gff {
    
    #-----------------------------+
    # SUBFUNCTION VARS            |
    #-----------------------------+
    my ($findmite_in, $gff_out, $append, $seqname) = @_;
    
    # GFF OUT VARS
    my $gff_seq_id;                 # Seq id for use in gff
    my $gff_start;                  # Start of the feature
    my $gff_end;                    # End of the feature
    my $gff_source;                 # Source ie findmite_at_11
    my $gff_name;                   # Name of the feature 

    # FINDMITE HEADER VARS
    my $fm_dir_repeats;             # Direct repeats
    my $fm_tir_len;                 # Lenght of the TIR
    my $fm_mismatch_max;            # Max number of mistmatches
    my $fm_filter_at;               # Were A/T strings filtered 
    my $fm_filter_cg;               # Were C/G strings filtered
    my $fm_filter_atta;             # Were AT/TA strings filtered
    my $fm_filter_2_base;           # Percent (0 to 100)
    my $fm_min_dist;                # Minimum distance used
    my $fm_max_dist;                # Maximum distance used
    my $fm_file_an;                 # Name of the file analyzed

    # FINDMATE INDIVIDUAL MITE VARS
    my $fm_pattern;                 # Pattern
    my $fm_seq_id;                  # Id of the query sequence
    my $fm_num_mismatch;            # Num of mismatches
    my $fm_seq = "";                # Sequence string as parsed from findmite
    my $mite_seq_string;            # Sequence string of the putatitve mite
    my $mite_context_5;             # 5' Context of the mite (40bp)
    my $mite_context_3;             # 3' Context of the mite (40bp)
    my $mite_start;                 # Start of the mite as parsed from FM
    my $mite_end;                   # End of the mite as parsed from FM

    # Counter
    my $mite_num = 0;               # Incremented ID number for the mite

    # BOOLEANS
    my $in_seq = 0;                 # Boolean, in seq data (past header info)

    #-----------------------------+
    # OPEN FILES                  |
    #-----------------------------+

    open (INFILE, "<$findmite_in") ||
	die "Can not open input file:\n$findmite_in\n";

    open (GFFOUT, ">$gff_out") ||
	die "Can not open output file:\n$gff_out\n";

    #-----------------------------+
    # PROCESS FINDMITE FILE       |
    #-----------------------------+
    while (<INFILE>) {
	chomp;
 	#print $_."\n";

	if(m/^Pattern  : (.*)/) {

	    # If we have previously loaded seq data then
	    # futher parse the sequence string and and
	    # print the results to the gff output
	    if ($in_seq) {
		                   
		if ($fm_seq =~ m/(.*)\((.*)\)\-{5}(.*)\-{5}\((.*)\)(.*)/) {
		    $mite_context_5 = $1;
		    $mite_start = $2;
		    $mite_seq_string = $3;
		    $mite_end = $4;
		    $mite_context_3 = $5;

		    print STDERR "\t5CON: $mite_context_5\n";
		    print STDERR "\tSTAR: $mite_start\n";
		    print STDERR "\tMITE: $mite_seq_string\n";
		    print STDERR "\tEND : $mite_end\n";
		    print STDERR "\t3CON: $mite_context_3\n";
		    print STDERR "\n\n";

		    #-----------------------------+
		    # PRINT TO GFF FILE           |
		    #-----------------------------+
		    # Parse seq id for shorter name
		    if ($fm_seq_id =~ m/^>(\S*)\s./) {
			$gff_seq_id = $1;
		    }
		    elsif ($fm_seq_id =~ m/^>(.*)/) {
			$gff_seq_id = $1;
		    }
		    else {
			$gff_seq_id = $fm_seq_id;
		    }

		    $gff_source = "findmite:".
			$fm_dir_repeats."_".$fm_tir_len;
		    $gff_name = $gff_seq_id."_".$fm_dir_repeats.
			"_".$fm_tir_len."_".$mite_num;
		    $gff_start = $mite_start;
		    $gff_end = $mite_end;

		    print GFFOUT 
			"$gff_seq_id\t".            # Seq name
			"$gff_source\t".            # Source
			"mite\t".                   # Feature type
			"$gff_start\t".             # Start
			"$gff_end\t".               # End
			".\t".                      # Score
			".\t".                      # Strand
			".\t".                      # Frame
			"$gff_name\n";              # Feature name

		}

		# Reset vals to null
		$fm_seq = "";
	    }

	    $in_seq = 1;
	    $mite_num++;
	    $fm_pattern = $1;
	    print STDERR "$fm_pattern\n";
	}
	elsif(m/^Sequence : (.*)/) {
	    $fm_seq_id = $1;
	    print STDERR "\t$fm_seq_id\n";
	}
	elsif(m/^Mismatch : (.*)/){
	    $fm_num_mismatch = $1;
	    print STDERR "\t$fm_num_mismatch\n";
	}
	elsif($in_seq) {
	    $fm_seq = $fm_seq.$_;
	}
	#-----------------------------+
	# HEADER INFORMATION          | 
	#-----------------------------+
	elsif(m/Direct repeats.{17}(.*)/) {
	    $fm_dir_repeats = $1;
	}
	elsif(m/Length of TIR.{18}(.*)/) {
	    $fm_tir_len = $1;
	}
	elsif(m/Number of mis.{18}(.*)/) {
	    $fm_num_mismatch = $1;
	}
	elsif(m/Filtering A\/T.{18}(.*)/) {
	    $fm_filter_at = $1;
	}
	elsif(m/Filtering C\/G.{18}(.*)/) {
	    $fm_filter_cg = $1;
	}
	elsif(m/Filtering AT\/.{18}(.*)/) {
	    $fm_filter_atta = $1;
	}
	elsif(m/Filtering 2.{20}(.*)/) {
	    $fm_filter_2_base = $1;
	}
	elsif(m/Minimum dist.{19}(.*)/) {
	    $fm_min_dist = $1;
	}
	elsif(m/Maximum dist.{19}(.*)/) {
	    $fm_max_dist = $1;
	}
	elsif(m/The results from the input.{17}(.*)/) {
	    $fm_file_an = $1;
	}

    }

    
    #-----------------------------+
    # PRINT OUT OF DATA IN VARS   |
    #-----------------------------+
    if ($fm_seq =~ m/(.*)\((.*)\)\-\-\-\-\-(.*)\-\-\-\-\-\((.*)\)(.*)/) {
	$mite_context_5 = $1;
	$mite_start = $2;
	$mite_seq_string = $3;
	$mite_end = $4;
	$mite_context_3 = $5;
	
	print STDERR "\t5CON: $mite_context_5\n";
	print STDERR "\tSTAR: $mite_start\n";
	print STDERR "\tMITE: $mite_seq_string\n";
	print STDERR "\tEND : $mite_end\n";
	print STDERR "\t3CON: $mite_context_3\n";
	print STDERR "\n\n";
	
    }
    # Reset vals to null
    $fm_seq = "";


    # Print while debu
    print STDERR "FILE  : $fm_file_an\n";
    print STDERR "DIRREP: $fm_dir_repeats\n";
    print STDERR "TIRLEN: $fm_tir_len\n";
    print STDERR "NUMMIS: $fm_num_mismatch\n";
    print STDERR "F_AT  : $fm_filter_at\n";
    print STDERR "F_CG  : $fm_filter_cg\n";
    print STDERR "F_ATTA: $fm_filter_atta\n";
    print STDERR "F_2BAS: $fm_filter_2_base\n";
    print STDERR "MIN   : $fm_min_dist\n";
    print STDERR "MAX   : $fm_max_dist\n";

    close INFILE;
    close GFFOUT;

}

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


=head1 HISTORY

STARTED: 08/30/2007

UPDATED: 08/30/2007

VERSION: $Rev$

=cut

#-----------------------------------------------------------+
# HISTORY                                                   |
#-----------------------------------------------------------+
# 
# 08/30/2007
# - Program started
