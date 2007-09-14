#!/usr/bin/perl -w
#-----------------------------------------------------------+
#                                                           |
# cnv_ltrfinder2gff.pl - Converts ltr_finder to gff         |
#                                                           |
#-----------------------------------------------------------+
#                                                           |
#  AUTHOR: James C. Estill                                  |
# CONTACT: JamesEstill_@_gmail.com                          |
# STARTED: 09/14/2007                                       |
# UPDATED: 09/14/2007                                       |
#                                                           |
# DESCRIPTION:                                              |
#  Converts the LTR_FINDER results to gff format.           |
#                                                           |
# VERSION: $Rev$                                            |
#                                                           |
#-----------------------------------------------------------+

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
my $infile;                    # Infile. textfile result from LTR_FINDER
my $outdir;

# Booleans
my $quiet = 0;
my $verbose = 0;
my $show_help = 0;
my $show_usage = 0;
my $show_man = 0;
my $show_version = 0;

# LTR_FINDER BOOLEANS
my $in_emp_details = 0;  # Exact match pairs details
my $in_ltr_align = 0;    # Details of the LTR alignment
my $in_pbs_align = 0;
my $in_ppt;

#
my $lf_prog_name;             # LTR Finder program name
my $lf_seq_id;                # Seq id
my $lf_seq_len;               # Sequence length
my $lf_version;               # Version of LTR finder being parsed
my $lf_trna_db;               # tRNA database used

my $lf_strand;                # Strand of the LTR
my $lf_span_start;            # Location start
my $lf_span_end;              # Location end
my $lf_length;                # Length

my $lf_score;                 # Score
my $lf_ltr_similarity;        # Similarity of the LTRs

# Status strings
my $has_5ltr_tg;                 # TG in 5' END of 5' LTR
my $has_5ltr_ca;                 # CA in 3' END of 5' LTR
my $has_3ltr_tg;                 # TG in 5' END of 3' LTR
my $has_3ltr_ca;                 # CA in 3' END of 3' LTR
my $has_tsr;                     # Has Target Site Replication
my $has_pbs;                     # Has Primer Binding Site
my $has_ppt;                     # Has Poly Purine Tract
my $has_rt;                      # Has Reverse Transcriptase
my $has_in_core;                 # Has Integrase Core
my $has_in_cterm;                # Has Integrase C-term
my $has_rh;                      # Has RNAseH

my $lf_ltr_id;                # Id number assigned to the LTR retrotransposon

# LTR COORDINATES
my $lf_5ltr_start;
my $lf_5ltr_end;
my $lf_5ltr_len;
my $lf_3ltr_start;
my $lf_3ltr_end;
my $lf_3ltr_len;

# TSR COORDINATES
my $lf_5tsr_start;             # Start of the 5' TSR
my $lf_5tsr_end;               # End of the 5' TSR
my $lf_3tsr_start;             # Start of the 3' TSR
my $lf_3tsr_end;               # End of the 3' TSR 
my $lf_tsr_string;             # String of bases in the TSR

# BOUNDARY SHARPNESS
my $lf_sharp_5;                # Sharpness of 5' Boundary
my $lf_sharp_3;                # Sharpness of 3' Boundary

#-----------------------------+
# COMMAND LINE OPTIONS        |
#-----------------------------+
my $ok = GetOptions(# REQUIRED OPTIONS
		    "i|infile=s"  => \$infile,
                    "o|outdir=s" => \$outdir,
		    # ADDITIONAL OPTIONS
		    "q|quiet"     => \$quiet,
		    "verbose"     => \$verbose,
		    # ADDITIONAL INFORMATION
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

#-----------------------------+
# OPEN INPUT FILE             |
#-----------------------------+
open (INFILE, "<$infile") ||
    die "ERROR: Can not open LTR_FINDER result file\n $infile\n";

while (<INFILE>) {
    chomp;
#    print $_."\n";



    # CHECK BOOLEANS


    # 
    if (m/Nothing Header(.*)/) {
	
    }
    
    # IN NEW REC, GET ID
    elsif (m/^\[(.*)\]/) {
	$lf_ltr_id = $1;
    }

    # SEQ ID AND LENGTH
    elsif (m/>Sequence: (.*) Len:(.*)/){
	$lf_seq_id = $1;
	$lf_seq_len = $2;
    }

    # SPAN LOCATION, LENGTH, AND STRAND
    elsif (m/^Location : (\d*) - (\d*) Len: (\d*) Strand:(.)/){
	$lf_span_start = $1;
	$lf_span_end = $2;
	$lf_length = $3;
	$lf_strand = $4;
    }

    # SCORE SIMILARITY
    elsif (m/^Score    : (.*) \[LTR region similarity:(.*)\]/){
	$lf_score = $1;
	$lf_ltr_similarity = $2;
    }

    # STATUS SET
    elsif (m/^Status   : (\d)(\d)(\d)(\d)(\d)(\d)(\d)(\d)(\d)(\d)(\d)/){
	# Since this is a binary string, it can be split as digits
	# and used to load the $has_* booleans
	$has_5ltr_tg = $1;
	$has_5ltr_ca = $2;
	$has_3ltr_tg = $3;
	$has_3ltr_ca = $4;
	$has_tsr = $5;
	$has_pbs = $6;
	$has_ppt = $7;
	$has_rt = $8;
	$has_in_core = $9;
	$has_in_cterm = $10;
	$has_rh = $11;
    }

    # 5' LTR
    elsif (m/^5\'-LTR   : (\d*) - (\d*) Len: (\d*)/){
	$lf_5ltr_start = $1;
	$lf_5ltr_end = $2;
	$lf_5ltr_len = $3;
    }

    # 3' LTR
    elsif (m/^3\'-LTR   : (\d*) - (\d*) Len: (\d*)/){
	$lf_3ltr_start = $1;
	$lf_3ltr_end = $2;
	$lf_3ltr_len = $3;
    }

    # TARGET SITE REPLICATION
    elsif (m/TSR      : (\d*) - (\d*) , (\d*) - (\d*) \[(.*)\]/){
	$lf_5tsr_start = $1;
	$lf_5tsr_end = $2;
	$lf_3tsr_start = $3;
	$lf_3tsr_end = $4;
	$lf_tsr_string = $5;
    }

    # SHARPNESS METRIC
    elsif (m/^Sharpness: (.*),(.*)/){
	$lf_sharp_5 = $1;
	$lf_sharp_3 = $2;
    }

    # PBS
    elsif (m//) {
	
    }

    #-----------------------------+
    # FILE HEADER INFORMATION     |
    #-----------------------------+

    # PROGRAM NAME
    elsif (m/^Program    : (.*)/) {
	$lf_prog_name = $1;
    }

    # PROGRAM VERSION
    elsif (m/^Version    : (.*)/) {
	$lf_version = $1;
    }

}

close INFILE;

print STDOUT "\n\n+-----------------------------+\n";
print STDOUT " RESULTS SUMMARY\n";
print STDOUT " $lf_prog_name\n";
print STDOUT " $lf_version\n";
print STDOUT "+-----------------------------+\n";

# ONLY SHOWS LAST PARSED\
print STDOUT "\n";
print STDOUT "SEQ_ID:\t\t$lf_seq_id\n";
print STDOUT "SEQ_LEN:\t$lf_seq_len\n";
print STDOUT "LTR_ID:\t\t$lf_ltr_id\n";
print STDOUT "LTRSTART:\t$lf_span_start\n";
print STDOUT "LTREND:\t\t$lf_span_end\n";
print STDOUT "LTRLEN:\t\t$lf_length\n";
print STDOUT "STRAND:\t\t$lf_strand\n";
print STDOUT "SCORE:\t\t$lf_score\n";
print STDOUT "SIMIL:\t\t$lf_ltr_similarity\n";

print STDOUT "\nLTR LOCATION:\n";
print STDOUT "5LTR START:\t$lf_5ltr_start\n";
print STDOUT "5LTR END:\t$lf_5ltr_end\n";
print STDOUT "5LTR LEN:\t$lf_5ltr_len\n";
print STDOUT "3LTR START:\t$lf_3ltr_start\n";
print STDOUT "3LTR END:\t$lf_3ltr_end\n";
print STDOUT "3LTR LEN:\t$lf_3ltr_len\n";

print STDOUT "\nTSR COORDINATES:\n";
print STDOUT "5-TSR START:\t$lf_5tsr_start\n";
print STDOUT "5-TSR END:\t $lf_5tsr_end\n";
print STDOUT "3-TSR START:\t$lf_3tsr_start\n";
print STDOUT "3-TSR END:\t$lf_3tsr_end\n";
print STDOUT "TSR-STRING:\t$lf_tsr_string\n";

print STDOUT "\nMETRICS:\n";
print STDOUT "5 SHARP: $lf_sharp_5\n";
print STDOUT "3 SHARP: $lf_sharp_3\n";

print STDOUT "\nCOMPONENT STATUS:\n";
print STDOUT "5LTR_TG\n" if $has_5ltr_tg;
print STDOUT "5LTR_CA\n" if $has_5ltr_ca;
print STDOUT "3LTR_TG\n" if $has_3ltr_tg;
print STDOUT "3LTR_CA\n" if $has_3ltr_ca;
print STDOUT "TSR\n" if $has_tsr;
print STDOUT "PBS\n" if $has_pbs;
print STDOUT "PPT\n" if $has_ppt;
print STDOUT "RT\n" if $has_rt;
print STDOUT "IN(core)\n" if $has_in_core;
print STDOUT "IN(c-term)\n" if $has_in_cterm;
print STDOUT "RH\n" if $has_rh;

exit;

#-----------------------------------------------------------+ 
# SUBFUNCTIONS                                              |
#-----------------------------------------------------------+

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

Name.pl - Short program description. 

=head1 VERSION

This documentation refers to program version $Rev$

=head1 SYNOPSIS

  USAGE:
    cnv_ltrfinder2gff.pl -i InFile -o OutDir
    
    --infile        # Path to the input file
    --outdir        # Base output dir

=head1 DESCRIPTION

Convert the ltrfinder 

=head1 COMMAND LINE ARGUMENTS

=head 2 Required Arguments

=over 2

=item -i,--infile

Path of the input file.

=item -o,--outfile

Path of the output file.

=back

=head1 Additional Options

=over

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

=head1 HISTORY

STARTED: 09/14/2007

UPDATED: 09/14/2007

VERSION: $Rev$

=cut

#-----------------------------------------------------------+
# HISTORY                                                   |
#-----------------------------------------------------------+
#
# 09/14/2007
# - Program started
