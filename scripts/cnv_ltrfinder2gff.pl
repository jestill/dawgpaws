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
# UPDATED: 09/17/2007                                       |
#                                                           |
# DESCRIPTION:                                              |
#  Converts the LTR_FINDER results to gff format.           |
#                                                           |
# VERSION: $Rev$                                      |
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
my $outfile;                   # Outfile.

# Booleans
my $quiet = 0;
my $verbose = 0;
my $show_help = 0;
my $show_usage = 0;
my $show_man = 0;
my $show_version = 0;
my $append = 0;                # Append gff output to existing file

my $suffix = "default";        # Suffix appended to the end of the gff 
                               # Source ie. 

#-----------------------------+
# COMMAND LINE OPTIONS        |
#-----------------------------+
my $ok = GetOptions(# REQUIRED OPTIONS
		    "i|infile=s"  => \$infile,
                    "o|outfile=s" => \$outfile,
		    # ADDITIONAL OPTIONS
		    "s|suffix=s"  => \$suffix,
		    "append"      => \$append,
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


# Throw error if required options not present
if ( (!$infile) || (!$outfile)  ) {
    print "\a";
    print "ERROR: Input file path required;" if !$infile;
    print "ERROR: Output file path required" if !$outfile;
    print_help("full");
}

#-----------------------------------------------------------+
# MAIN PROGRAM BODY                                         |
#-----------------------------------------------------------+

ltrfinder2gff ("HEXTEST", $infile, $outfile, $append, $suffix);

exit;

#-----------------------------------------------------------+ 
# SUBFUNCTIONS                                              |
#-----------------------------------------------------------+

sub print_help {

    # Print requested help or exit.
    # Options are to just print the full 
    my ($opt) = @_;

    my $usage = "USAGE:\n". 
	"cnv_ltrfinder2gff.pl -i InFile -o OutFile";
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


sub ltrfinder2gff {
    
    # seq_id could be extracted from the ltrfinder result
    # but passing it to the subfunction directly allows for cases
    # where the id assigned by ltr_struc differs from the way
    # the user is referring to the assembly
    # MAY WANT TO ALLOW FOR USING THE $ls_seq_id
    my ($seq_id, $lf_infile, $gffout, $do_append, $gff_suffix) = @_;
    
    # The gff src id
    my $gff_src = "ltr_finder:".$gff_suffix;
    my $gff_str_out;        # A single out string line of gff out

    # LTF FINDER COUTNERS/INDICES
    my $lf_id_num = 0;       # Incremented for every predicted model

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
    
    my $lf_ltr_id;    # Id number assigned to the LTR retrotransposon

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
    
    # PBS
    my $lf_pbs_num_match;          # Number of matched bases in the PBS
    my $lf_pbs_aln_len;            # PBS alignment length
    my $lf_pbs_start;              # Start of the PBS signal
    my $lf_pbs_end;                # End of the PBS signal
    my $lf_pbs_trna;               # PBS tRNA type and anti-codon

    # PPT 
    my $lf_ppt_num_match;          # Number of matched based in the PPT
    my $lf_ppt_aln_len;            # PPT alignment length
    my $lf_ppt_start;              # Start of the PPT
    my $lf_ppt_end;                # End of the PPT
    
    #-----------------------------+
    # DOMAIN DATA                 |
    #-----------------------------+
    
    # GENERAL DOMAIN VARS
    my $lf_domain_dom_start;
    my $lf_domain_dom_end;
    my $lf_domain_orf_start;
    my $lf_domain_orf_end;
    my $lf_domain_name;            # Type of the domain
    
    # INTEGRASE CORE
    #my $has_in_core = 0;           # Boolean for has integrase core
    my $lf_in_core_dom_start;
    my $lf_in_core_dom_end;
    my $lf_in_core_orf_start;
    my $lf_in_core_orf_end;
    
    # INTEGRASE C-TERM
    #my $has_in_cterm = 0;
    my $lf_in_cterm_dom_start;
    my $lf_in_cterm_dom_end;
    my $lf_in_cterm_orf_start;
    my $lf_in_cterm_orf_end;

    # RNASEH
    #my $has_rh = 0;
    my $lf_rh_dom_start;
    my $lf_rh_dom_end;
    my $lf_rh_orf_start;
    my $lf_rh_orf_end;
    
    # RT
    #my $has_rt = 0;
    my $lf_rt_dom_start;
    my $lf_rt_dom_end;
    my $lf_rt_orf_start;
    my $lf_rt_orf_end;
    

    #-----------------------------+
    # OPEN GFF OUTFILE            |
    #-----------------------------+
    if ($do_append) {
	open (GFFOUT, ">>$gffout") ||
	    die "ERROR: Can not open gff outfile:\n $gffout\n";
    }
    else {
	open (GFFOUT,">$gffout") ||
	    die "ERROR: Can not open gff outfile:\n $gffout\n";
    }

    #-----------------------------+
    # OPEN INPUT FILE             |
    #-----------------------------+
    open (INFILE, "<$lf_infile") ||
	die "ERROR: Can not open LTR_FINDER result file\n $lf_infile\n";

    while (<INFILE>) {
	chomp;
        #    print $_."\n";
	
	
	# CHECK BOOLEANS
	
	
	# 
	if (m/Nothing Header(.*)/) {
	    
	}
	
	#///////////////////////////////////////
	# PAY ATTENTION BELOW THIS MAY BE THE
	# PLACE TO PRINT OUT SAVED GFF DATA
	#\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

	# IN NEW REC, GET ID
	elsif (m/^\[(.*)\]/) {

	    #-----------------------------+
	    # PRINT STORED GFF OUTPUT     |
	    #-----------------------------+
	    unless ($1 == 1 ) {

		# Two alternatives here
		# keep appending to gff_lout and then print set
		# or print line at a time .. currently printing 
		# a line at a time. JCE 10/02/2007

		# FULL SPAN
		 $gff_str_out = "$lf_seq_id\t". # Seq ID
		    "$gff_src\t".                # Source
		    "LTR_retrotransposon\t".     # Data type
		    "$lf_span_start\t".          # Start
		    "$lf_span_end\t".            # End
		    "$lf_score\t".               # Score
		    "$lf_strand\t".              # Strand
		    ".\t".                       # Frame
		    "ltr_finder_$lf_ltr_id\n";   # Retro ID
		print STDOUT $gff_str_out;
		print GFFOUT $gff_str_out;
		

		if ($has_pbs) {
		    $gff_str_out = "$lf_seq_id\t".  # Seq ID
			"$gff_src\t".               # Source
			"primer_binding_site\t".    # Data type
			"$lf_pbs_start\t" .         # Start
			"$lf_pbs_end\t".            # End
			"$lf_score\t".              # Score
			"$lf_strand\t".             # Strand
			".\t".                      # Frame
			"ltr_finder_$lf_ltr_id\n";  # Retro ID
		    print STDOUT $gff_str_out;
		    print GFFOUT $gff_str_out;
		}

		
		if ($has_ppt) {
		    $gff_str_out = "$lf_seq_id\t".  # Seq ID
#		    print STDOUT "$lf_seq_id\t".    # Seq ID
			"$gff_src\t".               # Source
			"RR_tract\t".               # Data type
			"$lf_ppt_start\t".          # Start
			"$lf_ppt_end\t".            # End
			"$lf_score\t".              # Score
			"$lf_strand\t".             # Strand
			".\t".                      # Frame
			"ltr_finder_$lf_ltr_id\n";  # Retro ID
		    print STDOUT $gff_str_out;
		    print GFFOUT $gff_str_out;
		}

		
		 $gff_str_out = "$lf_seq_id\t".  # Seq ID
#		print STDOUT "$lf_seq_id\t".    # Seq ID
		    "$gff_src\t".               # Source
		    "five_prime_LTR\t".         # Data type
		    "$lf_5ltr_start\t".         # Start
		    "$lf_5ltr_end\t".           # End
		    "$lf_score\t".              # Score
		    "$lf_strand\t".             # Strand
		    ".\t".                      # Frame
		    "ltr_finder_$lf_ltr_id\n";  # Retro ID
		 print STDOUT $gff_str_out;
		 print GFFOUT $gff_str_out;

		 $gff_str_out = "$lf_seq_id\t".  # Seq ID
#		print STDOUT "$lf_seq_id\t".    # Seq ID
		    "$gff_src\t".               # Source
		    "three_prime_LTR\t".        # Data type
		    "$lf_3ltr_start\t".         # Start
		    "$lf_3ltr_end\t".           # End
		    "$lf_score\t".              # Score
		    "$lf_strand\t".             # Strand
		    ".\t".                      # Frame
		    "ltr_finder_$lf_ltr_id\n";  # Retro ID
		 print STDOUT $gff_str_out;
		 print GFFOUT $gff_str_out;

		if ($has_tsr) {
		    
		    $gff_str_out = "$lf_seq_id\t".  # Seq ID
#		    print STDOUT "$lf_seq_id\t".     # Seq ID
			"$gff_src\t".              # Source
			"target_site_duplication\t". # Data type
			"$lf_5tsr_start\t".          # Start
			"$lf_5tsr_end\t".            # End
			"$lf_score\t".               # Score
			"$lf_strand\t".              # Strand
			".\t".                       # Frame
			"ltr_finder_$lf_ltr_id\n";   # Retro ID
		    print STDOUT $gff_str_out;
		    print GFFOUT $gff_str_out;	    

		    $gff_str_out = "$lf_seq_id\t".     # Seq ID
			"$gff_src\t".              # Source
			"target_site_duplication\t". # Data type
			"$lf_3tsr_start\t".          # Start
			"$lf_3tsr_end\t".            # End
			"$lf_score\t".               # Score
			"$lf_strand\t".              # Strand
			".\t".                       # Frame
			"ltr_finder_$lf_ltr_id\n";   # Retro ID
		    print STDOUT $gff_str_out;
		    print GFFOUT $gff_str_out;
		    
		}


		# Integrase Core
		if ($has_in_core) {
		    #/////////
		    # NOT SONG
		    #\\\\\\\\\
		    $gff_str_out = "$lf_seq_id\t".     # Seq ID
			"$gff_src\t".              # Source
			"integrase_core_domain\t".  # Data type
			"$lf_in_core_dom_start\t".   # Start
			"$lf_in_core_dom_end\t".     # End
			"$lf_score\t".               # Score
			"$lf_strand\t".              # Strand
			".\t".                       # Frame
			"ltr_finder_$lf_ltr_id\n";   # Retro ID
		    print STDOUT $gff_str_out;
		    print GFFOUT $gff_str_out;

		    #/////////
		    # NOT SONG
		    #\\\\\\\\\
		    $gff_str_out = "$lf_seq_id\t".     # Seq ID
			"$gff_src\t".              # Source
			"integrase_core_orf\t".      # Data type
			"$lf_in_core_orf_start\t".   # Start
			"$lf_in_core_orf_end\t".     # End
			"$lf_score\t".               # Score
			"$lf_strand\t".              # Strand
			".\t".                       # Frame
			"ltr_finder_$lf_ltr_id\n";   # Retro ID
		    print STDOUT $gff_str_out;
		    print GFFOUT $gff_str_out;
		    
		} # End of has in_core


		if ($has_in_cterm) {

		    #/////////
		    # NOT SONG
		    #\\\\\\\\\
		    $gff_str_out = "$lf_seq_id\t".     # Seq ID
			"$gff_src\t".                # Source
			"integrase_cterm_domain\t".  # Data type
			"$lf_in_cterm_dom_start\t".  # Start
			"$lf_in_cterm_dom_end\t".    # End
			"$lf_score\t".               # Score
			"$lf_strand\t".              # Strand
			".\t".                       # Frame
			"ltr_finder_$lf_ltr_id\n";   # Retro ID
		    print STDOUT $gff_str_out;
		    print GFFOUT $gff_str_out;

		    $gff_str_out = "$lf_seq_id\t".     # Seq ID
			"$gff_src\t".                # Source
			"integrase_cterm_orf\t".     # Data type
			"$lf_in_cterm_orf_start\t".  # Start
			"$lf_in_cterm_orf_end\t".    # End
			"$lf_score\t".               # Score
			"$lf_strand\t".              # Strand
			".\t".                       # Frame
			"ltr_finder_$lf_ltr_id\n";   # Retro ID
		    print STDOUT $gff_str_out;
		    print GFFOUT $gff_str_out;
		    
		} # End of has_in_cterm

		if ($has_rh) {

		    #/////////
		    # NOT SONG
		    #\\\\\\\\\
		    $gff_str_out = "$lf_seq_id\t".     # Seq ID
			"$gff_src\t".                # Source
			"rnaseh_domain\t".           # Data type
			"$lf_rh_dom_start\t".        # Start
			"$lf_rh_dom_end\t".          # End
			"$lf_score\t".               # Score
			"$lf_strand\t".              # Strand
			".\t".                       # Frame
			"ltr_finder_$lf_ltr_id\n";   # Retro ID
		    print STDOUT $gff_str_out;
		    print GFFOUT $gff_str_out;

		    $gff_str_out = "$lf_seq_id\t".     # Seq ID
			"$gff_src\t".                # Source
			"rnaseh_orf\t".              # Data type
			"$lf_rh_orf_start\t".        # Start
			"$lf_rh_orf_end\t".          # End
			"$lf_score\t".               # Score
			"$lf_strand\t".              # Strand
			".\t".                       # Frame
			"ltr_finder_$lf_ltr_id\n";   # Retro ID
		    print STDOUT $gff_str_out;
		    print GFFOUT $gff_str_out;
		    
		}

		if ($has_rt) {

		    #/////////
		    # NOT SONG
		    #\\\\\\\\\
		    $gff_str_out = "$lf_seq_id\t".     # Seq ID
			"$gff_src\t".                # Source
			"rt_domain\t".               # Data type
			"$lf_rt_dom_start\t".        # Start
			"$lf_rt_dom_end\t".          # End
			"$lf_score\t".               # Score
			"$lf_strand\t".              # Strand
			".\t".                       # Frame
			"ltr_finder_$lf_ltr_id\n";   # Retro ID
		    print STDOUT $gff_str_out;
		    print GFFOUT $gff_str_out;

		    $gff_str_out = "$lf_seq_id\t".     # Seq ID
			"$gff_src\t".                # Source
			"rt_orf\t".                  # Data type
			"$lf_rt_orf_start\t".        # Start
			"$lf_rt_orf_end\t".          # End
			"$lf_score\t".               # Score
			"$lf_strand\t".              # Strand
			".\t".                       # Frame
			"ltr_finder_$lf_ltr_id\n";   # Retro ID
		    print STDOUT $gff_str_out;
		    print GFFOUT $gff_str_out;

		} # End of has_reverse_transcriptase

	    } # End of unless this is the first record


	    #-----------------------------+
	    # RESET VARS TO NULL          |
	    #-----------------------------+
	    # May not need to reset vars since existence of the vars
	    # is designated by the do_* variables.
	    
	    #-----------------------------+
	    # LOAD ID VAR                 |
	    #-----------------------------+
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
	elsif (m/PBS   : \[(\d*)\/(\d*)\] (\d*) - (\d*) \((.*)\)/) {
	    $lf_pbs_num_match = $1;
	    $lf_pbs_aln_len = $2;
	    $lf_pbs_start = $3;
	    $lf_pbs_end = $4;
	    $lf_pbs_trna = $5;
	}
	
	# PPT
	elsif (m/PPT   : \[(\d*)\/(\d*)\] (\d*) - (\d*)/) {
	    $lf_ppt_num_match = $1;
	    $lf_ppt_aln_len = $2;
	    $lf_ppt_start = $3;
	    $lf_ppt_end = $4;
	}
	
	# PROTEIN DOMAINS
	# This will need to be modified and checked after another run
	# using ps_scan to get the additional domains
	#
	#Domain: 56796 - 57326 [possible ORF:56259-59144, (IN (core))]
	elsif (m/Domain: (\d*) - (\d*) \[possible ORF:(\d*)-(\d*), \((.*)\)\]/) {
	    
	    $lf_domain_dom_start = $1;
	    $lf_domain_dom_end = $2;
	    $lf_domain_orf_start = $3;
	    $lf_domain_orf_end = $4;
	    $lf_domain_name = $5;
	    
	    # Temp while I work with this data
	    #print "DOMAIN:".$lf_domain_name."\n";
	    
	    if ($lf_domain_name =~ 'IN \(core\)') {
		
		$lf_in_core_dom_start = $lf_domain_dom_start;
		$lf_in_core_dom_end = $lf_domain_dom_end;
		$lf_in_core_orf_start = $lf_domain_orf_start;
		$lf_in_core_orf_end = $lf_domain_orf_end;
		
		# Temp while debug
		# This is to check the regexp vars can be fetched here
		#print "\tDom Start: $lf_in_core_dom_start\n";
		#print "\tDom End:   $lf_in_core_dom_end\n";
		#print "\tOrf Start: $lf_in_core_orf_start\n";
		#print "\tOrf End:   $lf_in_core_orf_end\n";
		
	    }
	    elsif ($lf_domain_name =~ 'IN \(c-term\)') {
		$lf_in_cterm_dom_start = $lf_domain_dom_start;
		$lf_in_cterm_dom_end = $lf_domain_dom_end;
		$lf_in_cterm_orf_start = $lf_domain_orf_start;
		$lf_in_cterm_orf_end = $lf_domain_orf_end;
	    }
	    elsif ($lf_domain_name =~ 'RH') {
		$lf_rh_dom_start = $lf_domain_dom_start;
		$lf_rh_dom_end = $lf_domain_dom_end;
		$lf_rh_orf_start = $lf_domain_orf_start;
		$lf_rh_orf_end = $lf_domain_orf_end;
	    }
	    elsif ($lf_domain_name =~ 'RT') {
		
		$lf_rt_dom_start = $lf_domain_dom_start;
		$lf_rt_dom_end = $lf_domain_dom_end;
		$lf_rt_orf_start = $lf_domain_orf_start;
		$lf_rt_orf_end = $lf_domain_orf_end;  
		
	    }
	    else {
		print "\a";
		print STDERR "Unknown domain type: $lf_domain_name\n";
	    }
	    
	    
	    
	} # End of elsif Domain
	
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
    
    #-----------------------------+
    # PRINT LAST GFFOUT           |
    #-----------------------------+
    
    # FULL SPAN
    $gff_str_out = "$lf_seq_id\t". # Seq ID
	"$gff_src\t".                # Source
	"LTR_retrotransposon\t".     # Data type
	"$lf_span_start\t".          # Start
	"$lf_span_end\t".            # End
	"$lf_score\t".               # Score
	"$lf_strand\t".              # Strand
	".\t".                       # Frame
	"ltr_finder_$lf_ltr_id\n";   # Retro ID
    print STDOUT $gff_str_out;
    print GFFOUT $gff_str_out;
    
    
    if ($has_pbs) {
	$gff_str_out = "$lf_seq_id\t".  # Seq ID
	    "$gff_src\t".               # Source
	    "primer_binding_site\t".    # Data type
	    "$lf_pbs_start\t" .         # Start
	    "$lf_pbs_end\t".            # End
	    "$lf_score\t".              # Score
	    "$lf_strand\t".             # Strand
	    ".\t".                      # Frame
	    "ltr_finder_$lf_ltr_id\n";  # Retro ID
	print STDOUT $gff_str_out;
	print GFFOUT $gff_str_out;
    }
    
    
    if ($has_ppt) {
	$gff_str_out = "$lf_seq_id\t".  # Seq ID
#		    print STDOUT "$lf_seq_id\t".    # Seq ID
	    "$gff_src\t".               # Source
	    "RR_tract\t".               # Data type
	    "$lf_ppt_start\t".          # Start
	    "$lf_ppt_end\t".            # End
	    "$lf_score\t".              # Score
	    "$lf_strand\t".             # Strand
	    ".\t".                      # Frame
	    "ltr_finder_$lf_ltr_id\n";  # Retro ID
	print STDOUT $gff_str_out;
	print GFFOUT $gff_str_out;
    }
    
    
    $gff_str_out = "$lf_seq_id\t".  # Seq ID
#		print STDOUT "$lf_seq_id\t".    # Seq ID
	"$gff_src\t".               # Source
	"five_prime_LTR\t".         # Data type
	"$lf_5ltr_start\t".         # Start
	"$lf_5ltr_end\t".           # End
	"$lf_score\t".              # Score
	"$lf_strand\t".             # Strand
	".\t".                      # Frame
	"ltr_finder_$lf_ltr_id\n";  # Retro ID
    print STDOUT $gff_str_out;
    print GFFOUT $gff_str_out;
    
    $gff_str_out = "$lf_seq_id\t".  # Seq ID
#		print STDOUT "$lf_seq_id\t".    # Seq ID
	"$gff_src\t".               # Source
	"three_prime_LTR\t".        # Data type
	"$lf_3ltr_start\t".         # Start
	"$lf_3ltr_end\t".           # End
	"$lf_score\t".              # Score
	"$lf_strand\t".             # Strand
	".\t".                      # Frame
	"ltr_finder_$lf_ltr_id\n";  # Retro ID
    print STDOUT $gff_str_out;
    print GFFOUT $gff_str_out;
    
    if ($has_tsr) {
	
	$gff_str_out = "$lf_seq_id\t".  # Seq ID
#		    print STDOUT "$lf_seq_id\t".     # Seq ID
	    "$gff_src\t".              # Source
	    "target_site_duplication\t". # Data type
	    "$lf_5tsr_start\t".          # Start
	    "$lf_5tsr_end\t".            # End
	    "$lf_score\t".               # Score
	    "$lf_strand\t".              # Strand
	    ".\t".                       # Frame
	    "ltr_finder_$lf_ltr_id\n";   # Retro ID
	print STDOUT $gff_str_out;
	print GFFOUT $gff_str_out;	    
	
	$gff_str_out = "$lf_seq_id\t".     # Seq ID
	    "$gff_src\t".              # Source
	    "target_site_duplication\t". # Data type
	    "$lf_3tsr_start\t".          # Start
	    "$lf_3tsr_end\t".            # End
	    "$lf_score\t".               # Score
	    "$lf_strand\t".              # Strand
	    ".\t".                       # Frame
	    "ltr_finder_$lf_ltr_id\n";   # Retro ID
	print STDOUT $gff_str_out;
	print GFFOUT $gff_str_out;
	
    }
    
    
    # Integrase Core
    if ($has_in_core) {
	#/////////
	# NOT SONG
	#\\\\\\\\\
	$gff_str_out = "$lf_seq_id\t".     # Seq ID
	    "$gff_src\t".              # Source
	    "integrase_core_domain\t".  # Data type
	    "$lf_in_core_dom_start\t".   # Start
	    "$lf_in_core_dom_end\t".     # End
	    "$lf_score\t".               # Score
	    "$lf_strand\t".              # Strand
	    ".\t".                       # Frame
	    "ltr_finder_$lf_ltr_id\n";   # Retro ID
	print STDOUT $gff_str_out;
	print GFFOUT $gff_str_out;
	
	#/////////
	# NOT SONG
	#\\\\\\\\\
	$gff_str_out = "$lf_seq_id\t".     # Seq ID
	    "$gff_src\t".              # Source
	    "integrase_core_orf\t".      # Data type
	    "$lf_in_core_orf_start\t".   # Start
	    "$lf_in_core_orf_end\t".     # End
	    "$lf_score\t".               # Score
	    "$lf_strand\t".              # Strand
	    ".\t".                       # Frame
	    "ltr_finder_$lf_ltr_id\n";   # Retro ID
	print STDOUT $gff_str_out;
	print GFFOUT $gff_str_out;
	
    } # End of has in_core
    
    
    if ($has_in_cterm) {
	
	#/////////
	# NOT SONG
	#\\\\\\\\\
	$gff_str_out = "$lf_seq_id\t".     # Seq ID
	    "$gff_src\t".                # Source
	    "integrase_cterm_domain\t".  # Data type
	    "$lf_in_cterm_dom_start\t".  # Start
	    "$lf_in_cterm_dom_end\t".    # End
	    "$lf_score\t".               # Score
	    "$lf_strand\t".              # Strand
	    ".\t".                       # Frame
	    "ltr_finder_$lf_ltr_id\n";   # Retro ID
	print STDOUT $gff_str_out;
	print GFFOUT $gff_str_out;
	
	$gff_str_out = "$lf_seq_id\t".     # Seq ID
	    "$gff_src\t".                # Source
	    "integrase_cterm_orf\t".     # Data type
	    "$lf_in_cterm_orf_start\t".  # Start
	    "$lf_in_cterm_orf_end\t".    # End
	    "$lf_score\t".               # Score
	    "$lf_strand\t".              # Strand
	    ".\t".                       # Frame
	    "ltr_finder_$lf_ltr_id\n";   # Retro ID
	print STDOUT $gff_str_out;
	print GFFOUT $gff_str_out;
	
    } # End of has_in_cterm
    
    if ($has_rh) {
	
	#/////////
	# NOT SONG
	#\\\\\\\\\
	$gff_str_out = "$lf_seq_id\t".     # Seq ID
	    "$gff_src\t".                # Source
	    "rnaseh_domain\t".           # Data type
	    "$lf_rh_dom_start\t".        # Start
	    "$lf_rh_dom_end\t".          # End
	    "$lf_score\t".               # Score
	    "$lf_strand\t".              # Strand
	    ".\t".                       # Frame
	    "ltr_finder_$lf_ltr_id\n";   # Retro ID
	print STDOUT $gff_str_out;
	print GFFOUT $gff_str_out;
	
	$gff_str_out = "$lf_seq_id\t".     # Seq ID
	    "$gff_src\t".                # Source
	    "rnaseh_orf\t".              # Data type
	    "$lf_rh_orf_start\t".        # Start
	    "$lf_rh_orf_end\t".          # End
	    "$lf_score\t".               # Score
	    "$lf_strand\t".              # Strand
	    ".\t".                       # Frame
	    "ltr_finder_$lf_ltr_id\n";   # Retro ID
	print STDOUT $gff_str_out;
	print GFFOUT $gff_str_out;
	
    }
    
    if ($has_rt) {
	
	#/////////
	# NOT SONG
	#\\\\\\\\\
	$gff_str_out = "$lf_seq_id\t".     # Seq ID
	    "$gff_src\t".                # Source
	    "rt_domain\t".               # Data type
	    "$lf_rt_dom_start\t".        # Start
	    "$lf_rt_dom_end\t".          # End
	    "$lf_score\t".               # Score
	    "$lf_strand\t".              # Strand
	    ".\t".                       # Frame
	    "ltr_finder_$lf_ltr_id\n";   # Retro ID
	print STDOUT $gff_str_out;
	print GFFOUT $gff_str_out;
	
	$gff_str_out = "$lf_seq_id\t".     # Seq ID
	    "$gff_src\t".                # Source
	    "rt_orf\t".                  # Data type
	    "$lf_rt_orf_start\t".        # Start
	    "$lf_rt_orf_end\t".          # End
	    "$lf_score\t".               # Score
	    "$lf_strand\t".              # Strand
	    ".\t".                       # Frame
	    "ltr_finder_$lf_ltr_id\n";   # Retro ID
	print STDOUT $gff_str_out;
	print GFFOUT $gff_str_out;
	
    } # End of has_reverse_transcriptase
    
    
    #-----------------------------+
    # PRINT ADDITIONAL DATA       |
    # IF REQUESTED                |
    #-----------------------------+
    if ($verbose) {
	print STDERR "\n\n+-----------------------------+\n";
	print STDERR " RESULTS SUMMARY\n";
	print STDERR " $lf_prog_name\n";
	print STDERR " $lf_version\n";
	print STDERR "+-----------------------------+\n";
	
	# ONLY SHOWS LAST PARSED\
	print STDERR "\n";
	print STDERR "SEQ_ID:\t\t$lf_seq_id\n";
	print STDERR "SEQ_LEN:\t$lf_seq_len\n";
	print STDERR "LTR_ID:\t\t$lf_ltr_id\n";
	print STDERR "LTRSTART:\t$lf_span_start\n";
	print STDERR "LTREND:\t\t$lf_span_end\n";
	print STDERR "LTRLEN:\t\t$lf_length\n";
	print STDERR "STRAND:\t\t$lf_strand\n";
	print STDERR "SCORE:\t\t$lf_score\n";
	print STDERR "SIMIL:\t\t$lf_ltr_similarity\n";
	
	print STDERR "\nLTR LOCATION:\n";
	print STDERR "5LTR START:\t$lf_5ltr_start\n";
	print STDERR "5LTR END:\t$lf_5ltr_end\n";
	print STDERR "5LTR LEN:\t$lf_5ltr_len\n";
	print STDERR "3LTR START:\t$lf_3ltr_start\n";
	print STDERR "3LTR END:\t$lf_3ltr_end\n";
	print STDERR "3LTR LEN:\t$lf_3ltr_len\n";
	
	print STDERR "\nTSR COORDINATES:\n";
	print STDERR "5-TSR START:\t$lf_5tsr_start\n";
	print STDERR "5-TSR END:\t$lf_5tsr_end\n";
	print STDERR "3-TSR START:\t$lf_3tsr_start\n";
	print STDERR "3-TSR END:\t$lf_3tsr_end\n";
	print STDERR "TSR-STRING:\t$lf_tsr_string\n";
    
	print STDERR "\nMETRICS:\n";
	print STDERR "5 SHARP: $lf_sharp_5\n";
	print STDERR "3 SHARP: $lf_sharp_3\n";
	
	print STDERR "\nPBS:\n";
	print STDERR "NUM MATCH:\t$lf_pbs_num_match\n" if $lf_pbs_num_match;
	print STDERR "ALN LENGT:\t$lf_pbs_aln_len\n" if $lf_pbs_aln_len;
	print STDERR "PBS START:\t$lf_pbs_start\n" if $lf_pbs_start;
	print STDERR "PBS END:\t$lf_pbs_end\n" if $lf_pbs_end;
	print STDERR "PBS tRNA:\t$lf_pbs_trna\n" if $lf_pbs_trna;
	
	print STDERR "\nPPT:\n";
	print STDERR "NUM MATCH:\t$lf_ppt_num_match\n";
	print STDERR "ALN LENT:\t$lf_ppt_aln_len\n";
	print STDERR "PPT START:\t$lf_ppt_start\n";
	print STDERR "PPT END:\t$lf_ppt_end\n";
	
	print STDERR "\nCOMPONENT STATUS:\n";
	print STDERR "5LTR_TG\n" if $has_5ltr_tg;
	print STDERR "5LTR_CA\n" if $has_5ltr_ca;
	print STDERR "3LTR_TG\n" if $has_3ltr_tg;
	print STDERR "3LTR_CA\n" if $has_3ltr_ca;
	print STDERR "TSR\n" if $has_tsr;
	print STDERR "PBS\n" if $has_pbs;
	print STDERR "PPT\n" if $has_ppt;
	print STDERR "RT\n" if $has_rt;
	print STDERR "IN(core)\n" if $has_in_core;
	print STDERR "IN(c-term)\n" if $has_in_cterm;
	print STDERR "RH\n" if $has_rh;
    }
    
    close GFFOUT;

} # End of ltrfinder2gff subfunction

=head1 NAME

cnv_ltrfinder2gff.pl - Converts LTR_Finder output to gff format

=head1 VERSION

This documentation refers to program version $Rev$

=head1 SYNOPSIS

  USAGE:
    cnv_ltrfinder2gff.pl -i InFile -o OutDir
    
    --infile        # Path to the input file
                    # Result from a single record fasta file
    --outdir        # Base output dir

=head1 DESCRIPTION

Convert the ltrfinder output to gff format. This assumes that the output
from ltr_finder correspons to a single BAC. A suffix can be passed with the
--suffix option to provide a suffix for the source column. For example
run default ltr_finder results with --suffix def to create 
ltr_finder:def. T

In addition to printing to the GFF file handle, this will also print
all results to STDOUT.

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

UPDATED: 10/02/2007

VERSION: $Rev$

=cut

#-----------------------------------------------------------+
# HISTORY                                                   |
#-----------------------------------------------------------+
#
# 09/14/2007
# - Program started
# 10/01/2007
# - Moved main body of program to subfunction ltrfinder2gff
#   This subfunction operates under the assumption that all
#   results in the ltrfinder output file are the results for
#   a single assembly. This will facilitate copying the results
#   to an individual directory for the contig being annotated.
# - Making gff output, making this song complient
# 10/02/2007
# - Finishing gff output, now saving to string to write to
#   both gffout and stdout.
