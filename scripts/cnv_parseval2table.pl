#!/usr/bin/perl -w
#-----------------------------------------------------------+
#                                                           |
# cnv_parseval2grid.pl
#                                                           |
#-----------------------------------------------------------+
#                                                           |
#  AUTHOR: James C. Estill                                  |
# CONTACT: JamesEstill_@_gmail.com                          |
# STARTED: 05/02/2012                                       |
# UPDATED: 05/03/2012                                       |
#                                                           |
# DESCRIPTION:                                              |
#  Convert a directory of parseeval outputs to a summary    |
#  file(s) that summarize the results among the different   |
#  gene prediction results.                                 |
#                                                           |
# USAGE:                                                    |
#  cnv_parseval2grid.pl -i dir_of_results/ -o summary_file  |
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

#-----------------------------+
# VARIABLE SCOPE              |
#-----------------------------+
my $in_path;
my $outfile;
my $format = "txt";            # Expect format to parse is a text file
                               # This allows for parseval to eventually
                               # produce xml format.

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
		    "i|infile|indir=s"   => \$in_path,
                    "o|outfile=s"        => \$outfile,
		    # ADDITIONAL OPTIONS
		    "q|quiet"            => \$quiet,
		    "verbose"            => \$verbose,
		    # ADDITIONAL INFORMATION
		    "usage"              => \$show_usage,
		    "test"               => \$do_test,
		    "version"            => \$show_version,
		    "man"                => \$show_man,
		    "h|help"             => \$show_help,);

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
# CHECK REQUIRED ARGS         |
#-----------------------------+
if ( (!$in_path)  ) {
    print "\a";
    print STDERR "\n";
    print STDERR "ERROR: An input directory or file was not specified at the".
	" command line\n" if (!$in_path);
    print_help ("usage", $0 );
}

#-----------------------------------------------------------+
# LOAD THE ARRAY OF FILE PATHS                              |
#-----------------------------------------------------------+
my @input_files;
if ($in_path) {
    if (-f $in_path) {
	print STDERR "Input path is a file\n"
	    if $verbose;
	push (@input_files, $in_path);
    }
    elsif (-d $in_path) {
	
	# NOTE: If other input formats are added, change the following to always
	# default to fasta format. Current here to allow for
	# input from other types of data.
	print STDERR "Input path is a directory\n" 
	    if $verbose;
	
	# GET THE DIRECTORY VAR
	my $in_dir = $in_path;
	# Add slash to indir if needed
	unless ($in_dir =~ /\/$/ ) {
	    $in_dir = $in_dir."/";
	}
	
	# LOAD FILES IN THE INTPUT DIRECTORY
	# First load to tmp array so that indir can be prefixed to inpath
	my @tmp_file_paths;
	# txt format is a tab delimited text file
	if ($format =~ "txt") {
	    opendir( DIR, $in_dir ) || 
		die "Can't open directory:\n$in_dir"; 
	    @tmp_file_paths = grep /\.text$|\.txt$/, readdir DIR ;
	    closedir( DIR );
	}
	elsif ($format =~ "xml") {
	    opendir( DIR, $in_dir ) || 
		die "Can't open directory:\n$in_dir"; 
	    @tmp_file_paths = grep /\.xml$|\.XML$/, readdir DIR ;
	    closedir( DIR );
	}
	
	# DIR directory to path of 
	foreach my $tmp_file_path (@tmp_file_paths ) {
	    push (@input_files, $in_dir.$tmp_file_path);
	}
	
	# If no files found matching expected extensions, may want to
	# just push all files in the directory
	
	
    } else {
	print STDERR "Input path is not a valid directory or file:\n";
	die;
    }
    
} else {
    print STDERR "\a";
    print STDERR "WARNING: A input directory or file has not been specified\n";
}

#-----------------------------------------------------------+
# OUTPUT FILE HANDLE                                        |
#-----------------------------------------------------------+
if ($outfile) {
    open (OUTFILE, ">$outfile") ||
	die "Can not open outfile $outfile"
}
else {
    open (OUTFILE, ">&STDOUT") ||
	die "Can not STDOUT";
}

# Print header for output file
print OUTFILE "Reference\t".
    "Prediction\t".
    "Shared Genes\t".
    "Unique to Reference\t".
    "Unique to Prediction\t".
    "Number Perfect Matches\t".
    "Percent Perfect Matches\t".
    "Number Non-Matches\t".
    "Percent Non-Matches\t".
    "Perfect Match Length\t".
    "Non-match Length\t".
    "Perfect Avgnumrefexons\t".
    "Perfect AvgNumRefExons\t".
    "Non-match Avgnumrefexons\t".
    "Non-match AvgNumRefExons\t".
    "CDS Sensitivity\t".
    "CDS Specificity\t".
    "CDS F1 Score\t".
    "CDS AED\t".
    "Exon Sensitivity\t".
    "Exon Specificity\t".
    "Exon F1 Score\t".
    "Exon AED\t".
    "Nuc CDS Match Coefficient\t".
    "Nuc UTR Match Coefficient\t".
    "Nuc Overall Match Coefficient\t".
    "Nuc CDS Specificity\t".
    "Nuc CDS Sensitivity\t".
    "Nuc CDS F1\t".
    "Nuc CDS AED\t".
    "\n";

# May want to print outfile header here

#-----------------------------------------------------------+
# MAIN PROGRAM BODY                                         |
#-----------------------------------------------------------+
foreach my $infile (@input_files) {

    my $ref_ann = "--";
    my $pred_ann = "--";

    my $num_genes_shared = "--";
    my $num_genes_u_ref = "--";
    my $num_genes_u_pred = "--";
    
    my $num_total_comparisons = "--";

    # Gene perfect matches
    my $percent_perfect_matches = "--";
    my $num_perfect_matches = "--";
    my $perfect_match_len = "--";
    my $perfect_match_ref_num_exons = "--";
    my $perfect_match_num_pred_exons = "--";
    my $perfect_match_ref_cds_len = "--";
    my $perfect_match_pred_cds_len = "--";

    # Gene non-match
    my $num_non_matches = "--";
    my $percent_non_matches = "--";
    my $non_match_len = "--";
    my $non_match_ref_num_exons = "--";
    my $non_match_num_pred_exons = "--";
    my $non_match_ref_cds_len = "--";
    my $non_match_pred_cds_len = "--";

    # File location booleans
    my $in_gene = 0;                      # Boolean for in gene information
    my $in_perfect_matches = 0;           # Boolean for in perfect match data
    my $in_perfect_matches_mis_utr = 0;   # Perfect match with mislabeled UTR
    my $in_cds_structure_match = 0;       # Perfect match with mislabeled UTR
    my $in_exon_structure_match = 0;      # Exon structure matches
    my $in_utr_structure_match = 0;       # UTR structure matches
    my $in_non_matches = 0;               # Boolean for in non-match data
    my $in_cds = 0;
    my $in_exon = 0;
    my $in_utr = 0;
    my $in_nuc = 0; # In nucleotide level comparison

    # CDS INFO
    my $cds_f1_score = "--";
    my $cds_aed = "--";
    my $cds_sensitivity = "--";
    my $cds_specificity = "--";
    
    # EXON INFO
    my $exon_f1_score = "--";
    my $exon_aed = "--";
    my $exon_sensitivity = "--";
    my $exon_specificity = "--";

    # UTR INFO
    my $utr_f1_score = "--";
    my $utr_aed = "--";
    my $utr_sensitivity = "--";
    my $utr_specificity = "--";

    # Nucleotide level info
    my $cds_match_coef = "--";
    my $utr_match_coef = "--";
    my $overall_match_coef = "--";

    my $cds_cor_coef = "--";
    my $utr_cor_coef = "--";
    my $overall_cor_coef = "--";
    
    my $cds_nuc_specificity = "--";
    my $utr_nuc_specificity = "--";
    my $overall_nuc_specificity = "--";

    my $cds_nuc_sensitivity = "--";
    my $utr_nuc_sensitivity = "--";
    my $overall_nuc_sensitivity = "--";

    my $cds_nuc_f1 = "--";
    my $utr_nuc_f1 = "--";
    my $overall_nuc_f1 = "--";

    my $cds_nuc_aed = "--";
    my $utr_nuc_aed = "--";
    my $overall_aed = "--";

    print STDERR "\n============================================\n"
	if $verbose;
    print STDERR "Processing:\n$infile\n"
	if $verbose;
    print STDERR "============================================\n"
	if $verbose;

    open(INFILE, "<$infile") ||
	die "Can not open infile $infile";

    while (<INFILE>) {
	chomp;


	# Skip all of the individual locus stuff
	if ( m/ ^\| / ) { next; }
	if (m /^     \|/) { next; }
	
	if ( m/Reference annotations:  (.*)/ ) {
	    print STDERR "\tReference: ".$1."\n"
		if $verbose;
	    $ref_ann = $1;
	}
	elsif ( m/Prediction annotations: (.*)/ ) {
	    print STDERR "\tPrediction: ".$1."\n"
		if $verbose;
	    $pred_ann=$1;
	}
	elsif ( m/  Gene loci/ ) {
	    print STDERR "\tin genes\n" 
		if $verbose;
	}
	elsif ( m /    shared\.+(.*)/ ) {
	    print STDERR "\t\tshared genes: ".$1."\n"
		if $verbose;
	    $num_genes_shared = $1;
	}
	elsif ( m/    unique to reference\.+(.*)/ ) {
	    print STDERR "\t\tunique to reference: ".$1."\n"
		if $verbose;
	    $num_genes_u_ref = $1;
	}
	elsif ( m/    unique to prediction\.+(.*)/ ) {
	    print STDERR "\t\tunique to prediction: ".$1."\n"
		if $verbose;
	    $num_genes_u_ref = $1;
	}
	elsif ( m/  Total comparisons\.+(.*)/ ) {
	    # This appears to be number of predicted genes
	    # minus the number unique to the predictin
	    # ir 15146
	    print STDERR "\tTotal Comparisions: ".$1."\n"
		if $verbose;
	     $num_total_comparisons = $1;
	}
	elsif ( m/    perfect matches\.+(.*) \((.*)\%\)/ ) {
	    print STDERR "\n\tNum perfect matches: ".$1."\n"
		if $verbose;
	    print STDERR "\tPercent perfect matches: ".$2."\n"
		if $verbose;
	    $num_perfect_matches = $1;
	    $percent_perfect_matches = $2;
	    $in_perfect_matches = 1;
	}
	elsif ( m/perfect matches with mislabeled UTRs\.+(.*) \((.*)\%\)/ ) {
	    $in_perfect_matches = 0;
	    $in_perfect_matches_mis_utr = 1;
	}
	elsif ( m/    non-matches\.+(.*) \((.*)\%\)/ ) {
	    print STDERR "\n\tNum non-matches: ".$1."\n"
		if $verbose;
	    print STDERR "\tPercent non-matches: ".$2."\n"
		if $verbose;
	    $num_non_matches = $1;
	    $percent_non_matches = $2;
	    $in_perfect_matches = 0;
	    $in_non_matches = 1;
	}
	elsif ( m/      avg. length\.+(.*)/ ) {
	    my $avg_len = $1;
	    if ($in_perfect_matches) {
		$perfect_match_len = $avg_len;
		print STDERR "\t\tPerfect match len: ".$perfect_match_len."\n"
		    if $verbose;
	    }
	    elsif ($in_non_matches) {
		$non_match_len = $avg_len;
		print STDERR "\t\tNon-match len: ".$non_match_len."\n"
		    if $verbose;
	    }
	}
	elsif ( m/      avg. # refr exons\.+(.*)/ ) {
	    my $num_exons = $1;
	    if ($in_perfect_matches) {
		$perfect_match_ref_num_exons = $num_exons;
		print STDERR "\t\tPerfect match ref exon: ".
		    $perfect_match_ref_num_exons."\n"
		    if $verbose;
	    }
	    elsif ($in_non_matches) {
		$non_match_ref_num_exons = $num_exons;
		print STDERR "\t\tNon-match ref exon: ".
		    $non_match_ref_num_exons."\n"
		    if $verbose;
	    }
	}
	elsif ( m/      avg. # pred exons\.+(.*)/ ) {
	    my $num_exon = $1;
	    if ($in_perfect_matches) {
		$perfect_match_num_pred_exons = $num_exon;
		print STDERR "\t\tPerfect match pred exon: ".
		    $perfect_match_num_pred_exons."\n"
		    if $verbose;
	    }
	    elsif ($in_non_matches) {
		$non_match_num_pred_exons = $num_exon;
		print STDERR "\t\tNon-match pred exon: ".
		    $non_match_num_pred_exons."\n"
		    if $verbose;
	    }
	}
	elsif ( m/      avg. refr CDS length\.+(.*)/ ) {
	    my $cds_len = $1;
	    if ($in_perfect_matches) {
		$perfect_match_ref_cds_len = $cds_len;
		print STDERR "\t\tPerfect match ref CDS len: ".
		    $perfect_match_ref_cds_len."\n"
		    if $verbose;
	    }
	    elsif ($in_non_matches) {
		$non_match_ref_cds_len = $cds_len;
		print STDERR "\t\tNon-match ref CDS len: ".
		    $non_match_ref_cds_len."\n"
		    if $verbose;
	    }
	}
	elsif ( m/avg. pred CDS length\.+(.*)/ ) {
	    my $cds_len = $1;
	    if ($in_perfect_matches) {
		$perfect_match_pred_cds_len = $cds_len;
		print STDERR "\t\tPerfect match pred CDS len: ".
		    $perfect_match_pred_cds_len."\n"
		    if $verbose;
	    }
	    elsif ($in_non_matches) {
		$non_match_pred_cds_len = $cds_len;
		print STDERR "\t\tNon-match pred CDS len: ".
		    $non_match_pred_cds_len."\n"
		    if $verbose;
	    }
	}
	elsif ( m/  CDS structure comparison/ ) {
	    $in_cds = 1;
	    $in_non_matches = 0;
	    print STDERR "\n\tIn CDS\n" 
		if $verbose;
	}
	elsif ( m/  Exon structure comparison/ ) {
	    print STDERR "\n\tIn Exon\n"
		if $verbose;
	    $in_exon = 1;
	    $in_cds = 0;
	}
	elsif ( m/  UTR structure comparison/ ) {
	    print STDERR "\n\tIn UTR\n"
		if $verbose;
	    $in_utr = 1;
	    $in_exon = 0;
	}
	elsif ( m/  Nucleotide-level comparison / ) {
	    print STDERR "\n\tIn Nucleotide Level Comparisions\n"
		if $verbose;
	    $in_utr = 0;
	    $in_nuc = 1;
	}
	# F1 Score is in both UTR structure and exon structure
	elsif ( m/   F1 Score\.+(.*)/ ) {
	    my $f1_score = $1;
	    if ($in_cds) {
		$cds_f1_score = $f1_score;
		print STDERR "\t\tCDS F1 Score: ".$cds_f1_score."\n"
		    if $verbose;
	    }
	    elsif ($in_exon) {
		$exon_f1_score = $f1_score;
		print STDERR "\t\tExon F1 Score: ".$exon_f1_score."\n"
		    if $verbose;
	    }
	    elsif ($in_utr) {
		$utr_f1_score = $f1_score;
		print STDERR "\t\tUTR F1 Score: ".$utr_f1_score."\n"
		    if $verbose;
	    }
#	    elsif ($in_nuc) {
#		$cds_nuc_f1 = substr ($_, 35, 5) ;
#		$utr_nuc_f1 = substr ($_, 48, 5);
#		$overall_nuc_f1 = substr ($_, 61, 5);
#		print STDERR "\t\tCDS Nuc spec: ".$cds_nuc_f1."\n"
#		    if $verbose;
#		print STDERR "\t\tUTR Nuc spec: ".$utr_nuc_f1."\n"
#		    if $verbose;
#		print STDERR "\t\tOverall Nuc spec: ".
#		    $overall_nuc_f1."\n"
#		    if $verbose;
#	    }
	}
	elsif ( m/    Annotation edit distance\.+(.*)/ ) { 
	    my $aed = $1;
	    if ($in_cds) {
		$cds_aed = $aed;
		print STDERR "\t\tCDS AED: ".$cds_aed."\n"
		    if $verbose;
	    }
	    elsif ($in_exon) {
		$exon_aed = $aed;
		print STDERR "\t\tExon AED: ".$exon_aed."\n"
		    if $verbose;
	    }
	    elsif ($in_utr) {
		$utr_aed = $aed;
		print STDERR "\t\tUTR AED: ".$utr_aed."\n"
		    if $verbose;
	    }
	}
	elsif ( m/    Sensitivity\.+(.*)/ ) { 
	    my $sensitivity = $1;
	    if ($in_cds) {
		$cds_sensitivity = $sensitivity;
		print STDERR "\t\tCDS sensitivity: ".$cds_sensitivity."\n"
		    if $verbose;
	    }
	    elsif ($in_exon) {
		$exon_sensitivity = $sensitivity;
		print STDERR "\t\tExon sensitivity: ".$exon_sensitivity."\n"
		    if $verbose;
	    }
	    elsif ($in_utr) {
		$utr_sensitivity =$sensitivity;
		print STDERR "\t\tUTR sensitivity: ".$utr_sensitivity."\n"
		    if $verbose;
	    }
	}
	elsif ( m/    Specificity\.+(.*)/ ) { 
	    my $specificity = $1;
	    if ($in_cds) {
		$cds_specificity = $specificity;
		print STDERR "\t\tCDS specificity: ".$cds_specificity."\n"
		    if $verbose;
	    }
	    elsif ($in_exon) {
		$exon_specificity = $specificity;
		print STDERR "\t\tExon sensitivity: ".$exon_specificity ."\n"
		    if $verbose;
	    }
	    elsif ($in_utr) {
		$utr_specificity = $specificity ;
		print STDERR "\t\tUTR sensitivity: ".$utr_specificity."\n"
		    if $verbose;
	    }
	}
	elsif ( m/    Matching coefficient:/) {
	    $cds_match_coef = substr ($_, 35, 5) ;
	    $utr_match_coef = substr ($_, 48, 5);
	    $overall_match_coef = substr ($_, 61, 5);
	    print STDERR "\t\tCDS match: ".$cds_match_coef."\n"
		if $verbose;
	    print STDERR "\t\tUTR match: ".$utr_match_coef."\n"
		if $verbose;
	    print STDERR "\t\tOverall match: ".$overall_match_coef."\n"
		if $verbose;
	}
	elsif ( m/    Correlation coefficient:/) {
	    $cds_cor_coef = substr ($_, 35, 5) ;
	    $utr_cor_coef = substr ($_, 48, 5);
	    $overall_cor_coef = substr ($_, 61, 5);
	    print STDERR "\t\tCDS corr: ".$cds_cor_coef."\n"
		if $verbose;
	    print STDERR "\t\tUTR corr: ".$utr_cor_coef."\n"
		if $verbose;
	    print STDERR "\t\tOverall corr: ".$overall_cor_coef."\n"
		if $verbose;
	}
	elsif ( m/    Sensitivity:/ ) { 
	    if ($in_nuc) {
		$cds_nuc_sensitivity = substr ($_, 35, 5) ;
		$utr_nuc_sensitivity = substr ($_, 48, 5);
		$overall_nuc_sensitivity = substr ($_, 61, 5);
		print STDERR "\t\tCDS Nuc sens: ".$cds_nuc_sensitivity."\n"
		    if $verbose;
		print STDERR "\t\tUTR Nuc sens: ".$utr_nuc_sensitivity."\n"
		    if $verbose;
		print STDERR "\t\tOverall Nuc sens: ".
		    $overall_nuc_sensitivity."\n"
		    if $verbose;
	    }
	}
	elsif ( m/    Specificity:/ ) { 
	    if ($in_nuc) {
		$cds_nuc_specificity = substr ($_, 35, 5) ;
		$utr_nuc_specificity = substr ($_, 48, 5);
		$overall_nuc_specificity = substr ($_, 61, 5);
		print STDERR "\t\tCDS Nuc spec: ".$cds_nuc_specificity."\n"
		    if $verbose;
		print STDERR "\t\tUTR Nuc spec: ".$utr_nuc_specificity."\n"
		    if $verbose;
		print STDERR "\t\tOverall Nuc spec: ".
		    $overall_nuc_specificity."\n"
		    if $verbose;
	    }
	}
	elsif (m/    F1 Score:/) {
	    if ($in_nuc) {
		$cds_nuc_f1 = substr ($_, 35, 5);
		$utr_nuc_f1 = substr ($_, 48, 5);
		$overall_nuc_f1 = substr ($_, 61, 5);
		print STDERR "\t\tCDS Nuc F1: ".$cds_nuc_f1."\n"
		    if $verbose;
		print STDERR "\t\tUTR Nuc F1: ".$utr_nuc_f1."\n"
		    if $verbose;
		print STDERR "\t\tOverall Nuc F1: ".
		    $overall_nuc_f1."\n"
		    if $verbose;
	    }
	}
	elsif (m/    Annotation edit distance:/) {
	    if ($in_nuc) {
		$cds_nuc_aed = substr ($_, 35, 5);
		$utr_nuc_aed = substr ($_, 48, 5);
		$overall_aed = substr ($_, 61, 5);
		print STDERR "\t\tCDS Nuc AED: ".$cds_nuc_aed."\n"
		    if $verbose;
		print STDERR "\t\tUTR Nuc AED: ".$utr_nuc_aed."\n"
		    if $verbose;
		print STDERR "\t\tOverall AED: ".
		    $overall_aed."\n"
		    if $verbose;
	    }

	}
    }

    close(INFILE);

    


    # Write to outfile
print OUTFILE $ref_ann."\t".
    $pred_ann."\t".
    $num_genes_shared."\t".
    $num_genes_u_ref."\t".
    $num_genes_u_pred."\t".
    $num_perfect_matches."\t".
    $percent_perfect_matches."\t".
    $num_non_matches."\t".
    $percent_non_matches."\t".
    $perfect_match_len."\t".
    $non_match_len."\t".
    $perfect_match_num_pred_exons."\t".
    $perfect_match_ref_num_exons."\t".
    $non_match_num_pred_exons."\t".
    $non_match_ref_num_exons."\t".
    $cds_sensitivity."\t".
    $cds_specificity."\t".
    $cds_f1_score."\t".
    $cds_aed."\t".
    $exon_sensitivity."\t".
    $exon_specificity."\t".
    $exon_f1_score."\t".
    $exon_aed."\t".
    $cds_match_coef."\t".
    $utr_match_coef."\t".
    $overall_match_coef."\t".
    $cds_nuc_specificity."\t".
    $cds_nuc_sensitivity."\t".
    $cds_nuc_f1."\t".
    $cds_nuc_aed.
    "\n";
}

close (OUTFILE);

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

Name.pl - Short program description. 

=head1 VERSION

This documentation refers to program version 0.1

=head1 SYNOPSIS

  USAGE:
    Name.pl -i InDir -o OutDir

    --infile        # Path to the input file
    --outfie        # Path to the output file

=head1 DESCRIPTION

This is what the program does

=head1 REQUIRED ARGUMENTS

=over 2

=item -i,--indir

Path of the directory containing the sequences to process.

=item -o,--outdir

Path of the directory to place the program output.

=item -c, --config

Path to a config file. This is a tab delimited text file
indicating the required information for each of the databases to blast
against. Lines beginning with # are ignored.

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

=item --verbose

Run the program with maximum output.

=item -q,--quiet

Run the program with minimal output.

=item --test

Run the program without doing the system commands. This will
test for the existence of input files.

=back

=head1 Additional Options

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

Error messages generated by this program and possible solutions are listed
below.

=over 2

=item ERROR: No fasta files were found in the input directory

The input directory does not contain fasta files in the expected format.
This could happen because you gave an incorrect path or because your sequence 
files do not have the expected *.fasta extension in the file name.

=item ERROR: Could not create the output directory

The output directory could not be created at the path you specified. 
This could be do to the fact that the directory that you are trying
to place your base directory in does not exist, or because you do not
have write permission to the directory you want to place your file in.

=back

=head1 CONFIGURATION AND ENVIRONMENT

=head2 Configuration File

The location of the configuration file is indicated by the --config option
at the command line.
This is a tab delimited text file
indicating required information for each of the databases to blast
against. Lines beginning with # are ignored, and data are in six 
columns as shown below:

=over 2

=item Col 1. Blast program to use [ tblastx | blastn | blastx ]

The blastall program to use. DAWG-PAWS will support blastn,
tblastx, and blastx format.

=item Col 2. Extension to add to blast output file. (ie. bln )

This is the suffix which will be added to the end of your blast
output file. You can use this option to set different extensions
for different types of blast. For example *.bln for blastn
output and *.blx for blastx output.

=back

An example config file:

 #-----------------------------+
 # BLASTN: TIGR GIs            |
 #-----------------------------+
 blastn	bln	8	1e-5	TaGI_10	-a 2 -U
 blastn	bln	8	1e-5	AtGI_13	-a 2 -U
 blastn	bln	8	1e-5	ZmGI_17	-a 2 -U
 #-----------------------------+
 # TBLASTX: TIGR GIs           |
 #-----------------------------+
 tblastx	blx	8	1e-5	TaGI_10	-a 2 -U
 tblastx	blx	8	1e-5	AtGI_13	-a 2 -U
 tblastx	blx	8	1e-5	ZmGI_17	-a 2 -U

=head2 Environment

This program does not make use of variables in the user environment.

=head1 DEPENDENCIES

=head2 Required Software

=over

=item * Software Name

Any required software will be listed here.

=back

=head2 Required Perl Modules

=over

=item * Getopt::Long

This module is required to accept options at the command line.

=back

=head1 BUGS AND LIMITATIONS

Any known bugs and limitations will be listed here.

=head2 Bugs

=over 2

=item * No bugs currently known 

If you find a bug with this software, file a bug report on the DAWG-PAWS
Sourceforge website: http://sourceforge.net/tracker/?group_id=204962

=back

=head2 Limitations

=over

=item * Known Limitation

If this program has known limitations they will be listed here.

=back

=head1 SEE ALSO

The program is part of the DAWG-PAWS package of genome
annotation programs. See the DAWG-PAWS web page 
( http://dawgpaws.sourceforge.net/ )
or the Sourceforge project page 
( http://sourceforge.net/projects/dawgpaws ) 
for additional information about this package.

=head1 LICENSE

GNU General Public License, Version 3

L<http://www.gnu.org/licenses/gpl.html>

THIS SOFTWARE COMES AS IS, WITHOUT ANY EXPRESS OR IMPLIED
WARRANTY. USE AT YOUR OWN RISK.

=head1 AUTHOR

James C. Estill E<lt>JamesEstill at gmail.comE<gt>

=head1 HISTORY

STARTED:

UPDATED:

VERSION: $Rev$

=cut

#-----------------------------------------------------------+
# HISTORY                                                   |
#-----------------------------------------------------------+
#
