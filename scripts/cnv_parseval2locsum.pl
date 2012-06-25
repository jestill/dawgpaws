#!/usr/bin/perl -w
#-----------------------------------------------------------+
#                                                           |
# cnv_parseval2locsum.pl - Convert parseval to locus summary|
#                                                           |
#-----------------------------------------------------------+
#                                                           |
#  AUTHOR: James C. Estill                                  |
# CONTACT: JamesEstill_@_gmail.com                          |
# STARTED: 06/18/2012                                       |
# UPDATED: 06/21/2012                                       |
#                                                           |
# DESCRIPTION:                                              |
#  Convert the parseval result to a locus summary in a      |
#  format that can be parsed and loaded to a database.      |
#  Initially this will be a tab delimited text file or      |
#  set of tab delimited text files. Providing the -o option |
#  will create multiple outfiles named with root, otherwise |
#  all output will be written to STDOUT.                    |
#                                                           |
# USAGE:                                                    |
#  cnv_parseval2locus.pl -i parsetout.txt -o out_root       |
#                                                           |
# VERSION: $Rev$                                            |
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

#-----------------------------+
# COMMAND LINE OPTIONS        |
#-----------------------------+
my $ok = GetOptions(# REQUIRED OPTIONS
		    "i|infile=s"  => \$infile,
                    "o|outfile=s" => \$outfile,
		    # ADDITIONAL OPTIONS
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
    print "\ncnv_parseval2locsum.pl:\n".
	"Version: $VERSION\n\n";
    exit;
}

#-----------------------------+
# MAIN PROGRAM BODY           |
#-----------------------------+

if ($infile) {
    open(INFILE, "<$infile") ||
	die "Can not open infile $infile";
}
else {
    print STDERR "Expecting input from STDIN\n";
    open (INFILE, "<&STDIN") ||
	die "Can not open STDIN for input";
}


# If outfile root name given write outfiles
# otherwise sending evertyhing to STDOUT
if ($outfile) {
    my $tab_out = $outfile."_transcript_vals.txt";
    my $trans_out = $outfile."_transcript_list.txt";
    my $loc_out = $outfile."_locus_list.txt";
    my $loc_val = $outfile."_locus_vals.txt";

    open (TABOUT, ">$tab_out") ||
	die "Can not open outfile $tab_out";
    open (TRANSOUT, ">$trans_out") ||
	die "Can not open outfile $trans_out";
    open (LOCOUT, ">$loc_out") ||
	die "Can not open outfile $loc_out";
    open (LOCVAL, ">$loc_val") ||
	die "Can not open outfile $loc_val";

}
else {
    open (TABOUT, ">&STDOUT") ||
	die "Can not open STDOUT for output\n";
}

# Print header to tabout here

# BOOLEANS
my $in_locus;
my $in_reference_gene;
my $in_prediction_gene;
my $in_reference_transcript;
my $in_prediction_transcript;
my $in_reference_gff;
my $in_prediction_gff;
my $in_cds;
my $in_exon;
my $in_utr;
my $in_nuc;
my $in_comparison;
my $end_comparison;
my $in_novel_transcript;
my $in_loc_splice_complexity;

# A novel locus is a predicted gene
# without a reference gene
my $in_novel_locus;

# Locus vars
my $seq_id;
my $start;
my $end;
my $locus;
my $ref_splice_complexity;
my $pred_splice_complexity;

# Transcript booleans
my $perfect_cds_match;
my $perfect_exon_match;
my $perfect_utr_match;
my $perfect_gene_match;

# Transcript cds vars
my $cds_sensitivity;
my $cds_specificity;
my $cds_f1_score;
my $cds_aed;

# Transcript exon vars
my $exon_sensitivity;
my $exon_specificity;
my $exon_f1_score;
my $exon_aed;

# Transcript UTR vars
my $utr_sensitivity;
my $utr_specificity;
my $utr_f1_score;
my $utr_aed;

my $nuc_cds_match_coefficient;
my $nuc_cds_cor_coefficient;
my $nuc_cds_sensitivity;
my $nuc_cds_specificity;
my $nuc_cds_f1_score;
my $nuc_cds_aed;

my $nuc_utr_match_coefficient;
my $nuc_utr_cor_coefficient;
my $nuc_utr_sensitivity;
my $nuc_utr_specificity;
my $nuc_utr_f1_score;
my $nuc_utr_aed;

my $nuc_all_match_coefficient;
my $nuc_all_cor_coefficient;
my $nuc_all_sensitivity;
my $nuc_all_specificity;
my $nuc_all_f1_score;
my $nuc_all_aed;

my @ref_transcripts;
my @pred_transcripts;
my @ref_genes; 
my @pred_genes; 

my $comp_count = 0;

while (<INFILE>) {

    chomp;

    print STDERR $_."\n" if $verbose;

    # Determine where in the file we are
    if ( $_ =~ m/Locus: sequence \'(.*)\' from (.*) to (.*)/ ) {
	print STDERR "In Locus\n"
	    if $verbose;
	$seq_id = trim($1);
	$start = trim($2);
	$end = trim($3);
	$locus = $seq_id."_".$start."_".$end;
	print STDERR "Locus: $locus\n" if $verbose;
	print STDERR "\t1: $seq_id\n" if $verbose;
	print STDERR "\t2: $start\n" if $verbose;
	print STDERR "\t3: $end\n" if $verbose;
	if ($in_novel_transcript) {
	    print STDERR "...Out of novel transcript\n"
		if $verbose;
	}
	$in_locus = 1;
	$in_novel_transcript = 0;
    }
    elsif ( $_ =~ m/reference genes/ ) {
	print STDERR "\tIn reference genes\n"
	    if $verbose;
	$in_reference_gene = 1;
    }
    elsif ( $_ =~ m/prediction genes/ ) {
	print STDERR "\tIn predictions genes\n"
	    if $verbose;
	$in_reference_gene = 0;
	$in_prediction_gene = 1;
    }
    elsif ( $_ =~ m/locus splice complexity/ ) {
	$in_prediction_gene = 0;
	$in_loc_splice_complexity = 1;
    }
    elsif ( $_ =~ m/reference transcripts/ ) {
	$in_reference_transcript = 1;
    }
    elsif ( $_ =~ m/prediction transcripts/ ) {
	$in_prediction_transcript = 1;
	$in_reference_transcript = 0;
	if ( $_ =~ m/novel prediction transcripts/ ) {
	    $in_prediction_transcript = 0;
	    $in_novel_transcript = 1;
	    print STDERR "Novel transcript\n"
		if $verbose;
	}
    }
    elsif ( $_ =~ m/novel prediction transcripts/ ) {
	$in_reference_transcript = 0;
	$in_prediction_transcript = 0;
    }
    elsif ( $_ =~ m/reference GFF3/ ) {
	$in_prediction_transcript = 0;
    }
    elsif ( $_ =~ m/prediction GFF3/ ) {

    }
    elsif ( $_ =~ m/CDS structure comparison/ ) {
	print STDERR "\tCDS structure comparison\n"
	    if $verbose;
	$in_cds = 1;
    }
    elsif ( $_ =~ m/Exon structure comparison/ ) {
	print STDERR "\tExon structure comparison\n"
	    if $verbose;
	$in_exon = 1;
	$in_cds = 0;
    }
    elsif ( $_ =~ m/UTR structure comparison/ ) {
	print STDERR "\tUTR structure comparison\n"
	    if $verbose;
	$in_utr = 1;
	$in_exon = 0;
    }
    elsif ( $_ =~ m/Nucleotide-level comparison/ ) {
	$in_nuc = 1;
	$in_utr = 0;
    }
    elsif ( $_ =~ m/No comparisons were performed for this locus/ ) {
	#$in_loc_splice_complexity = 0;
    }
    elsif ( $_ =~ m/Begin Comparison/ ) {
	#$in_loc_splice_complexity = 0;
	print STDERR "\n>>> In Comparison >>>\n\n"
	    if $verbose;
    }
    elsif ( $_ =~ m/End Comparison/ ) {
	$end_comparison = 1;
	print STDERR "\n<<< Finished Comparision <<<\n\n"
	    if $verbose;
	$in_nuc = 0;
    }

    #-----------------------------+
    # PERFECT MATCHES             |
    #-----------------------------+
    # FOR PERFECT MATCHES SET VALUES TO MAX
    if ($_ =~ m/CDS structures match perfectly/ ) {

	print STDERR "\n\tPerfect CDS match\n"
	    if $verbose;

	$perfect_cds_match = 1;

	$cds_sensitivity = "1.000";
	$cds_specificity = "1.000";
	$cds_f1_score = "1.000";
	$cds_aed = "0.000";

	$nuc_cds_match_coefficient = "1.000";
	$nuc_cds_cor_coefficient = "1.000";
	$nuc_cds_sensitivity = "1.000";
	$nuc_cds_specificity = "1.000";
	$nuc_cds_f1_score = "1.000";
	$nuc_cds_aed = "0.000";
    }

    if ($_ =~ m/Exon structures match perfectly/ ) {
	print STDERR "\tPerfect exon match\n"
	    if $verbose;

	$perfect_exon_match = 1;
	
	$exon_sensitivity = "1.000";
	$exon_specificity = "1.000";
	$exon_f1_score = "1.000";
	$exon_aed = "0.000";

    }

    # Perfect UTR structure match
    if ($_ =~ m/UTR structures match perfectly/ ) {
	print STDERR "\tPerfect UTR match\n"
	    if $verbose;
	
	$perfect_utr_match = 1;

	$utr_sensitivity = "1.000";
	$utr_specificity = "1.000";
	$utr_f1_score = "1.000";
	$utr_aed = "0.000";

	$nuc_utr_match_coefficient = "1.000";
	$nuc_utr_cor_coefficient = "1.000";
	$nuc_utr_sensitivity = "1.000";
	$nuc_utr_specificity = "1.000";
	$nuc_utr_f1_score = "1.000";
	$nuc_utr_aed = "0.000";
	
    }

    # UTRs not predicteid
    if ($_ =~ m/No UTRs annotated for this locus/ ) {
	print STDERR "NO UTRS Present\n"
	    if $verbose;

	$perfect_utr_match = 0;

	$utr_sensitivity = "--";
	$utr_specificity = "--";
	$utr_f1_score = "--";
	$utr_aed = "--";
	
	$nuc_utr_match_coefficient = "--";
	$nuc_utr_cor_coefficient = "--";
	$nuc_utr_sensitivity = "--";
	$nuc_utr_specificity = "--";
	$nuc_utr_f1_score = "--";
	$nuc_utr_aed = "--";
	
    }

    # Perfect gene structure match
    if ($_ =~ m/Gene structures match perfectly/ ) {
	print STDERR "\tPerfect gene structure match\n"
	    if $verbose;

	$perfect_gene_match = 1;

	$nuc_all_match_coefficient = "1.000";
	$nuc_all_cor_coefficient = "1.000";
	$nuc_all_sensitivity = "1.000";
	$nuc_all_specificity = "1.000";
	$nuc_all_f1_score = "1.000";
	$nuc_all_aed = "0.000";
	
    }


    #-----------------------------+
    # Work with components        |
    #-----------------------------+

    # Splice complexity
    # This also used to determine
    # when finished processing locus information
    if ( $in_loc_splice_complexity ) {
	if ($_ =~ m/^.*reference\:(.*)/ ) {
	    $ref_splice_complexity = trim($1);
	    print STDERR "\nRef--$ref_splice_complexity\n"
		if $verbose;

	}
	elsif ($_ =~ m/^.*prediction\:(.*)/ ) {
	    $pred_splice_complexity = trim($1);
	    print STDERR "Pred--$pred_splice_complexity\n"
		if $verbose;
	}
	elsif ( $_ =~ m/locus splice complexity/ ) {
	    # Do nothing this is the header
	}
	elsif ( $_ =~ m/\s|.*/) {

	    $in_loc_splice_complexity = 0;

	    my $num_pred_genes = @pred_genes;
	    my $num_ref_genes = @ref_genes;



	    foreach my $ref_gene (@ref_genes) {
		if ($outfile) {
		    print LOCOUT $locus."\t".
			"ref_gene_id\t".
			$ref_gene."\n";
		}
		else {
		    print TABOUT "## ".$locus."\t".
			"ref_gene_id\t".
			$ref_gene."\n";
		}

		if ( $ref_gene =~ "None" ) {
		    $num_ref_genes = int(0);
		}
	    }

	    foreach my $pred_gene (@pred_genes) {

		if ($outfile) {
		    print LOCOUT $locus."\t".
			"pred_gene_id\t".
			$pred_gene."\n";
		}
		else {
		    print TABOUT "## ".$locus."\t".
			"pred_gene_id\t".
			$pred_gene."\n";
		}
		if ( $pred_gene =~ "None" ) {
		    $num_pred_genes = int(0);
		}


	    }
	    

	    # Print locus summary values information
	    if ($outfile) {
		print LOCVAL $locus."\t".
		    $seq_id."\t".
		    $start."\t".
		    $end."\t".
		    $num_ref_genes."\t".
		    $num_pred_genes."\t".
		    $ref_splice_complexity."\t".
		    $pred_splice_complexity."\n";
	    }
	    else {
		print TABOUT "## ".$locus."\t".
		    $seq_id."\t".
		    $start."\t".
		    $end."\t".
		    $num_ref_genes."\t".
		    $num_pred_genes."\t".
		    $ref_splice_complexity."\t".
		    $pred_splice_complexity."\n";
	    }

	    # Reset gene arrays to zero
	    @pred_genes = ();
	    @ref_genes = ();
	    

	}
	
    } # End of in_loc_splice_complexity

    # Novel transcripts ...
    if ( $in_novel_transcript ) {
	if ($_ =~ m/\s\|\s(.*)/ ) {
	    my $novel_transcript = trim($1);
	    unless ( $_ =~ m/novel prediction transcripts/ ) {
		$comp_count++;
		print STDERR "Novel pred transcript -- $novel_transcript \n"
		    if $verbose;
		
		if ( $outfile ) {
		    print TRANSOUT $comp_count."\tpred_novel\t".
			$novel_transcript."\n";
		}
		else {
		    print TABOUT "# ".$comp_count."\tpred_novel\t".
			$novel_transcript."\n";
		}

		print TABOUT $comp_count."0\t1\n";
		# May want to add columns of null values here
		# for AED values etc ...

	    }
	}
    }

    # Predicted transcript IDs
    if ( $in_prediction_transcript ) {
	if ($_ =~ m/\s\|\s(.*)/ ) {
	    unless ( $_ =~ m/prediction transcripts/ ) {
		my $pred_transcript = trim($1);
		push ( @pred_transcripts, $pred_transcript);
		print STDERR "\tPred match:".$pred_transcript."\n"
		    if $verbose;	    
	    }
	}
    }

    # Reference transcript IDs
    if ( $in_reference_transcript ) {
	if ($_ =~ m/\s\|\s(.*)/ ) {
	    unless ($_ =~ m/reference transcripts/ ) {
		my $ref_transcript = trim($1);
		push (@ref_transcripts, $ref_transcript);
		print STDERR "\tRef match:".$ref_transcript."\n"
		    if $verbose;
	    }
	}
    }


    # Load sensitivy values
    # Use presence absence of sensitivity value
    # to set perfect match booleans
    if ($_ =~ m/Sensitivity:\s(.*)$/ ) {
	my $sensitivity = trim($1);
	$perfect_gene_match = 0;

	print STDERR "\n"
	    if $verbose;
	print STDERR "\tSensitivity:--".$sensitivity."\n"
	    if $verbose;
	if ($in_cds) {
	    $cds_sensitivity = $sensitivity;
	    $perfect_cds_match = 0;
	}
	elsif ($in_exon) {
	    $exon_sensitivity =  $sensitivity;
	    $perfect_exon_match = 0;
	}
	elsif ($in_utr) {
	    $utr_sensitivity = $sensitivity;
	    $perfect_utr_match = 0;
	}
	elsif ($in_nuc) {
	    ($nuc_cds_sensitivity,
	     $nuc_utr_sensitivity,
	     $nuc_all_sensitivity) = split (/\s+/,  $sensitivity);
	}
    }

    # Load specificity values
    if ($_ =~ m/Specificity:\s(.*)$/ ) {
	my $specificity = trim($1);
	print STDERR "\tSpecificity:--".$specificity."\n"
	    if $verbose;

	if ($in_cds) {
	    $cds_specificity = $specificity;
	}
	elsif ($in_exon) {
	    $exon_specificity = $specificity;
	}
	elsif ($in_utr) {
	    $utr_specificity = $specificity;
	}
	elsif ($in_nuc) {
	    ($nuc_cds_specificity,
	     $nuc_utr_specificity,
	     $nuc_all_specificity) = split (/\s+/, $specificity);
	}
    }

    # Load F1 Score
    if ($_ =~ m/F1 Score:\s(.*)$/ ) {
	my $f1_score = trim($1);
	print STDERR "\tF1 Score:--".$f1_score."\n"
	    if $verbose;

	if ($in_cds) {
	    $cds_f1_score = $f1_score;
	}
	elsif ($in_exon) {
	    $exon_f1_score = $f1_score;
	}
	elsif ($in_utr) {
	    $utr_f1_score = $f1_score;
	}
	elsif ($in_nuc) {
	    ($nuc_cds_f1_score,
	     $nuc_utr_f1_score,
	     $nuc_all_f1_score) = split (/\s+/, $f1_score);
	}
    }

    # Annotation edit distance
    if ($_ =~ m/Annotation edit distance:\s(.*)$/ ) {
	my $aed = trim($1);

	print STDERR "\tAED:--".$aed."\n"
	    if $verbose;

	if ($in_cds) {
	    $cds_aed = $aed;
	}
	elsif ($in_exon) {
	    $exon_aed = $aed;
	}
	elsif ($in_utr) {
	    $utr_aed = $aed;
	}
	elsif ($in_nuc) {
	    ($nuc_cds_aed,
	     $nuc_utr_aed,
	     $nuc_all_aed) = split (/\s+/, $aed);
	}
    }

    # Matching Coefficient
    if ($_ =~ m/Matching coefficient:\s(.*)$/ ) {
	my $match_coefficient = trim($1);
	($nuc_cds_match_coefficient,
	 $nuc_utr_match_coefficient,
	 $nuc_all_match_coefficient) = split (/\s+/, $match_coefficient);
    }

    #  Correlation coefficient
    if ($_ =~ m/Correlation coefficient:\s(.*)$/ ) {
	my $cor_coefficient = trim($1);
	($nuc_cds_cor_coefficient,
	 $nuc_utr_cor_coefficient,
	 $nuc_all_cor_coefficient) = split (/\s+/, $cor_coefficient);
    }

    # Get Reference Gene
    if ( $in_reference_gene ) {
	if ($_ =~ m/^\|\s(.*)/ ) {
	    unless ($_ =~ m/reference genes/ ) {
		my $ref_gene = trim($1);
		push (@ref_genes, $ref_gene);
		print STDERR "\tRef gene".$ref_gene."\n"
		    if $verbose;
		
		if ($ref_gene =~ "None") {
		    print STDERR "NO REF GENE MATCH\n"
			if $verbose;

		    # May want to load some dummy values here
		    # for the tab out table
		}


	    }
	}

    }

    if ( $in_prediction_gene ) {
	if ($_ =~ m/\|\s(.*)/ ) {
	    unless ($_ =~ m/prediction genes/ ) {
		my $pred_gene = trim($1);
		push (@pred_genes, $pred_gene);
		print STDERR "\tPred gene match:".$1."\n"
		    if $verbose;

		# /////////////
		if ($_ =~ m/None/) {
		    print STDERR "NO PRED GENE MATCH\n"
			if $verbose;
		}

	    }
	}
    }
    
    #-----------------------------+
    # Write to outfile            |
    #-----------------------------+
    # When comparisons ended .. write transcript output
    if ( $end_comparison ) {

	$comp_count++;

	my $num_pred_transcripts =  @pred_transcripts;
	my $num_ref_transcripts = @ref_transcripts;

	my $pred_transcripts;
	my $ref_transcripts;

	foreach my $transcript ( @ref_transcripts ) {
	    if ( $outfile ) {
		print TRANSOUT $comp_count."\tref_transcript_id\t".
		    $transcript."\n";
	    }
	    else {
		print TABOUT "# ".$comp_count."\tref_transcript_id\t".
		    $transcript."\n";
	    }
	}	

	
	
	foreach my $transcript ( @pred_transcripts ) {
	    if ( $outfile ) {
		print TRANSOUT $comp_count."\tpred_transcript_id\t".
		    $transcript."\n";
	    }
	    else {
		print TABOUT "# ".$comp_count."\tpred_transcrirpt_id\t".
		    $transcript."\n";
	    }
	}

	# Write the values determined
	print TABOUT $comp_count."\t";
	print TABOUT $locus."\t".
	    $seq_id."\t".
	    $start."\t".
	    $end."\t";
	print TABOUT $num_ref_transcripts."\t".
	    $num_pred_transcripts."\t";


	# Booleans as Y/N

	# Perfect CDS match
	if ( $perfect_cds_match ) {
	    print TABOUT "Y\t";
	}
	else {
	    print TABOUT "N\t";
	}

	# Perfect exon match
	if ( $perfect_exon_match ) {
	    print TABOUT "Y\t";
	}
	else {
	    print TABOUT "N\t";
	}

	# Perfect UTR match
	if ( $perfect_utr_match ) {
	    print TABOUT "Y\t";
	}
	else {
	    print TABOUT "N\t";
	}

	print TABOUT $cds_sensitivity."\t".
	    $cds_specificity."\t".
	    $cds_f1_score."\t".
	    $cds_aed."\t".
	    $nuc_cds_match_coefficient."\t".
	    $nuc_cds_cor_coefficient."\t".
	    $nuc_cds_sensitivity."\t".
	    $nuc_cds_specificity."\t".
	    $nuc_cds_f1_score."\t".
	    $nuc_cds_aed."\t".
	    # exon values
	    $exon_sensitivity."\t".
	    $exon_specificity."\t".
	    $exon_f1_score."\t".
	    $exon_aed."\t".
	    # UTR VALUES
	    $utr_sensitivity."\t".
	    $utr_specificity."\t".
	    $utr_f1_score."\t".
	    $utr_aed."\t".
	    $nuc_utr_match_coefficient."\t".
	    $nuc_utr_cor_coefficient."\t".
	    $nuc_utr_sensitivity."\t".
	    $nuc_utr_specificity."\t".
	    $nuc_utr_f1_score."\t".
	    $nuc_utr_aed."\t".
	    "\n";

	# Reset values
	@pred_transcripts = ();
	@ref_transcripts = ();
	$end_comparison = 0;

    }

}

close (INFILE);
close (LOCOUT);
close (LOCVAL);
close (TABOUT);
close (TRANSOUT);

exit 0;

#-----------------------------------------------------------+ 
# SUBFUNCTIONS                                              |
#-----------------------------------------------------------+

sub trim {
    my $string = shift;
    $string =~ s/^\s+//;
    $string =~ s/\s+$//;
    return $string;
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

1;
__END__

=head1 NAME

cnv_parseval2locsum.pl -  Convert parseval to locus summary

=head1 VERSION

This documentation refers to program version 0.1

=head1 SYNOPSIS

=head2 Usage

    cnv_parseval2locus.pl -i parsetout.txt -o out_root

=head2 Required Arguments

    --infile        # Path to the parseval result file
    --outfile       # Root name for outfiles

=head1 DESCRIPTION

Convert the parseval result to a locus summary in a
format that can be parsed and loaded to a database.
Initially this will be a tab delimited text file or
set of tab delimited text files. Providing the -o option
will create multiple outfiles named with root, otherwise
all output will be written to STDOUT.  

=head1 REQUIRED ARGUMENTS

=over 2

=item -i,--infile

Path of the input file. This is the result file parseval.

=item -o,--outfile

Root name for the outfiles produced. If no outfile root name is given,
all output is written to STDOUT.

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

=head2 Write output to output files

The --outfile or -o option is used as a root identiftier to build names for 
output files:

 cnv_parseval2locsum.pl -i athal_test.txt -o athal_out

This will produce two files. One file named root_name_vals.txt with
similarity values between prediction and reference transcripts, and
root_name_loc.txt with the list of the loci invovled. The first column
in each file is an unique integer ID that allows for matching between
the transcript pairings and the similarity values file.

=head3 Transcript Output File

The file with reference and predicted transcript pairings 
(ie athal_out_loc.txt):

 1	ref_transcript_id	AT1G01010.1
 1	pred_transcript_id	AT1G01010.1
 2	ref_transcript_id	AT1G01020.1
 2	pred_transcript_id	AT1G01020.1
 3	ref_transcript_id	AT1G01020.2
 3	pred_transcript_id	AT1G01020.2
 4	ref_transcript_id	AT1G01030.1
 4	pred_transcript_id	AT1G01030.1
 5	ref_transcript_id	AT1G01040.1
 5	pred_transcript_id	AT1G01040.1

=head3 Values Output File

The file with the (athal_out_vals.txt):

 1	1	1	Y	Y	Y	1.000	1.000	1.000
 2	1	1	Y	Y	Y	1.000	1.000	1.000
 3	1	1	Y	Y	Y	1.000	1.000	1.000
 4	1	1	Y	Y	Y	1.000	1.000	1.000
 5	1	1	Y	Y	Y	1.000	1.000	1.000

=head2 Write output to STDOUT

To write all output to the standard output stream:

 cnv_parseval2locsum.pl -i athal_test.txt

The values that would have been written to the locus pairings output file
are written to the values output file preceded by #. 

An example of produce results:

  # 1	ref_transcript_id	AT1G01010.1
  # 1	pred_transcrirpt_id	AT1G01010.1
  1	1	1	Y	Y	Y	1.000	1.000	1.000	0.000
  # 2	ref_transcript_id	AT1G01020.1
  # 2	pred_transcrirpt_id	AT1G01020.1
  2	1	1	Y	Y	Y	1.000	1.000	1.000	0.000
  # 3	ref_transcript_id	AT1G01020.2
  # 3	pred_transcrirpt_id	AT1G01020.2
  3	1	1	Y	Y	Y	1.000	1.000	1.000	0.000

=head1 DIAGNOSTICS

=over 2

=item * Expecting input from STDIN

If you see this message, it may indicate that you did not properly specify
the input sequence with -i or --infile flag. 

=back

=head1 CONFIGURATION AND ENVIRONMENT

This program does not make use of variables defined in the user environment.

=head1 DEPENDENCIES

=head2 Required Software

This program is dependent on the following program

=over

=item ParsEval

ParsEval is a program for comparing two sets of gene structure 
annotations (provided as GFF3 files) for the same genomic region. 

Source: http://parseval.sourceforge.net/

=item GenomeTools

The ParsEval program is dependend on the GenomeTools program.

Source: http://genometools.org/

=back

=head1 BUGS AND LIMITATIONS

Any known bugs and limitations will be listed here.

=head1 REFERENCE

Your use of the DAWGPAWS programs should reference the following manuscript:

JC Estill and JL Bennetzen. 2009. 
"The DAWGPAWS Pipeline for the Annotation of Genes and Transposable 
Elements in Plant Genomes." Plant Methods. 5:8.

=head1 LICENSE

GNU General Public License, Version 3

L<http://www.gnu.org/licenses/gpl.html>

=head1 AUTHOR

James C. Estill E<lt>JamesEstill at gmail.comE<gt>

=head1 HISTORY

STARTED: 06/18/2012

UPDATED: 06/22/2012

VERSION: $Rev$

=cut

#-----------------------------------------------------------+
# HISTORY                                                   |
#-----------------------------------------------------------+
#
# The following is what needs to be parsed
#     |
#     |  CDS structure comparison
#     |    8 reference CDS segments
#     |      8 match prediction
#     |      0 don't match prediction
#     |    28 prediction CDS segments
#     |      8 match reference
#     |      20 don't match reference
#     |    Sensitivity:                   1.000     
#     |    Specificity:                   0.286     
#     |    F1 Score:                      0.444     
#     |    Annotation edit distance:      0.357     
#     |
#     |  Exon structure comparison
#     |    9 reference exons
#     |      9 match prediction
#     |      0 don't match prediction
#     |    29 prediction exons
#     |      9 match reference
#     |      20 don't match reference
#     |    Sensitivity:                   1.000     
#     |    Specificity:                   0.310     
#     |    F1 Score:                      0.474     
#     |    Annotation edit distance:      0.345     
#     |
#     |  UTR structure comparison
#     |    2 reference UTRs
#     |      2 match prediction
#     |      0 don't match prediction
#     |    4 prediction UTRs
#     |      2 match reference
#     |      2 don't match reference
#     |    Sensitivity:                   1.000     
#     |    Specificity:                   0.500     
#     |    F1 Score:                      0.667     
#     |    Annotation edit distance:      0.250     
#     |
#     |  Nucleotide-level comparison      CDS          UTRs         Overall   
#     |    Matching coefficient:          0.427        0.986        0.230
#     |    Correlation coefficient:       0.197        0.831        --        
#     |    Sensitivity:                   1.000        1.000        --        
#     |    Specificity:                   0.388        0.985        --        
#     |    F1 Score:                      0.182        0.824        --        
#     |    Annotation edit distance:      0.306        0.007        --        
#     |
#     |--------------------------
#     |----- End Comparison -----
#     |--------------------------
#
# Complete matches are as follows
#     |  CDS structure comparison
#     |    7 reference CDS segments
#     |    7 prediction CDS segments
#     |    CDS structures match perfectly!
#     |
#     |  Exon structure comparison
#     |    10 reference exons
#     |    10 prediction exons
#     |    Exon structures match perfectly!
#     |
#     |  UTR structure comparison
#     |    5 reference UTRs
#     |    5 prediction UTRs
#     |    UTR structures match perfectly!
#     |
#     |  Gene structures match perfectly!
#     |
#     |--------------------------
#     |----- End Comparison -----
#     |--------------------------
#     |
#     |---------
#
# Novel prediction transcripts
#     |--------------------------
#     |----- End Comparison -----
#     |--------------------------
#     |
#     |  novel prediction transcripts (or transcript sets)
#     |    AT5G67580.2
#
#|-------------------------------------------------
#|---- Locus: sequence 'Chr5' from 26958001 to 26959557
#|-------------------------------------------------
#|
#
# It is also possible for there to be no predicted genes
# for ra locs

#|-------------------------------------------------
#|---- Locus: sequence 'AmTr_v1.0_scaffold00001' from 185581 to 198384
#|-------------------------------------------------
#|
#|  reference genes:
#|    evm_15.TU.AmTr_v1.0_scaffold00001.13
#|
#|  prediction genes:
#|    None!
#|
#|  locus splice complexity:
#|    reference:   0.000
#|    prediction:  0.000
#|
#|
#|----------
#     |
#     |  No comparisons were performed for this locus
#     |
#
#|-------------------------------------------------
#|---- Locus: sequence 'AmTr_v1.0_scaffold00001' from 201972 to 206735
#|-------------------------------------------------
#|
#|  reference genes:
#|    evm_15.TU.AmTr_v1.0_scaffold00001.14
#|
#|  prediction genes:
#|    None!
#|
#|  locus splice complexity:
#|    reference:   0.000
#|    prediction:  0.000
#|
#|
#|----------
#     |
#     |  No comparisons were performed for this locus
#     |
#
