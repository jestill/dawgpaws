#!/usr/bin/perl -w
#-----------------------------------------------------------+
#                                                           |
# cnv_bsseeker2wiggle.pl - Convert bsseeker to wiggle       |
#                                                           |
#-----------------------------------------------------------+
#                                                           |
#  AUTHOR: James C. Estill                                  |
# CONTACT: JamesEstill_@_gmail.com                          |
# STARTED: 05/16/2013                                       |
# UPDATED: 06/18/2013                                       |
#                                                           |
# DESCRIPTION:                                              |
#  ASSUMES ALL SCAFFOLS HAVE UNIQUE NUMERICAL IDS THAT CAN  |
#  BE COERCED INTO INTEGERS !!!!!!!                         |
#  Conversion to WIGGLE is a modification of a previous     |
#  output format.                                           |
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
my ($VERSION) = q$Rev$ =~ /(\d+)/ || "pre-release";
# Get GFF version from environment, GFF2 is DEFAULT
my $gff_ver = uc($ENV{DP_GFF}) || "GFF2";

my $start_time = time;

#-----------------------------+
# VARIABLE SCOPE              |
#-----------------------------+
my $infile;
my $outdir;                    # Parent dir to hold processed output files
my $summary_out;               # Summary output is optional, tab 
                               # Delimited text file suitable for
                               # Use in R analysis

# BOOLEANS
my $quiet = 0;
my $verbose = 0;
my $show_help = 0;
my $show_usage = 0;
my $show_man = 0;
my $show_version = 0;
my $do_test = 0;                  # Run the program in test mode

my $search_scaffold;
my $bin_size = 100000;
my $do_circos = 0;
my $scaffold_prefix = "scaffold";
# Will have to get max bin location from the file


# Index vars
my $i = 0;
my $j = 0;

#-----------------------------+
# COMMAND LINE OPTIONS        |
#-----------------------------+
my $ok = GetOptions(# REQUIRED OPTIONS
		    "i|infile=s"  => \$infile,
                    "o|outdir=s"  => \$outdir,
		    "summary=s"   => \$summary_out,
		    # ADDITIONAL OPTIONS
		    "bin-size=s"  => \$bin_size,
		    "scaffold=s"  => \$search_scaffold,
		    "prefix=s"    => \$scaffold_prefix,
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
    print "\nbatch_mask.pl:\n".
	"Version: $VERSION\n\n";
    exit;
}


#-----------------------------+
# FILEHANDLES                 |
#-----------------------------+

if ($summary_out) {
    open (SUMOUT, ">$summary_out") ||
	die "Can not open summmary out: $summary_out";
    # Write summary file headers
    print SUMOUT "scaffold\tstart\tend\t".
	"X_count\t".
	"xX_Count\t".
	"X_ratio\t".
	"Y_count\t".
	"yY_count\t".
	"Y_ratio\t".
	"Z_count\t".
	"zZ_count\t".
	"Z_ratio\n";
}

if ($infile) {
    open (INFILE, "<$infile") ||
	die "Can not open infile: $infile";
}
else {
    open (INFILE, "<&STDIN") ||
	die "Can not opn STDIN for intput";
}


unless ($outdir =~ /\/$/ ) {
    $outdir = $outdir."/";
}

unless (-e $outdir) {
    print "Creating output dir ...\n" if $verbose;
    mkdir $outdir ||
	die "Could not create the output directory:\n$outdir";
}

#
# Create filehands for various output files
#
my $cov_out = $outdir.$scaffold_prefix."_bs_coverage.txt";
open (COVOUT, ">$cov_out" ) ||
    die "Could not open $cov_out for output";

my $gc_count_out = $outdir.$scaffold_prefix."_bs_gc_count.txt";
open (GCOUT, ">$gc_count_out" ) ||
    die "Could not open $gc_count_out for output";

my $metgc_count_out = $outdir.$scaffold_prefix."_bs_metgc_count.txt";
open (METGCOUT, ">$metgc_count_out" ) ||
    die "Could not open $metgc_count_out for output";


my $xX_count_out = $outdir.$scaffold_prefix."_bs_CpG_site_count.txt";
open (XXOUT, ">$xX_count_out" ) ||
    die "Could not open $xX_count_out for output";

my $X_count_out = $outdir.$scaffold_prefix."_bs_mCpG_count.txt";
open (XOUT, ">$X_count_out" ) ||
    die "Could not open $X_count_out for output";

my $X_ratio_out = $outdir.$scaffold_prefix."_bs_mCpG_ratio.txt";
open (XROUT, ">$X_ratio_out" ) ||
    die "Could not open $X_ratio_out for output";

my $yY_count_out = $outdir.$scaffold_prefix."_bs_CpHpG_site_count.txt";
open (YYOUT, ">$yY_count_out" ) ||
    die "Could not open $yY_count_out for output";

my $Y_count_out = $outdir.$scaffold_prefix."_bs_mCpHpG_count.txt";
open (YOUT, ">$Y_count_out" ) ||
    die "Could not open $Y_count_out for output";

my $Y_ratio_out = $outdir.$scaffold_prefix."_bs_mCpHpG_ratio.txt";
open (YROUT, ">$Y_ratio_out" ) ||
    die "Could not open $Y_ratio_out for output";

my $zZ_count_out = $outdir.$scaffold_prefix."_bs_CpHpH_site_count.txt";
open (ZZOUT, ">$zZ_count_out" ) ||
    die "Could not open $zZ_count_out for output";

my $Z_count_out = $outdir.$scaffold_prefix."_bs_mCpHpH_count.txt";
open (ZOUT, ">$Z_count_out" ) ||
    die "Could not open $Z_count_out for output";

my $Z_ratio_out = $outdir.$scaffold_prefix."_bs_mCpHpH_ratio.txt";
open (ZROUT, ">$Z_ratio_out" ) ||
    die "Could not open $Z_ratio_out for output";


# Global counts to include in final report
my $gc_x = 0;
my $gc_X = 0;
my $gc_y = 0;
my $gc_Y = 0;
my $gc_z = 0;
my $gc_Z = 0;

my $gc_line_count = 0;

# Will need to genrate a bin array
my @bin_coverage_count;  # Total mapping to bins
my @bin_gc_count;        # Total target gc
my @bin_metgc_count;     # Total methylated gc

# Broken down by the three types
my @bin_xX_count;
my @bin_X_count;    # methylated CG = mCpG
my @bin_yY_count;
my @bin_Y_count;    # methylated CHG = mCpHpG (where H = A, T, or C)
my @bin_zZ_count;
my @bin_Z_count;    # methylated CHH = mCpHpH (where H = A, T, or C)  

# Load bins with zero counts data using the max bin location
# as the place to end

#-----------------------------+
# MAIN PROGRAM BODY           |
#-----------------------------+
# Currently works on single scaffold, but using
# 2d arrays would make this work on all scaffolds at once
#  could also use arrays of hashes for text values of scaffolds

while (<INFILE>) {
    chomp;


    $gc_line_count++;
    my @bs_parts = split(/\t/, $_);

    # The following assumes positive strand only
    #my ($chrom, $start) = split (/\+/, $bs_parts[3]);
    # While the following allows for positive or negative strand
    my ($chrom, $start) = split (/[+-]/, $bs_parts[3]);
    $start = int($start);

    # ci is chromosome integer
    my $ci = int($chrom);              # ASSUMES ALL SCAFFOLS HAVE NUMBER IDS!

    my $bin_math = $start/$bin_size;

    # bi is bin integer
    my $bi = int($bin_math);

    #///////////////////////////////
    # NEED TO ADD TRIM STEP HERE
    # to take into account sequences
    # that are split among bins.
    # loook for features that need
    # to be split into to multiple
    # bins (may limit to two bins)
    # and then increment multiple
    # bins as needed. May also
    # need to consider strandedness
    # when splitting up the trimmed
    # sequence since the negative
    # strand sequence will probably
    # need to be added on the opposite
    # side of the bin.
    #\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
    
    # Get length of the appropriate bs_part

    # Get endpoint
    # this is 
    # can take int of endpoint
    # and see if endpoint does not equal $bi
    # this is $start + $len_seq
    # if endpoint is greater than then 


    print STDERR "\t\t".$chrom."\t".$start."\t".
	$bin_math."\t".$bi."\n"
	if $verbose;

    #X = methylated CG = mCpG
    #Y = methylated CHG = mCpHpG (where H = A, T, or C)
    #Z = methylated CHH = mCpHpH (where H = A, T, or C)
    #x = UNmethylated CG
    #y = UNmethylated CG
    #z = UNmethylated CG

    # count x,y,z and X,Y,Z characters
    my $c_x = ($bs_parts[6] =~ tr/x//);
    my $c_y = ($bs_parts[6] =~ tr/y//);
    my $c_z = ($bs_parts[6] =~ tr/z//);
    my $c_X = ($bs_parts[6] =~ tr/X//);
    my $c_Y = ($bs_parts[6] =~ tr/Y//);
    my $c_Z = ($bs_parts[6] =~ tr/Z//);


    # Increment global counters
    $gc_x = $gc_x + $c_x;
    $gc_X = $gc_X + $c_X;
    $gc_y = $gc_y + $c_y;
    $gc_Y = $gc_Y + $c_Y;
    $gc_z = $gc_z + $c_z;
    $gc_Z = $gc_Z + $c_Z;
    
    my $total_cg = $c_x + $c_y + $c_z + $c_X + $c_Y + $c_Z;
    my $methyl_cg = $c_X + $c_Y + $c_Z;

    # Increment bin coverage count
    if ( $bin_coverage_count[$ci][$bi] ) {
	$bin_coverage_count[$ci][$bi] = $bin_coverage_count[$ci][$bi] + 1;
    }
    else { 
	# First time in tin this bin set to one
	$bin_coverage_count[$ci][$bi] = 1;
    }

    # Increment gc count
    if ( $bin_gc_count[$ci][$bi] ) {
	$bin_gc_count[$ci][$bi] = $bin_gc_count[$ci][$bi] + $total_cg;
    }
    else { 
	# First time in tin this bin set to one
	$bin_gc_count[$ci][$bi] = $total_cg;
    }

    # Increment methyl gc count
    if ( $bin_metgc_count[$ci][$bi] ) {
	$bin_metgc_count[$ci][$bi] = $bin_metgc_count[$ci][$bi] + $methyl_cg;
    }
    else { 
	# First time in tin this bin set to one
	$bin_metgc_count[$ci][$bi] = $methyl_cg;
    }
    
    # Increment  @bin_xX_count
    if ( $bin_xX_count[$ci][$bi] ) {
	$bin_xX_count[$ci][$bi] = $bin_xX_count[$ci][$bi] + $c_x + $c_X;
    }
    else { 
	# First time in tin this bin set to one
	$bin_xX_count[$ci][$bi] = $c_x + $c_X;
    } 

    # Increment @bin_X_count
    if ( $bin_X_count[$ci][$bi] ) {
	$bin_X_count[$ci][$bi] = $bin_X_count[$ci][$bi] + $c_X;
    }
    else { 
	# First time in tin this bin set to one
	$bin_X_count[$ci][$bi] = $c_X;
    } 

    # incrment  @bin_yY_count
    if ( $bin_yY_count[$ci][$bi] ) {
	$bin_yY_count[$ci][$bi] = $bin_yY_count[$ci][$bi] + $c_y + $c_Y;
    }
    else { 
	# First time in tin this bin set to one
	$bin_yY_count[$ci][$bi] = $c_y + $c_Y;
    }
    
    # increment @bin_Y_count
    if ( $bin_Y_count[$ci][$bi] ) {
	$bin_Y_count[$ci][$bi] = $bin_Y_count[$ci][$bi] + $c_Y;
    }
    else { 
	# First time in tin this bin set to one
	$bin_Y_count[$ci][$bi] = $c_Y;
    }
    
    # Increment @bin_zZ_count
    if ( $bin_zZ_count[$ci][$bi] ) {
	$bin_zZ_count[$ci][$bi] = $bin_zZ_count[$ci][$bi] + $c_z + $c_Z;
    }
    else { 
	# First time in tin this bin set to one
	$bin_zZ_count[$ci][$bi] = $c_z + $c_Z;
    }

    # Increment @bin_Z_count
    if ( $bin_Z_count[$ci][$bi] ) {
	$bin_Z_count[$ci][$bi] = $bin_Z_count[$ci][$bi] + $c_Z;
    }
    else { 
	# First time in tin this bin set to one
	$bin_Z_count[$ci][$bi] = $c_Z;
    }


}


#-----------------------------+
# COVERAGE OUTPUT             |
#-----------------------------+
# This assumes that every bin contains some data

my $scaffold_num = 0;
for my $bc_ref (@bin_coverage_count) {
    # You can iterate through the array:
    
    # Get scaffold num as formatted for circos
    my $pad_s = sprintf( "%05d", $scaffold_num );
    
    my $out_i = 0;
    for my $element (@$bc_ref) {
	
	my $out_bin_start = $out_i * $bin_size + 1;
	my $out_bin_end = $out_bin_start + $bin_size - 1;
	
	my $sc_id = $scaffold_prefix.$pad_s;
	# the following is the coverage of matched strings
	print COVOUT $sc_id."\t".
	    $out_bin_start."\t".
	    $out_bin_end."\t".
	    $element."\n";
	
	$out_i++;
    }

    $scaffold_num++;

}


#-----------------------------+
# GC COVERAGE OUTPUT          |
#-----------------------------+
# Reset scaffold number to zero
$scaffold_num = 0;
for my $bc_ref (@bin_gc_count) {
    # You can iterate through the array:
    
    # Get scaffold num as formatted for circos
    my $pad_s = sprintf( "%05d", $scaffold_num );
    
    my $out_i = 0;
    for my $element (@$bc_ref) {
	
	my $out_bin_start = $out_i * $bin_size + 1;
	my $out_bin_end = $out_bin_start + $bin_size - 1;
	
	my $sc_id = $scaffold_prefix.$pad_s;
	# the following is the coverage of matched strings
	print GCOUT $sc_id."\t".
	    $out_bin_start."\t".
	    $out_bin_end."\t".
	    $element."\n";
	
	$out_i++;
    }

    $scaffold_num++;

}


#-----------------------------+
# METHYLGC COVERAGE OUTPUT    |
#-----------------------------+
# Reset scaffold number to zero
$scaffold_num = 0;
for my $bc_ref (@bin_metgc_count) {
    # You can iterate through the array:
    
    # Get scaffold num as formatted for circos
    my $pad_s = sprintf( "%05d", $scaffold_num );
    
    my $out_i = 0;
    for my $element (@$bc_ref) {
	
	my $out_bin_start = $out_i * $bin_size + 1;
	my $out_bin_end = $out_bin_start + $bin_size - 1;
	
	my $sc_id = $scaffold_prefix.$pad_s;
	# the following is the coverage of matched strings
	print METGCOUT $sc_id."\t".
	    $out_bin_start."\t".
	    $out_bin_end."\t".
	    $element."\n";
	
	$out_i++;
    }

    $scaffold_num++;

}



#-----------------------------+
# INDIVIDUAL SITE CATEGORIES  |
# AND RATIOS                  |
#-----------------------------+
# 
# Reset scaffold number to zero
$scaffold_num = 0;
for my $bc_ref (@bin_metgc_count) {
    # You can iterate through the array:
    
    # Get scaffold num as formatted for circos
    my $pad_s = sprintf( "%05d", $scaffold_num );
    
    print STDERR "Processing scaffold:\t".$pad_s."\n";
    
    my $out_i = 0;
    for my $element (@$bc_ref) {

	my $out_bin_start = $out_i * $bin_size + 1;
	my $out_bin_end = $out_bin_start + $bin_size - 1;
	
	print STDERR "\t\tProcessing array ref:".$out_i.
	    "\tBin-ID:".$out_bin_start."--".$out_bin_end."\t";
	
	my $sc_id = $scaffold_prefix.$pad_s;

	#-----------------------------+
	# X Ratio
	#-----------------------------+
	# mCpG
	my $cur_X_ratio = 0;
	my $cur_X_count = 0;
	my $cur_xX_count = 0;
	# Only do the math if the numerator has a value
	# Since we generated the denominator above, we do not
	if ( $bin_X_count[$scaffold_num][$out_i] ) {
	    if ( $bin_xX_count[$scaffold_num][$out_i] ) {
		$cur_X_ratio =  $bin_X_count[$scaffold_num][$out_i] /
		    $bin_xX_count[$scaffold_num][$out_i];
		$cur_X_count = $bin_X_count[$scaffold_num][$out_i];
		$cur_xX_count = $bin_xX_count[$scaffold_num][$out_i];
	    }
	    
	}
	$cur_X_ratio = sprintf( "%.5f", $cur_X_ratio );
	print STDERR $cur_X_ratio."\t";
	print XOUT $sc_id."\t".
	    $out_bin_start."\t".
	    $out_bin_end."\t".
	    $cur_X_count."\n";
	print XXOUT $sc_id."\t".
	    $out_bin_start."\t".
	    $out_bin_end."\t".
	    $cur_xX_count."\n";
	print XROUT $sc_id."\t".
	    $out_bin_start."\t".
	    $out_bin_end."\t".
	    $cur_X_ratio."\n";

	#-----------------------------+
	# Y Ratio
	#-----------------------------+
	# mCpHpG
	my $cur_Y_ratio = 0;
	my $cur_Y_count = 0;
	my $cur_yY_count = 0;
	# Only do the math if the numerator has a value
	# Since we generated the denominator above, we do not
	if ( $bin_Y_count[$scaffold_num][$out_i] ) {
	    if ( $bin_yY_count[$scaffold_num][$out_i] ) {
		$cur_Y_ratio =  $bin_Y_count[$scaffold_num][$out_i] /
		    $bin_yY_count[$scaffold_num][$out_i];
		$cur_Y_count =  $bin_Y_count[$scaffold_num][$out_i];
		$cur_yY_count = $bin_yY_count[$scaffold_num][$out_i];
	    }
	}
	$cur_Y_ratio = sprintf( "%.5f", $cur_Y_ratio );
	print STDERR $cur_Y_ratio."\t";

	print YOUT $sc_id."\t".
	    $out_bin_start."\t".
	    $out_bin_end."\t".
	    $cur_Y_count."\n";
	print YYOUT $sc_id."\t".
	    $out_bin_start."\t".
	    $out_bin_end."\t".
	    $cur_yY_count."\n";
	print YROUT $sc_id."\t".
	    $out_bin_start."\t".
	    $out_bin_end."\t".
	    $cur_Y_ratio."\n";

	#-----------------------------+
	# Z Ratio
	#-----------------------------+
	# mCpHpH
	my $cur_Z_ratio = 0;
	my $cur_Z_count = 0;
	my $cur_zZ_count = 0;
	# Only do the math if the numerator has a value
	# Since we generated the denominator above, we do not
	if ( $bin_Z_count[$scaffold_num][$out_i] ) {
	    if ( $bin_zZ_count[$scaffold_num][$out_i] ) {
		$cur_Z_ratio =  $bin_Z_count[$scaffold_num][$out_i] /
		    $bin_zZ_count[$scaffold_num][$out_i];
		$cur_Z_count = $bin_Z_count[$scaffold_num][$out_i];
		$cur_zZ_count = $bin_zZ_count[$scaffold_num][$out_i];
	    }
	    
	}
	$cur_Z_ratio = sprintf( "%.5f", $cur_Z_ratio );

	print STDERR $cur_Z_ratio."\t";	

	print ZOUT $sc_id."\t".
	    $out_bin_start."\t".
	    $out_bin_end."\t".
	    $cur_Z_count."\n";
	print ZZOUT $sc_id."\t".
	    $out_bin_start."\t".
	    $out_bin_end."\t".
	    $cur_zZ_count."\n";
	print ZROUT $sc_id."\t".
	    $out_bin_start."\t".
	    $out_bin_end."\t".
	    $cur_Z_ratio."\n";



	#-----------------------------+
	# PRINT SUMMARY OUTPUT        |
	#-----------------------------+
	if ($summary_out) {
	    print SUMOUT $sc_id."\t".
		$out_bin_start."\t".
		$out_bin_end."\t".
		$cur_X_count."\t".
		$cur_xX_count."\t".
		$cur_X_ratio."\t".
		$cur_Y_count."\t".
		$cur_yY_count."\t".
		$cur_Y_ratio."\t".
		$cur_Z_count."\t".
		$cur_zZ_count."\t".
		$cur_Z_ratio."\n";
	}



	print STDERR "\n";

	$out_i++;
    }

    $scaffold_num++;

}


print STDERR "\n";
print STDERR "Mapped Reads Processed:\t".$gc_line_count."\n";

my $count_obs_scaffolds = @bin_coverage_count;
$count_obs_scaffolds = $count_obs_scaffolds - 1;
# There may be fewer scaffolds than this, assumes that
# Each scaffold has a unique and incremented id
print STDERR "Highest scaffold number id: ".$count_obs_scaffolds."\n";
print STDERR "\n";
print STDERR "uCpG sites:\t".$gc_x."\n";
print STDERR "mCpG sites:\t".$gc_X."\n";
my $gc_X_ratio = $gc_X/($gc_x + $gc_X);
$gc_X_ratio = sprintf( "%.5f", $gc_X_ratio );
print STDERR "mCpG ratio:\t".$gc_X_ratio."\n";
print STDERR "\n";

print STDERR "uCpHpG sites:\t".$gc_y."\n";
print STDERR "mCpHpG sites:\t".$gc_Y."\n";
my $gc_Y_ratio = $gc_Y/($gc_y + $gc_Y);
$gc_Y_ratio = sprintf( "%.5f", $gc_Y_ratio );
print STDERR "mCpHpG ratio:\t".$gc_Y_ratio."\n";
print STDERR "\n";

print STDERR "uCpHpH sites:\t".$gc_z."\n";
print STDERR "mCpHpH sites:\t".$gc_Z."\n";
my $gc_Z_ratio = $gc_Z/($gc_z + $gc_Z);
$gc_Z_ratio = sprintf( "%.5f", $gc_Z_ratio );
print STDERR "mCpHpH ratio:\t".$gc_Z_ratio."\n";
print STDERR "\n";

my $end_time = time;
my $total_time = $end_time - $start_time;
print STDERR "Total Elapsed Time: ".$total_time."\n";
print STDERR "\n";

# Test print the coverage count

close (INFILE);
close (COVOUT);
close (GCOUT);
close (METGCOUT);

close(XXOUT);
close(XOUT);
close(XROUT);
close(YYOUT);
close(YOUT);
close(YROUT);
close(ZZOUT);
close(ZOUT);
close(ZROUT);

if ($summary_out) {
    close (SUMOUT);
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

Name.pl - Short program description. 

=head1 VERSION

This documentation refers to program version 0.1

=head1 SYNOPSIS

=head2 Usage

    Name.pl -i InFile -o OutFile

=head2 Required Arguments

    --infile        # Path to the input file
    --outfie        # Path to the output file

=head1 DESCRIPTION

This is what the program does

=head1 REQUIRED ARGUMENTS

=over 2

=item -i,--infile

Path of the input file.

=item -o,--outfile

Path of the output file.

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

This is a typcial use case.

time cnv_bsseeker2bincount.pl -i test-42.txt -o this_out --bin-size 10000 --summary summary.txt --prefix AmTr_v1.0_scaffold

=head1 DIAGNOSTICS

=over 2

=item * Expecting input from STDIN

If you see this message, it may indicate that you did not properly specify
the input sequence with -i or --infile flag. 

=back

=head1 CONFIGURATION AND ENVIRONMENT

This program does not use variables in the user environment, and does
not make use of a configuration file.

=head1 DEPENDENCIES

Other modules or software that the program is dependent on.

=head1 BUGS AND LIMITATIONS

Assumes that all sequence identifiers are unique numbers that can
be coerced into integers. If this is a problem, the program
can be rewritten to use a hash table instead of an array.


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

STARTED:

UPDATED:

VERSION: $Rev$

=cut

#-----------------------------------------------------------+
# HISTORY                                                   |
#-----------------------------------------------------------+
#
