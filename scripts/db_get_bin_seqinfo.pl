#!/usr/bin/perl -w
#-----------------------------------------------------------+
#                                                           |
# db_get_bin_seq_info.pl - Sequence base information        |
#                                                           |
#-----------------------------------------------------------+
#                                                           |
#  AUTHOR: James C. Estill                                  |
# CONTACT: JamesEstill_@_gmail.com                          |
# STARTED: 05/13/2013                                       |
# UPDATED: 05/20/2013                                       |
#                                                           |
# DESCRIPTION:                                              |
#  Given a bin size summarize sequence information for      |
#  that bin.                                                |
#  - all nucleotides                                        |
#  - all 16 dinucleotides                                   |
#  - all potential methylation sites                        |
#    as (methylation potential)                             |
#                                                           |
# USAGE:                                                    |
#                                                           |
# VERSION: $Rev$                                            |
#                                                           |
#                                                           |
# LICENSE:                                                  |
#  GNU General Public License, Version 3                    |
#  http://www.gnu.org/licenses/gpl.html                     |  
#                                                           |
#-----------------------------------------------------------+
# TO DO:
# Sub selection of each scaffold in the GFF file.
# This is using GFF file and object
# "Bio::SeqFeature::Annotated
# Methods documented at
# http://doc.bioperl.org/bioperl-live/Bio/SeqFeature/Annotated.html
# It would be better to load these to a GFF data store and use features
# to fetch the overlapping features. This will be more scalable
# http://search.cpan.org/~cjfields/BioPerl-1.6.1/Bio/DB/SeqFeature/Segment.pm
# #
#
#http://www.bioperl.org/wiki/Module:Bio::Db::SeqFeature::Store
# Methods at
# Http://search.cpan.org/~cjfields/BioPerl-1.6.1/Bio/DB/SeqFeature/Store.pm
# Example of adding new exon
# http://doc.bioperl.org/bioperl-live/Bio/DB/SeqFeature/NormalizedFeature.html
package DAWGPAWS;

#-----------------------------+
# INCLUDES                    |
#-----------------------------+
use Graph;
use strict;
use warnings;
use Getopt::Long;
# The following needed for printing help
use Pod::Select;               # Print subsections of POD documentation
use Pod::Text;                 # Print POD doc as formatted text file
use IO::Scalar;                # For print_help subfunction
use IO::Pipe;                  # Pipe for STDIN, STDOUT for POD docs
use File::Spec;                # Convert a relative path to an abosolute path
#use Set::IntervalTree;
#use Data::Range::Compare;      #
#use Bio::Range;
#use Bio::FeatureIO;

# This uses a database 
use Bio::DB::SeqFeature::Store;

#-----------------------------+
# Program VARIABLES           |
#-----------------------------+
my ($VERSION) = q$Rev$ =~ /(\d+)/ || "pre-release";
# Get GFF version from environment, GFF2 is DEFAULT
my $gff_ver = uc($ENV{DP_GFF}) || "GFF2";
my %all_sources = ();          # Sources

# The parent feature to fetch overlaps for, others 
# that may make use of this include features such as gene or transposable_element
# that may have child features you want to fetch. This may require
# moving to bioperl GFF feature objects.
# parent feature will be target DDE matches
#jamie-macbook-pro:scripts jestill$ ./db_get_target_repseek_overlap.pl --feature "match:TARGET"
# flanking feature is 
# direct_repeat:repseek_l24

#-----------------------------+
# VARIABLE SCOPE              |
#-----------------------------+
my $infile;
my $outfile;
my $summary_outfile;

# BOOLEANS
my $quiet = 0;
my $verbose = 0;
my $show_help = 0;
my $show_usage = 0;
my $show_man = 0;
my $show_version = 0;
my $do_test = 0;                  # Run the program in test mode

# DB variables
my $adaptor = 'DBI::mysql';
my $user = 'dawgpaws';
my $password ='f34tur3s';
my $host = 'localhost';
my $database = 'amtr_ltr';
my $fasta_out;
my $overlap_feature = "gene";
my $min_feature_length = 100;
#my $min_feature_length;
my $bin_length = 100000;

# GBROWSE
my $gb_url = "http://localhost/cgi-bin/gb2/gbrowse/amborella_te/?name=";

#-----------------------------+
# COMMAND LINE OPTIONS        |
#-----------------------------+
my $ok = GetOptions(# REQUIRED OPTIONS
		    "i|infile=s"   => \$infile,
                    "o|outfile=s"  => \$outfile,
		    "summary=s"    => \$summary_outfile,
#		    "f|feature=s"  => \$parent_feature,
		    "overlap=s"    => \$overlap_feature,
		    "min-length=s" => \$min_feature_length,
		    "bin-length=s" => \$bin_length,
		    "fasta=s"      => \$fasta_out, 
		    # Seq data store database options
#		    "parent"       => \$parent_feature,
		    "p|pass=s"     => \$password,
		    "u|user=s"     => \$user,
		    "host"         => \$host,
		    "adaptor"      => \$adaptor,
		    "database"     => \$database,
		    # ADDITIONAL OPTIONS
		    "gff-ver=s"    => \$gff_ver,
		    "q|quiet"      => \$quiet,
		    "verbose"      => \$verbose,
		    # ADDITIONAL INFORMATION
		    "usage"        => \$show_usage,
		    "test"         => \$do_test,
		    "version"      => \$show_version,
		    "man"          => \$show_man,
		    "h|help"       => \$show_help,);


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
    # User perldoc to generate the man documentation.    system ("perldoc $0");
    exit($ok ? 0 : 2);
}

if ($show_version) {
    print "\nsql_cluster_features.pl:\n".
	"Version: $VERSION\n\n";
    exit;
}


#if ($fasta_out) {
#    open (FASTAOUT,">$fasta_out") ||
#	die "Can not open fasta output file. $fasta_out\n";
#}

if ($outfile) {
    open (OUTFILE, ">$outfile") ||
	die "Can not open output file at $outfile";
}
else {
    open (OUTFILE, ">&STDOUT") ||
	die "Can not write output to STDOUT\n";
}

if ($summary_outfile) {
    open (OUTSUM, ">$summary_outfile") ||
	die "Can not open summary output file\n";
}


# SUMMARY FILE HEADER
if ($summary_outfile) {

	    print OUTSUM "id\t".
		"bin_start\t".
		"bin_end\t".
		"seq_length\t".
		# Cg count and ratio
		"count_GC\t".
		"ratio_GC\t".     # denome is just gcat
		# Cytosine methylation sites
		"count_CpG\t".
		"count_CpHpG\t".
		"count_CpHpH\t".
		# all single nucleotides
		"count_G\t".
		"count_C\t".
		"count_A\t".
		"count_T\t".
		"count_N\t".
		# dinucelotides
		"count_GG\t".
		"count_GC\t".
		"count_GA\t".
		"count_GT\t".
		"count_CG\t".
		"count_CC\t".
		"count_CA\t".
		"count_CT\t".
		"count_AG\t".
		"count_AC\t".
		"count_AA\t".
		"count_AT\t".
		"count_TG\t".
		"count_TC\t".
		"count_TA\t".
		"count_TT".
		"\n";
}

# DB Connectoin
 # Open the sequence database
my $dsn = "$database:$host";
my $db = Bio::DB::SeqFeature::Store->new( -adaptor => $adaptor,
					  -dsn     => $dsn,
					  -user    => $user,
					  -pass    =>$password,
					  -write   => 1
    );

 
# Can cycle through objects and see what you have with
#print ref($seq_object);
my $count = 0;
my $tot_num_features = 0;
my $tot_num_clusters = 0;
my $tot_num_exact_match = 0;
print STDERR "Sequence IDS:\n"
    if $verbose;  
my @ids = $db->seq_ids();





my $global_bin_count = 0;  # Allows for ordering all bins with unique integer
for my $id (@ids) {

    print STDERR "Processing $id \n";
    
    
    $count++;
    
    # Get the scaffold length
    # and skip it if length is less than bin
#    my $scaffold_length = $id->length;
    my @scaffolds = $db->features( -seqid=> $id,
				   -type =>  "scaffold");
    my $scaffold_length;
    my $num_scaffolds = 0;
    for my $scaffold (@scaffolds) {
	$num_scaffolds++;
	$scaffold_length = $scaffold->length();
    }
    if ($num_scaffolds > 1) {
	print "\a";
	print STDERR "More than one scaffold assiged to this seqid\n";
    }
    print STDERR "\tScaffold Length: ".$scaffold_length."\n";

    
    if ($scaffold_length < $bin_length) {
	print STDERR "\tScaffold legnth less than bin too small ".
	    $id."\n";
	next;
    }

    #-----------------------------+
    # FOR EACH BIN ON SCAFFOLD    |
    #-----------------------------+
    my $bin_start = 1;
    my $bin_end = $bin_start  + $bin_length - 1;


    while ($bin_start < $scaffold_length) {
	
	$global_bin_count++;
	
	# for the last bin have to use the length of the scaffold
	if ($bin_end > $scaffold_length) {
	    $bin_end = $scaffold_length;
	}

	print STDERR "\tBin: ".$bin_start."--".$bin_end."\n"
	    if $verbose;



	#-----------------------------+
	# FETCH & ANALYZE SUBSEQUENCE |
	#-----------------------------+

        #
	my $sequence = $db->fetch_sequence(-seq_id  => $id,
					   -start   => $bin_start,
					   -end     => $bin_end);


	# Make all sequence string characters uppercase
	$sequence = uc($sequence); 

	my $seq_length = length($sequence);
	print STDERR "\t\tLen:\t".$seq_length."\n" if $verbose;


	# Counts can use perl tr function
	# these will all be positive strand counts
	#-----------------------------+
	# Count nucleotides
	#-----------------------------+
	my $c_g = ($sequence =~ tr/G//);
	my $c_c = ($sequence =~ tr/C//);
	my $c_a = ($sequence =~ tr/A//);
	my $c_t = ($sequence =~ tr/T//);
	my $c_n = ($sequence =~ tr/N//);

	my $c_c_g = $c_c + $c_g;

	#-----------------------------+
	# Count dinucleotides
	#-----------------------------+
	# Could also do an insta-cast to scalar using the = () + idiom
	#my $c_gg = () = $sequence =~ /GG/g;
	# pulled code from 
	my $c_gg = scalar( @{[ $sequence =~ /GG/g ]} );
	my $c_gc = scalar( @{[ $sequence =~ /GC/g ]} );
	my $c_ga = scalar( @{[ $sequence =~ /GA/g ]} );
	my $c_gt = scalar( @{[ $sequence =~ /GT/g ]} );

	my $c_cg = scalar( @{[ $sequence =~ /CG/g ]} );
	my $c_cc = scalar( @{[ $sequence =~ /CC/g ]} );
	my $c_ca = scalar( @{[ $sequence =~ /CA/g ]} );
	my $c_ct = scalar( @{[ $sequence =~ /CT/g ]} );

	my $c_ag = scalar( @{[ $sequence =~ /AG/g ]} );
	my $c_ac = scalar( @{[ $sequence =~ /AC/g ]} );
	my $c_aa = scalar( @{[ $sequence =~ /AA/g ]} );
	my $c_at = scalar( @{[ $sequence =~ /AT/g ]} );

	my $c_tg = scalar( @{[ $sequence =~ /GG/g ]} );
	my $c_tc = scalar( @{[ $sequence =~ /GC/g ]} );
	my $c_ta = scalar( @{[ $sequence =~ /GA/g ]} );
	my $c_tt = scalar( @{[ $sequence =~ /GT/g ]} );

	#-----------------------------+
	# Count trinucleotides
	#-----------------------------+
	# Not currently doing this, but could add by
	# emulating the code above for dinucelotides
	# Other option is to do single translate of GCAT to 0123
	# and increment a 3d matrix as we roll through the sequence


	#-----------------------------+
	# Count Potential Positive 
        # StrandCytosine Methylation
	# Sites
	#-----------------------------+
	# Look for potential cytosine methylation sites
	#  CpG, CpHpG, and CpHpH sites, 
        # where H represents any nucleotide but G
	my $c_CpG = scalar( @{[ $sequence =~ /CG/g ]} );
	my $c_CpHpG = scalar( @{[ $sequence =~ /C[CAT]G/g ]} );
	my $c_CpHpH = scalar( @{[ $sequence =~ /C[CAT][CAT]/g ]} );


	#-----------------------------+
	# A COMPLETELY DIFFERENT 
	# APPROACH
	# BASED ON MATRIX COUNTING
	#-----------------------------+
#
#	# Inititalize matrix
#	my @seq_matrix;
#	for (my $i = 0; $i<5; $i++) {
#	    for (my $j = 0; $j<5; $j++) {
#		for (my $k = 0; $k<5; $k++) {
#		    $seq_matrix[$i][$j][$k] = 0;
#		}
#	    }
#	}
#
#	# Translate sequence letters to seuqennce integers
#	# this makes it quicker to increment count matrix
#	# The integer sequence is the integer translated nuc seq
#	my $int_seq = $sequence;
#	$int_seq =~ tr/GCATN/01234/;
#
#	
#	my $ss_start = 0;
#	my $ss_stop = $seq_length - 2;
#	for ($ss_start = 1; $ss_start < $ss_stop; $ss_start++) {
#	    my $ss = substr ($int_seq, $ss_start, 3);
#	    my $i = substr ($ss, 0, 1);
#	    my $j = substr ($ss, 1, 1);
#	    my $k = substr ($ss, 2, 1);
#	    
#	    
#	    $seq_matrix[$i][$j][$k] =  $seq_matrix[$i][$j][$k] + 1;
#	    
#	}
#
#	my $matrix_chg_count = $seq_matrix[1][1][0] +
#	    $seq_matrix[1][2][0] +
#	    $seq_matrix[1][3][0];
#	
#	print STDERR "\t\t\t\tMatCHG:\t".$matrix_chg_count."\n";
#	print STDERR "\t\t\t\t\t110:\t".$seq_matrix[1][1][0]."\n";
#	print STDERR "\t\t\t\t\t120:\t".$seq_matrix[1][2][0]."\n";
#	print STDERR "\t\t\t\t\t130:\t".$seq_matrix[1][3][0]."\n";

	# End of matrix approach

	my $sum_nuc = $c_g + $c_c + $c_a + $c_t;
	my $sum_all = $c_g + $c_c + $c_a + $c_t  + $c_n;

	# The ratio of gc
	my $ratio_gc;
	if ($sum_nuc == 0) {
	    $ratio_gc = 0;
	}
	else {
	    $ratio_gc = $c_c_g / $sum_nuc;
	}
	$ratio_gc = sprintf( "%.5f", $ratio_gc );


	# Summary output file matrix
	if ($summary_outfile) {
	    print OUTSUM $id."\t".
		$bin_start."\t".
		$bin_end."\t".
		$seq_length."\t".
		# Cg count and ratio
		$c_c_g."\t".
		$ratio_gc."\t".     # denome is just gcat
		# Cytosine methylation sites
		$c_CpG."\t".
		$c_CpHpG."\t".
		$c_CpHpH."\t".
		# all single nucleotides
		$c_g."\t".
		$c_c."\t".
		$c_a."\t".
		$c_t."\t".
		$c_n."\t".
		# dinucelotides
		$c_gg."\t".
		$c_gc."\t".
		$c_ga."\t".
		$c_gt."\t".
		$c_cg."\t".
		$c_cc."\t".
		$c_ca."\t".
		$c_ct."\t".
		$c_ag."\t".
		$c_ac."\t".
		$c_aa."\t".
		$c_at."\t".
		$c_tg."\t".
		$c_tc."\t".
		$c_ta."\t".
		$c_tt.
		"\n";
	}

	if ($verbose) {
	    print STDERR "\n";
	    print STDERR "\t\t\t\tCpG:\t".$c_CpG."\n";
	    print STDERR "\t\t\t\tCpHpG:\t".$c_CpHpG."\n";
	    print STDERR "\t\t\t\tCpHpH:\t".$c_CpHpH."\n";
	}

	# Setting the next bin_start and bin_end
	$bin_start = $bin_start + $bin_length;
	$bin_end = $bin_end + $bin_length;


    } # End of for each BIN segment in the scaffold

    # For each feature

    # Temp exit on couting contigs
    if ($do_test) {
	if ($count ==5) {
	    exit;
	}
    }

    my $num_qualifying_features = 0;

    


#    for my $feature (@features) {
#	
#	# Can process name here
#	my @name_parts = split(/\-/, $feature->name);
#	
#	my $family;
#	if ($name_parts[1] =~ "single") {
#	    #$family = $feature->name;
#	    $family = $parent_feature."_single";
#	}
#	elsif ( $name_parts[0] =~ m/.*R/ ) {
#	    $family = $feature->name;
#	}
#	else {
#	    $family = $name_parts[0];
#	}
#
#	# Can add logic to reject feature types here
#	# for eample by length.
#	# This will require a separate set of 
#
#	print STDERR "\t\tFeature: ".$feature."\n"
#	    if $verbose;
#
#	if ( $feature->length < $min_feature_length) {
#	    next;
#	}
#	$num_qualifying_features++;
#
#	print OUTFILE $id."\t".
#	    $parent_feature."\t".
#	    $family."\t".
#	    $feature->name."\t".
##	    $feature->start."\t".     # Feature starte (ie. TEstart)
##	    $feature->end."\t".       # Feature end (ie. TE end)
#	    $feature->length."\t".    # Length of the feature
#	    "\n";		
#
#
#    } # End of for each feature (the parent featuer in the input query) 
    

} # End of for each major sequence id


close (OUTSUM) if ($summary_outfile);
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

=over

=item Data::Range::Compare

This perl module is required to find overlaps.

=back

Other modules or software that the program is dependent on.

=head1 BUGS AND LIMITATIONS

Any known bugs and limitations will be listed here.

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
