#!/usr/bin/perl -w
#-----------------------------------------------------------+
#                                                           |
# db_get_bin_details.pl - Get details for each spatial bin  |
#                                                           |
#-----------------------------------------------------------+
#                                                           |
#  AUTHOR: James C. Estill                                  |
# CONTACT: JamesEstill_@_gmail.com                          |
# STARTED: 05/13/2013                                       |
# UPDATED: 05/13/2013                                       |
#                                                           |
# DESCRIPTION:                                              |
#  Given a bin size generate summary information for each
#  bin for the feature. 
#  Composition AND count for most features such as:
#     - gaps
#     - 
#   
#    ie 
#    gene .. gene length, 
#    
#  this will produce a table that can be used to make a 
#  figure in another program such as circos to summarize
#  the data. Will include logic to check that the scaffolds
#  are at least the length of the bin, and will not process |
#  scaffolds that are less than the bin lenght.             |
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

# the list of features to search for in each bin
my @search_features = ( "gap",    # Scaffold gaps
			"gene:EVM_27",   # Number of genes
			"exon:EVM_27",   # Exons 
			"intron:EVM_WI",   # introns
			"MITE",   # MITES
			"DHH",    # Helitron
			"DTA",    # Hat
			"DTM",    # Mutator
			"DTH",    # PIF-Harbinger
			"RLEPRV", # EPRVs
			"RLG",    # Gypsy LTR Retrotransposons
			"RLC",    # Copia LTR Retrotransposons
			"RIL",    # L1 LINEs
			"RIT",    # RTE LINEs
    );

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


my $num_search_features = @search_features;

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
    print OUTSUM "unique_bin_id\tscaffold_id\tbin_start\tbin_end\t";
    my $feature_header_count = 0;
    for my $feature_header (@search_features) {
	$feature_header_count++;
	
	if ($feature_header_count == $num_search_features) {
	    print OUTSUM $feature_header."_count\t";
	    print OUTSUM $feature_header."_composition\t";
	    print OUTSUM $feature_header."_ratio\n";
	}
	else {
	    print OUTSUM $feature_header."_count\t";
	    print OUTSUM $feature_header."_composition\t";
	    print OUTSUM $feature_header."_ratio\t";
	}
    
    }
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




# For testing
#$parent_feature = $search_features[0];
my $search_feature = $search_features[0];






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

	my $feature_num = 0;   # To keep location in search_features list
	for my $search_feature (@search_features) {
	    
	    $feature_num++;
	    
	    #-----------------------------+
	    # FOR EACH FEATURE IN LIST    |
	    #-----------------------------+
	    # Will also need to add logic below to trim features
	    # that range partially outside of the bin
	    
	    # Fetch some features by chromosome
	    # and the parent feature type passed to the program
	    # it may make sense to test for features above
	    # before proceeding to cycle through all scaffolds
	    my @features = $db->features( -seqid=> $id,
					  -start=> $bin_start,
					  -end=>   $bin_end,
					  -type => $search_feature);
	    
	    my $num_features = @features;
	    
	    # For each individual feature in the bin
	    # we need to assess the $start and end
	    # to see if fully within the area of interest
	    # and trim if not when calculating composition
	    my $feature_composition = 0;  # Compositiion (length
	    for my $ind_feature (@features) {
		my $feature_start = $ind_feature->start;
		my $feature_end = $ind_feature->end;

		# Trim end if necessary
		if ($feature_end > $bin_end) {
		    $feature_end = $bin_end;
#		    print STDERR "Trimming end\n";
		}
		
		if ($feature_start < $bin_start) {
		    $feature_start = $bin_start;
#		    print STDERR "Trimming start";
		}

		my $feature_length = $feature_end - $feature_start + 1;
		
		$feature_composition = $feature_composition + $feature_length;
	    }


	    my $bin_length = $bin_end - $bin_start + 1;
	    my $feature_ratio = $feature_composition / $bin_length;

	    # The following does not include logic to ignore
	    # features below a threshold
	    print STDERR "\t\t ".$search_feature.
		" feature_num:\t".$num_features."\t".
		" feature_comp:\t".$feature_composition.
		" feature_ratio:\t".$feature_ratio.
		"\n"
		if $verbose;

	    # Can also add a diversity metric here!
	    
	    # Print information to output file
	    # This will be done to allow for cross-tab type
	    # queries in a database or analyis package such as R
	    # Name the bin as 
#	    print OUTFILE $global_bin_count."\t".
#		$id.":".$bin_start."..".$bin_end."\t".
#		$search_feature."\t".
#		"feature_num\t".$num_features.
#		"\n";

	    # If we have a summary outfile
	    # could possibly have all information... in columns
	    if ($summary_outfile) {

		if ($feature_num == 1) {
		    # First feature in list 
		    print OUTSUM $global_bin_count."\t";
#		    print OUTSUM $id.":".$bin_start."..".$bin_end."\t";
		    print OUTSUM $id."\t";
		    print OUTSUM $bin_start."\t";
		    print OUTSUM $bin_end."\t";
		    print OUTSUM $num_features."\t".
			$feature_composition."\t".
			$feature_ratio."\t";
		}
		elsif ($feature_num == $num_search_features) {
		    # Last feature in list include line ending
		    print OUTSUM $num_features."\t".
			$feature_composition."\t".
			$feature_ratio."\n";
		}
		else {
		    print OUTSUM $num_features."\t".
			$feature_composition."\t".
			$feature_ratio."\t";
		    
		}

	    }
	    

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
