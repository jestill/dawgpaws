#!/usr/bin/perl -w
#-----------------------------------------------------------+
#                                                           |
# db_get_te_gene_overlap.pl - Get overlap between genes & TE|
#                                                           |
#-----------------------------------------------------------+
#                                                           |
#  AUTHOR: James C. Estill                                  |
# CONTACT: JamesEstill_@_gmail.com                          |
# STARTED: 05/08/2013                                       | 
# UPDATED: 05/08/2013                                       |
#                                                           |
# DESCRIPTION:                                              |
#   Goal is to find repeatmasker hits located within any    |
#   of the potential LTR retrotransposons. Goal will be     |
#   to look for encasement of EPRV domain within LTRs.      |
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
my $parent_feature = "LTR_retrotransposon";

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

# GBROWSE
my $gb_url = "http://localhost/cgi-bin/gb2/gbrowse/amborella_te/?name=";

#-----------------------------+
# COMMAND LINE OPTIONS        |
#-----------------------------+
my $ok = GetOptions(# REQUIRED OPTIONS
		    "i|infile=s"   => \$infile,
                    "o|outfile=s"  => \$outfile,
		    "summary=s"    => \$summary_outfile,
		    "f|feature=s"  => \$parent_feature,
		    "overlap=s"    => \$overlap_feature,
		    "min-length=s" => \$min_feature_length,
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
    # User perldoc to generate the man documentation.    system ("perldoc $0");
    exit($ok ? 0 : 2);
}

if ($show_version) {
    print "\nsql_cluster_features.pl:\n".
	"Version: $VERSION\n\n";
    exit;
}


if ($fasta_out) {
    open (FASTAOUT,">$fasta_out") ||
	die "Can not open fasta output file. $fasta_out\n";
}

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

# DB Connectoin
 # Open the sequence database
my $dsn = "$database:$host";
my $db = Bio::DB::SeqFeature::Store->new( -adaptor => $adaptor,
					  -dsn     => $dsn,
					  -user    => $user,
					  -pass    =>$password,
					  -write   => 1
    );

## Are sub-features being indexed?
print "\nIndexed sub-features? : ",
  $db->index_subfeatures, "\n";

## What serializer is being used?
print "\nThe serializer is : ",
  $db->serializer, "\n";


## List all the feature types in the database:
print "\nFeature types:\n";
print "\t", join("\n\t", $db->types), "\n";


## List the feature attributes (tags) in the database:
print "\nAttributes:\n";
print "\t", join("\n\t", $db->attributes), "\n";
 
## How many sequence ids in the database:
print "\nThere are : ", scalar($db->seq_ids), " sequences in the database\n";

# Can cycle through objects and see what you have with
#print ref($seq_object);
my $count = 0;
my $tot_num_features = 0;
my $tot_num_clusters = 0;
my $tot_num_exact_match = 0;
print STDERR "Sequence IDS:\n";  
my @ids = $db->seq_ids();

print STDERR "Searching for parent feature type: ".$parent_feature."\n";
print STDERR "Searching for overlaps with".$overlap_feature."\n";

for my $id (@ids) {

    print STDERR "Processing $id \n";

    $count++;


    # Temp exit on couting contigs
    if ($do_test) {
	if ($count ==5) {
	    exit;
	}
    }

    # Fetch some features by chromosome
    # and the parent feature type passed to the program
    # it may make sense to test for features above
    # before proceeding to cycle through all scaffolds
    my @features = $db->features( -seqid=> $id,
				  -type =>  $parent_feature);
    
    my $num_features = @features;
    # The features that qualify based on length etc.
    my $num_qualifying_features = 0;

    # The following count does not include logic to ignore short repeats
    print STDERR "\tParent feature count:\t".$num_features."\n";
    
    # cycle through matching parent_features and look for overlaps
    # in the seqid of interest

    # The following used to sum overlapping features for the feature
    my $seq_overlapping_features = 0;

    for my $feature (@features) {
	

	# Can add logic to reject feature types here
	# for eample by length.
	# This will require a separate set of 

	print STDERR "\t\tFeature: ".$feature."\n"
	    if $verbose;
#	print STDERR "\t\tLength".$feature->length."\n";

	if ( $feature->length < $min_feature_length) {
#	    print STDERR "This TE is Too Small\t".
#		$feature->length." is less than".$min_feature_length."\n";
	    next;
	}
	$num_qualifying_features++;

	# It is possible to fetch the feature start and end
	print STDERR "\t\t\tPrimary-ID:".$feature->primary_id."\n"
	    if $verbose;
	print STDERR "\t\t\t\t".$feature->start."--".$feature->end."\n"
	    if $verbose;

	my $segment  = $db->segment($id, $feature->start, $feature->end);

	my @overlapping_features = $segment->features ( -type=> $overlap_feature,
							-seqid=> $id);

	my $num_overlapping_features =  @overlapping_features;
	
	if ($num_overlapping_features > 0) {

	    $seq_overlapping_features++;


	    # for EACH overlapping feature
	    # since in some situations there will be more than
	    # on overlapping feature.
	    for my $ov_feature (@overlapping_features) {
		# id is the scaffold
		# Try to produce URL outfile as
		# http://localhost/cgi-bin/gb2/gbrowse/amborella_te/?name=AmTr_v1.0_scaffold00078%3A2807578..2831520

		# Details for this individual overlap
		my $ov_url = $gb_url.$id."%3A".$ov_feature->start."..".$ov_feature->end;

		if ($ov_feature->name) {
		    print OUTFILE $id."\t".
			$feature->name."\t".
			$feature->start."\t".     # Feature starte (ie. TEstart)
			$feature->end."\t".       # Feature end (ie. TE end)
			$feature->length."\t".    # Length of the feature
			$ov_feature->name."\t". 
			$ov_feature->start."\t".  # Start of the overlapping feature (ie gene start)
			$ov_feature->end."\t".    # End of the overlapping feature (ie. gene end)
			$ov_feature->length."\t". # Overlapping feature length
			$ov_url."\t".
			"\n";		
		}
		else {
		    print OUTFILE $id."\t".
			$feature->name."\t".
			$feature->start."\t".     # Feature starte (ie. TEstart)
			$feature->end."\t".       # Feature end (ie. TE end)
			$feature->length."\t".    # Length of the feature
#			$ov_feature->name."\t". 
			"unnamed_feature\t". 
			$ov_feature->start."\t".  # Start of the overlapping feature (ie gene start)
			$ov_feature->end."\t".    # End of the overlapping feature (ie. gene end)
			$ov_feature->length."\t". # Overlapping feature length
			$ov_url."\t".
			"\n";		
		}

		# Could possible drill down overlap feature here to find that these are in introns

	    }

	} # end of if over lap feature more than 0
	

    } # End of for each feature (the parent featuer in the input query) 
    
    # Report the number of overlapping features
    print STDERR "\tQualifying feature:\t".$num_qualifying_features."\n";
    print STDERR "\tOverlapping features:\t".$seq_overlapping_features."\n";


    # Print summary output if path provided
    my $qual_info = "none";
    if ( $min_feature_length ) {
	$qual_info = "min_len=".$min_feature_length;
    }
    if ($summary_outfile) {
	print OUTSUM $id."\t".
	    $parent_feature."\t".
	    $overlap_feature."\t".
	    $qual_info."\t".
	    $num_features."\t".                 # all feature count
	    $num_qualifying_features."\t".      # qual feature count
	    $seq_overlapping_features."\t".     # overlaps
	    "\n";
    }

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
