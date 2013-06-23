#!/usr/bin/perl -W
#-----------------------------------------------------------+
#                                                           |
# db_add_intro_order.pl - Add intron order                  |
#                                                           |
#-----------------------------------------------------------+
#                                                           |
#  AUTHOR: James C. Estill                                  |
# CONTACT: JamesEstill_@_gmail.com                          |
# STARTED: 06/14/2013                                       |
# UPDATED: 06/14/2013                                       |
#                                                           |
# DESCRIPTION:                                              |
#  Similar code can be used to push other values to the     |
#  attributes of features already in the database. For      |
#  example to add GO ontology reference IDs to the database |
#  or gene family information.                              |
#  This would need to be given a list of GO values for      |
#  each gene using a gene_id or name that could be looked   |
#  up in the database.                                      |
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
my $parent_feature = "gene";
my $overlap_feature = "test";
my $child_feature = "intron";  # The child feature of interest that will
                               # be VARIABLE

#-----------------------------+
# modified SCOPE              |
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
		    "child=s"      => \$child_feature,
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
					  -pass    => $password,
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


    print STDERR "Working on each feature\n";
    for my $feature (@features) {
	
	print STDERR "\t\tFeature: ".$feature."\n"
	    if $verbose;

	$num_qualifying_features++;

	# It is possible to fetch the feature start and end
	print STDERR "\t\t\tPrimary-ID:".$feature->primary_id."\n"
	    if $verbose;
	print STDERR "\t\t\t\t".$feature->start."--".$feature->end."\n"
	    if $verbose;


	# Get the child features

	my @child_features = $db->fetch_SeqFeatures($feature);
	my $num_child_features =  @child_features;
	
	if ($num_child_features > 0) {
	    # for EACH overlapping feature
	    # since in some situations there will be more than
	    # on overlapping feature.
	    for my $child (@child_features) {
		# Details for this individual overlap
#		print STDERR "\t\t\t".$id."\t".
#		    $feature->name."\t".
#		    $feature->start."\t".     # Feature starte (ie. TEstart)
#		    $feature->end."\t".       # Feature end (ie. TE end)
#		    $feature->length."\t".    # Length of the feature
#		    $child->type."\t".
#		    $child->start."\t".  # Start of the overlapping feature (ie gene start)
#		    $child->end."\t".    # End of the overlapping feature (ie. gene end)
#		    $child->length."\t". # Overlapping feature length
#		    "\n";		

		
		# The child of gene is mRNA so need to go one more level down for intron
		my @grandchild_features = $db->fetch_SeqFeatures($child);
		my $num_grandchild_features = @grandchild_features;



		# Temp array to hold information on introns
		my @intron_features;
		my @introns_sorted;
		my $intron_num = 0;
		
		for my $grandchild (@grandchild_features) {


		    # Need to load intron to array
		    # then sort start (up or down based on strand)
		    # and add intron order.
		    print STDERR "\t\t\t\t".$grandchild->type."\n"
			if $verbose;
		    if ($grandchild->type =~ $child_feature) {

			print STDERR $feature."\t".$feature->strand."\t".
			    "ID".$grandchild->primary_id."\t".
			    $grandchild->start."--".$grandchild->end."\n";

			# This will load starting at the 0 position
			$intron_features[$intron_num][0] = $grandchild->primary_id;
			$intron_features[$intron_num][1] = $grandchild->start;
			$intron_num++;
			    
		    }



		}

		#////////// WORKING ON SORTING HERE
		# Now we have all of the introns in the intron_features array
		if ($intron_num > 0) {
		    if ($feature->strand == 1) {
			print STDERR "Positive strand features\n";
			# sort by ascending order on start column
			@introns_sorted = sort { $a->[1] <=> $b->[1] } @intron_features;
		    }
		    elsif ($feature->strand == -1) {
			print STDERR "Negative strand features\n";
			# sort by descending order on start column
			@introns_sorted = sort { $b->[1] <=> $a->[1] } @intron_features;
			
		    }
		} # End of if we have introns
		
		print STDERR "Sorted Introns:\n"
		    if $intron_num > 0;

		# Need to make sure this stops at the correct place
		my $intron_position = 0;
		for (my $i = 0; $i<$intron_num; $i++) {
		    $intron_position++;
		    print STDERR $introns_sorted[$i][0]."\t".
			$introns_sorted[$i][1]."\t".
			$intron_position."\n";
		    
		    # the current intron feature as a database object
		    # fetched using the 
		    my $intron_feature=$db->fetch($introns_sorted[$i][0]);


		    # Trying treating this as a direct hash
		    # Using this seems to work since it does fetch the size
		    my %intron_feature_attributes = $intron_feature->attributes();
		    print "size of hash:  " . keys( %intron_feature_attributes ) . ".\n";
		    # Temp show the key and values in the hash
		    for my $key ( keys %intron_feature_attributes ) {
			my $value = $intron_feature_attributes{$key};
			print STDERR "$key => $value\n";
		    }

		    # If this does not work can print an output file that gives
		    # the key values to update the database then make a text.
		    
		    # file to manually update the intron ids numbes in the database.
#  -attributes   a hashref of tag value attributes, in which the key is the tag
#                  and the value is an array reference of values
# Need to possbily push the value into an array reference, and then set it
# to be matched to the desired key.
# array reference information at 
# http://www.thegeekstuff.com/2010/06/perl-array-reference-examples/
		    # Temp array to hold the intron position to push to the
		    # has that will be loaded to the BioDB
		    my @intron_pos_array;
		    $intron_pos_array[0] = $intron_position;
		    my $intron_pos_array_ref = \@intron_pos_array;
		    # Next need to push the array reference to the 
		    # %intron_feature_attributes hash using a value
		    # already in the database which is intron_order ...
		    # other tags that could be used accodring to MISO are
		    #  - interior_intron
                    #  - UTR_intron
                    #  - five_prime_intron
                    #  - three_prime_intron
		    # should push as:
		    # 
		    $intron_feature_attributes{intron_order} = $intron_pos_array_ref;

		    # Then test that the array reference was pushed
		    print STDERR "Newly updated!\n";
		    for my $key ( keys %intron_feature_attributes ) {
			my $value = $intron_feature_attributes{$key};
			print STDERR "$key => $value\n";
			if ($key =~ "intron_order") {
			    my @returned_intron_order_array = @$value;
			    my $returned_first_value =  $returned_intron_order_array[0];
			    print STDERR "The intron order value stored is:".
				$returned_first_value."\n";
			}
		    }
		    print STDERR "\n";


		    

# Then need to update the feature as
# 
# Title   : update
# Usage   : $flag = $feature->update()
# Function: Update feature in the database
# Returns : true if successful

		    # If everything checks out
		    # then need to updated the feature to include the new tag value pair
		    # hopefully 
		    # Do not update the feature if we are in test mode.
#		    unless ($do_test) {
		    my $intron_feature_hash_ref = \%intron_feature_attributes;
		    # trying this with a has reference
#		    $intron_feature->attributes( $intron_feature_hash_ref );
		    $intron_feature->attributes( "$intron_feature_hash_ref" );
# The following did not work		    
#		    $intron_feature->attributes() = $intron_feature_hash_ref;
# The following also does not work
#		    my $flag = $intron_feature->update();
#		    $db->update($intron_feature);
#		    my $flag = $db->update($intron_feature);
		    # Should be able to update as shown at
		    # http://search.cpan.org/~cjfields/BioPerl-1.6.1/Bio/DB/SeqFeature/Store/DBI/mysql.pm
#		    my $flag = $db->update($intron_feature);
# The following throws an error ...
#		    $db->update($intron_feature) || die "Could not update";
#		    unless ($flag) {
#			print STDERR "This record did not update\n";
#		    }
#		    else {
#			print STDERR "Update successfull\n" if $verbose;
#		    }
#			if ($flag)
#		    }

		    #////////////////////////// AT WORK HERE

		    
		    
#		    # print the existing attributes
#		    for (keys %$intron_feature_attributes) {
#			delete $href->{$_};
#		    }

		}
		print STDERR "\n" if $intron_num > 0;

		# Will need to add the attribute as needed similar to below
		  # change the feature and update it
#		$f->start(100);
#		$db->update($f) or die "Couldn't update!";
		# THe existing features would be returned by a hash ... so
		# may want to first establish an empty hash
		# load existing features, push new features
		# and then to a save
		#@features = $db->get_features_by_attribute({description => 'protein kinase'})
		
	    }



	    

	} # end of if over lap feature more than 0	

    } # End of for each feature (the parent featuer in the input query) 
    
    # Report the number of overlapping features


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
