#!/usr/bin/perl -w
#-----------------------------------------------------------+
#                                                           |
# sql_cluster_features.pl - Spatially cluster features      |
#                                                           |
#-----------------------------------------------------------+
#                                                           |
#  AUTHOR: James C. Estill                                  |
# CONTACT: JamesEstill_@_gmail.com                          |
# STARTED: 08/17/2012                                       |
# UPDATED: 08/17/2012                                       |
#                                                           |
# DESCRIPTION:                                              |
#  Short Program Description                                |
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

# DB variables
my $adaptor = 'DBI::mysql';
my $user = 'dawgpaws';
my $password ='f34tur3s';
my $host = 'localhost';
my $database = 'amtr_ltr';
my $fasta_out;
#-----------------------------+
# COMMAND LINE OPTIONS        |
#-----------------------------+
my $ok = GetOptions(# REQUIRED OPTIONS
		    "i|infile=s"  => \$infile,
                    "o|outfile=s" => \$outfile,
		    "f|feature=s" => \$parent_feature,
		    "fasta=s"     => \$fasta_out, 
		    # Seq data store database options
		    "p|pass=s"    => \$password,
		    "u|user=s"    => \$user,
		    "host"        => \$host,
		    "adaptor"     => \$adaptor,
		    "database"    => \$database,
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
    print "\nsql_cluster_features.pl:\n".
	"Version: $VERSION\n\n";
    exit;
}


if ($fasta_out) {
    open (FASTAOUT,">$fasta_out") ||
	die "Can not open fasta output file. $fasta_out\n";
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
for my $id (@ids) {

    $count++;

    # Fetch some featurs by chromosome
    # and the parent feature type passed to the program
    # it may make sense to test for features above
    # before proceeding to cycle through all scaffolds
    my @features = $db->features( -seqid=> $id,
				  -type =>  $parent_feature);

    # ESTABLISH THE GRAPH OBJECT
    my $feature_graph = Graph->new (directed => 0);

    # LOAD PRIMARY IDS OF THE FEATURES AS NODES TO THE GRAPH
    for my $feature (@features) {
	$feature_graph->add_vertex( $feature->primary_id );
    }

    my $num_features = @features;

    # cycle through matching parent_features and look for overlaps
    # in the seqid of interest
    for my $feature (@features) {
	print STDERR "\t\tFeature: ".$feature."\n"
	    if $verbose;
	# It is possible to fetch the feature start and end
	print STDERR "\t\t\tPrimary-ID:".$feature->primary_id."\n"
	    if $verbose;
	print STDERR "\t\t\t\t".$feature->start."--".$feature->end."\n"
	    if $verbose;
#	my $feature_segment = $db->segment
	# Get feature location
	# use location to get segments
	my $segment  = $db->segment($id, $feature->start, $feature->end);
	
	# Get the features that match this segment
	# Using type from above to select the same features type
	# however this would be the ploce to load a sepearate feature name
	my @overlapping_features = $segment->features ( -type=> $parent_feature,
							-seqid=> $id);
	# count the number of overlapping features and report to STDERR
	my $num_overlapping_features =  @overlapping_features;
	print STDERR "\t\t\t\t Overlaping features:".$num_overlapping_features."\n"
	    if $verbose;

	# This is where the overlapping features edges are loaded into the network
	# THis is using the primary_id from the database, so the db
	# will be needed later to fetch information for the graph
	for my $overlapping_feature (@overlapping_features) {
	    print STDERR "\t\t\t\t" if $verbose;
	    print STDERR $feature->primary_id."-->".
		$overlapping_feature->primary_id."\n" if $verbose;
	    $feature_graph->add_edge($feature->primary_id,
				     $overlapping_feature->primary_id);





	    # FIND EXACT MATCHES WITH DIFFERENT IDS FROM THE DATABASE
	    # Test for exact overlap from a different source
	    # this may just be the caes of the same feautre
	    if ( $feature->equals($overlapping_feature) ) {
		if ($feature->primary_id ne $overlapping_feature->primary_id) {
		    print STDERR "\t\t\tEXACT MATCH\n" if $verbose;
		    $tot_num_exact_match++;
		}
	    }




	}

	# Find overlapping features and load to array
	# This does not appear to work and throws error
	# MSG: Can't locate object method "features" via package "Bio::DB::SeqFeature"
#	my @overlapping_features = $feature->features;
	
#	my @overlaps = $feature->features( -type => $parent_feature,
#					   -seqid=> $id);
    } # End of looking for overlaping features

    # Once we have graph we can cluster the graph here

    # once we have clusters we can cycle through the clusters
    # and add spatial_cluster_id to the database if desired
    # this is also where we will cycle through all nodes in the cluster
    # and get the longest LTR to be representative of that cluster
    # if the LTRs are in alignment. otherwise choose mulitple LTRs
    
    
    # This is where we find overlaps for each feature
    # and load to graph for connected compoennts spatial clustering
    # could also use this to generate sequences for ABA blast
    # this would use the fetch_sequence function.

    if ($num_features > 0 ) {
	print STDERR "\t".$id."\n";
	print STDERR "\t\tFeatures of type ".$parent_feature." : ".$num_features."\n";
	$tot_num_features = $tot_num_features + $num_features;

	my @cc = $feature_graph->connected_components();
	
	my $num_cc = @cc;
	$tot_num_clusters = $tot_num_clusters + $num_cc;
	print STDERR "\t\tTotal Clusters: ".$num_cc."\n";


	# Do something with each cluster something like
	# Will name as seq_id:feature_name:cluster_num
	my $clust_num = 0;
	foreach my $comp_list (@cc) {
	    $clust_num++;
	    my $num_clust_nodes = @$comp_list;
	    print STDERR "\t\t\tCluster_".$clust_num.":".$num_clust_nodes."\n";
	    
	    my $min_start;
	    my $max_end;
	    # For each compent

	    # Diagnostic LTR info
	    my $max_len_ltr = 0;
	    my $max_len_ltr_seq;
	    my $max_len_ltr_id;

	    # The $ind_comp is the ID of hte LTR_retrotransposon in the database
	    for my $ind_comp (@$comp_list) {
		print STDERR "\t\t\t\t".$ind_comp."\n";
#		    $ind_comp->start."-".$ind_comp->end."\n";
		# Get the feature using the primary id
		my $ind_feature = $db->get_feature_by_primary_id($ind_comp);


		

		print STDERR "\t\t\t\t".$ind_feature->start."-".$ind_feature->end."\n";
		if ($min_start) {
		    if ($ind_feature->start < $min_start) {
			$min_start = $ind_feature->start;
		    };
		}
		else {
		    $min_start = $ind_feature->start;
		}

		if ($max_end) {
		    if ($ind_feature->end > $max_end) {
			$max_end = $ind_feature->end;
		    }
		}
		else {
		    $max_end = $ind_feature->end;
		}

		# Fetch the children of the parent feature
		my @children = $db->fetch_SeqFeatures($ind_feature);
		
		
		# This prints all the children features
		# here we are looking for either
		# 1. Long_terminal_repeat
		# 2. five_prime_LTR
		# 3. three_prime_LTR
		# and then get the longest LTR of the cluster
		
		#// MAY NEED TO INITIALIZE THIS ABOVE


		for my $child (@children) {
#		    print $child->
		    #print STDERR "\t\t\t\t\t".$child."\t".$child->start."--".
		    #$child->end.":".$len."\n"
		    # If child is an LTR
		    if ($child =~ "long_terminal_repeat" || 
			$child =~ "five_prime_LTR" ||
			$child =~ "three_prime_LTR" ) {
			my $len = $child->end - $child->start + 1;
			$child->end.":".$len."\n";
			if ($len > $max_len_ltr) {
			    $max_len_ltr = $len;
			    $max_len_ltr_seq = $db->fetch_sequence($id, $child->start, $child->end);
			    $max_len_ltr_id = $ind_comp;
			}
		    }
		}


	    }



	    # THE FOLLOWING PRINTS THE REPRESENTATIVE SEQUENCE FOR THE
	    # CLUSTER ID
	    my $cluster_id = $id.":".$min_start."..".$max_end;
	    print STDERR "\t\t\tRepresentative:\n";
	    my $url = "http://localhost/cgi-bin/gb2/gbrowse/amborella_te/?name=".
		$cluster_id;

#	    print STDERR "\t\t\t\t".$max_len_ltr_id."|".$id.":".$min_start."..".$max_end.
	    print STDERR "\t\t\t\t".$max_len_ltr_id."|".$cluster_id.
		"\t".$max_len_ltr_seq."\n";
	    print STDERR "\t\t\tCluster Span:".$min_start."::".$max_end."\n\n";
	    print STDERR "\t\t\tCluster URL: ".$url."\n";
	    # The can format GBrowse URL as

	    if ($fasta_out) {
		print FASTAOUT ">".$max_len_ltr_id."|".$cluster_id."\n";
		print FASTAOUT $max_len_ltr_seq."\n";
	    }
	    
	    # Then to add cluster ID to each would be to cycle back through
	    # and add using the id from ind_comp in the graph
	    # This would need to add an attribute to an existing features
	    # 1. make sure item exists in attributelist and get attribute_id
	    # 2. load attribute table as
	    #     id=id of node
            #     attribute_id=id from above(ltr_clust_id)
	    #     attribute_value=$cluster_id
	    # could do this as follows of as separate
	    # DBI connection to the database with direct SQL
	    for my $ind_comp (@$comp_list) {
             	my $ind_feature = $db->get_feature_by_primary_id($ind_comp);
#		my ($attributes) = $ind_feature->attributes();
#		my $ltr_similarity = $ind_feature
		if (my $t = ($ind_feature->attributes('LTR_similarity'))[0]) {
		    print STDOUT "HAS LTR SIMILARITY\t";
		    print STDOUT "\t".$t."\n";
		}

		# The following sets the spatial cluster values of the object
		# take FROM the database but does not update the database itself
		if (my $t = ($ind_feature->attributes('TE_spatial_cluster'))[0]) {
		    print STDERR "Has existing spatial cluster. Will not update\n";
		}
		else {
		    my @spatial_cluster_values;
		    $spatial_cluster_values[0] = $cluster_id;
		    $ind_feature->{attributes}{'TE_spatial_cluster'} = \@spatial_cluster_values;
		    # Only update when spatial cluster values are blank
#		    print STDERR "\t\tSetting TE spatial cluster to: ".$set_cluster."\n";
		    $ind_feature->update ||
			die "Can not update this feature\n";
		    
		}

#		# Test that it was set to the attribute outside of the database
#		if (my $set_cluster = ($ind_feature->attributes('TE_spatial_cluster'))[0]) {
#		}

		# To update the database would require

		# could possibly add new information as
		#$ind_comp->{attributes}{TE_spatial_cluster}[0] = $cluster_id;

#		print STDERR $attributes->('LTR)similarity');
#		my $attributes = $ind_feature->attributes();
#		my ($targetdesc) = ($ind_comp->attributes('LTR_similarity'))[0];
#		my ($targetdesc) = $ind_comp->attributes('LTR_similarity');
#		my ($ltr_similarity) = ($ind_comp>attributes('LTR_similarity'))[0];

#		my @attributes = $ind_feature->attributes();
#		for my $ind_attribute (@attributes) {
#		    print STDERR "\t".$ind_attribute."\n";
#		}
#		print STDERR "Attributes:".$attributes."\n";
		
	    }

# It looks like it is possible to change the feature object and then update it
# as # get the feature back out:
# http://search.cpan.org/~cjfields/BioPerl-1.6.901/Bio/DB/SeqFeature/Store/DBI/mysql.pm
#  my $f  = $db->fetch($id);
#  # change the feature and update it
#  $f->start(100);
#  $f->update($f) or die "Couldn't update!";
# In this case we already have the feature as
# $ind_feature
# so the goal is to add to the attributes has of $ind_feature
#
	    # Send Full sequence span on to dotter if desired
	}
	
	
    }

    # exit while testing
    if ($do_test) {
	if ($count == 5) {
	    exit;
	}
    }

}
# It is possible to fetch sequence as
# $sequence = $db->fetch_sequence(-seq_id=>$seqid,-start=>$start,-end=>$end)


print STDERR "\nTotal matching features: ".$tot_num_features."\n";
print STDERR "Total clusters: ".$tot_num_clusters."\n";
print STDERR "Total exact matches: ".$tot_num_exact_match."\n";

## The follow working code gets the features
## and counts the numbe of features of this type
#my $count = 0;
#for( $db->types ){
#  print "Feature $_\n";
#  my @feats =
#    $db->get_features_by_type( $_ );
#  print "got ", scalar(@feats), " $_ features\n";
#}
if ($fasta_out) {
    close (FASTAOUT);
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
