#!/usr/bin/perl -w
#-----------------------------------------------------------+
#                                                           |
# db_get_ltrs_and_internal
#   - Get both ltr sequences and internal sequences
#                                                           |
#-----------------------------------------------------------+
#                                                           |
#  AUTHOR: James C. Estill                                  |
# CONTACT: JamesEstill_@_gmail.com                          |
# STARTED: 09/14/2012                                       |
# UPDATED: 09/14/2012                                       |
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

# Note it is possible to narrow down feature to a single source by using
# jamie-macbook-pro:ltr_db_work jestill$ db_assign_ltr_domains.pl --verbose --feature "LTR_retrotransposon:ltr_harvest_65_5000"
# jamie-macbook-pro:ltr_db_work jestill$ db_assign_ltr_domains.pl --verbose --feature "LTR_retrotransposon:LTR_Struc"

# This is using GFF file and object
# "Bio::SeqFeature::Annotated
# Methods documented at
# http://doc.bioperl.org/bioperl-live/Bio/SeqFeature/Annotated.html
# It would be better to load these to a GFF data store and use features
# to fetch the overlapping features. This will be more scalable
# http://search.cpan.org/~cjfields/BioPerl-1.6.1/Bio/DB/SeqFeature/Segment.pm
# #
#
#Http://www.bioperl.org/wiki/Module:Bio::Db::SeqFeature::Store
# Methods at
# Http://search.cpan.org/~cjfields/BioPerl-1.6.1/Bio/DB/SeqFeature/Store.pm
# Example of adding new exon
# http://doc.bioperl.org/bioperl-live/Bio/DB/SeqFeature/NormalizedFeature.html
package DAWGPAWS;

#-----------------------------+
# INCLUDES                    |
#-----------------------------+
use strict;
use warnings;
use Text::Wrap;                # Allows word wrapping and hanging indents
use Getopt::Long;
# The following needed for printing help
use Pod::Select;               # Print subsections of POD documentation
use Pod::Text;                 # Print POD doc as formatted text file
use IO::Scalar;                # For print_help subfunction
use IO::Pipe;                  # Pipe for STDIN, STDOUT for POD docs
use File::Spec;                # Convert a relative path to an abosolute path

# Use hmmer to search for domain in database
use Bio::Tools::HMMER::Results;# Bioperl HMMER results parser

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
my $hmmer_model_dir;
my $tmp_out_dir = "temp";

#-----------------------------+
# COMMAND LINE OPTIONS        |
#-----------------------------+
my $ok = GetOptions(# REQUIRED OPTIONS
		    "i|infile=s"  => \$infile,
                    "o|outfile=s" => \$outfile,
		    "f|feature=s" => \$parent_feature,
		    "fasta=s"     => \$fasta_out, 
		    "hmmer-dir=s" => \$hmmer_model_dir,
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

if ($outfile) {
    open (OUTFILE, ">$outfile") ||
	die "Can not open the output file\n";
}
else {
    open (OUTFILE, ">&STDOUT") ||
	die "Can not open the output file\n";
}
print OUTFILE "feature_id\tsource\n";


# DB Connectoin
# Open the sequence database
my $dsn = "$database:$host";
my $db = Bio::DB::SeqFeature::Store->new( -adaptor => $adaptor,
					  -dsn     => $dsn,
					  -user    => $user,
					  -pass    =>$password,
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
for my $id (@ids) {

    $count++;

    # Fetch some featurs by chromosome
    # and the parent feature type passed to the program
    # it may make sense to test for features above
    # before proceeding to cycle through all scaffolds
    my @features = $db->features( -seqid=> $id,
				  -type =>  $parent_feature);

    my $num_features = @features;

    for my $feature (@features) {
	print STDERR "\t\tFeature: ".$feature."\n"
	    if $verbose;
	# It is possible to fetch the feature start and end
	print STDERR "\t\t\tPrimary-ID:".$feature->primary_id."\n"
	    if $verbose;
	print STDERR "\t\t\t\t".$feature->start."--".$feature->end."\n"
	    if $verbose;


	# Get all child features
	my @children = $db->fetch_SeqFeatures($feature);
	
	my $line_out = $feature->primary_id;
	
	

	print OUTFILE $line_out."\t".$feature->source()."\n";

    } # End of for each ltr retro on the scaffold 
    
    
    if ($do_test) {
	if ($count == 1) {
	    exit;
	}
    }
    
}

if ($fasta_out) {
    close (FASTAOUT);
}

close (OUTFILE);

exit 0;

#-----------------------------------------------------------+ 
# SUBFUNCTIONS                                              |
#-----------------------------------------------------------+

sub run_hmmer {
# Run the HMMER program for a set of models

    #-----------------------------+
    # VARS PASSED TO THE SUBFUN   |
    #-----------------------------+
    my ($seq_name, $seq_path, $model_dir, $hmm_suffix,
	$out_dir, $par_name) = @_;

    # $seq_name   # name of the seq record id HEX0014K09
    # $seq_path   # Path to the seq 
    # $model_dir  # Path to the dir containing the seq models 
    # $hmm_suffix # Options to add to the hmmsearch command line
    #             # May need to chomp the suffix
    # $out_dir    # Out dir for the location of the hmm out files
    # $par_name   # Name given to the hmm parameter set

    # Append slash to end of model_dir if needed
    unless ($model_dir =~ /\/$/ ) {
	$model_dir = $model_dir."/";
    }

    unless ($out_dir =~ /\/$/ ) {
	$out_dir = $out_dir."/";
    }

    #-----------------------------+
    # VAR SCOPE AND INITIALIZE    |
    #-----------------------------+
    my $sub_proc_num = 0;           # Process number starts at zero

    my $hmm_cmd;                   # hmmsearch command line 
#    my $out_path = $out_dir.$seq_name."_hmm_".$par_name.".hmmout";
    
    # Will need to delete an already existing outdir	

    #-----------------------------+
    # OPEN MODEL DIR              |
    #-----------------------------+
    opendir( MODDIR, $model_dir ) ||
	die "Can not open the hmm model directory:\n $model_dir\n";
    # The following does not search for the hmm extension
    # and assumes that all files in this directory are hmm models

    # The following will only use files with the hmm extension
    my @models = grep /hmm?$/, readdir MODDIR;

    closedir( MODDIR );    

    my $num_models = @models;   # Number of models
    
    for my $ind_model (@models) {
	
	$sub_proc_num++; # Increment the process number

	# Model root name will be used in the naming of files
	# default is just the name of the individual model file
	# will try to trim off the hmm when it is present
	my $model_root_name = $ind_model;

	if ($ind_model =~ m/(.*)\.hmm$/) {
	    $model_root_name = $1;
	}

 	# HMMER MODEL PATH
	my $model_path = $model_dir.$ind_model;
	my $out_path = $out_dir.$seq_name."_hmm_".$par_name."_".
	    $model_root_name.".hmmout";

	#-----------------------------+
	# DOES THE HMM MODEL EXIST    |
	#-----------------------------+
	if ( (-e $model_path) && (-e $seq_path) ) {
	    

	    print "\nMOD PATH: $model_path\n" if $verbose;
	    print "SEQ PATH: $seq_path\n" if $verbose;
	    print "OUT PATH: $out_path\n" if $verbose;

#	    $hmm_cmd = "hmmsearch " . 
#		"--domT 2 $model_path $seq_path >$out_path";
	    
	    # Alternative may be to set a single $out_path
	    # However, the hmmer parser may have a problem with this
	    $hmm_cmd = "hmmsearch " . 
		"--domT 2 $model_path $seq_path >>$out_path";

	    #-----------------------------+
	    # PRINT STATUS OF EACH HMM    |
	    # PROCESS                     |
	    #-----------------------------+	
	    if ($verbose) {
		print "HMM model ".$sub_proc_num." of ".$num_models."\n";
		print wrap("\t", "\t", $hmm_cmd );
		print "\n";
	    }
	 
	    # Run the hmm command if this is not a test
	    # commented out, test is used above in a different context
#	    system ($hmm_cmd) if (!$do_test);
	    system ($hmm_cmd);

	}
	else {
	    
	    # Show error if the model_path or input sequence
	    # could not be located
	   
	    unless (-e $model_path) {
		print STDERR "ERROR: Can not find model:\n$ind_model\n";
	    }
	    
	    unless (-e $seq_path) {
		print STDERR "ERROR: Can not find sequense:\n$seq_path\n";
	    }

	} # End of test for existence of model_path and seq_path
       
	
    } # End of for each database loop


} # End of run hmmer subfunction


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
