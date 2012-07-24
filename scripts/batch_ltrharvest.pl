#!/usr/bin/perl -w
#-----------------------------------------------------------+
#                                                           |
# batch_ltrharvest.pl - Run LTRHarvest/Digest in batch mode |
#                                                           |
#-----------------------------------------------------------+
#                                                           |
#  AUTHOR: James C. Estill                                  |
# CONTACT: JamesEstill_@_gmail.com                          |
# STARTED: 07/20/2012                                       |
# UPDATED: 07/20/2012                                       |
#                                                           |
# DESCRIPTION:                                              |
#  Run the LTRHarvest program in batch mode.  This will     |
#  the following basic steps.                               |
#   1) Index the sequence                                   |
#   2) Run LTRHarvest on the Index                          |
#   3) Run LTRDigest on the LTRHarvest GFF Results          |
#   4) Convert GFF results from Genome Tools format to      |
#      standard sequence ontology complient GFF.            |
#   5) Move or remove the intermediate files.               |
#  Genome tools can be a bit buggy, so need to check that   |
#  all expected files and report errors if needed. Then     |
#  skip on to next scaffold as needed.                      |
#                                                           |
# USAGE:                                                    |
#  batch_ltrharvest.pl -i InFasta.fasta -o annotations/     |
#                      -m pHMM_Models/                      |
#                      -t tRNAFile.fasta                    |
#                                                           |
# VERSION: $Rev$                                            |
#                                                           |
# LICENSE:                                                  |
#  GNU General Public License, Version 3                    |
#  http://www.gnu.org/licenses/gpl.html                     |  
#                                                           |
#-----------------------------------------------------------+
#
# Usage as: 
#  batch_ltrharvest.pl -i test_in/ -o test_out/ --hmm hmm/ --verbose
#  batch_ltrharvest.pl -i test_in/ -o test_out/ --hmm hmm/ -t tRNA_aragorn_default.fasta --verbose
#
# Run index on fasta file as
#  gt suffixerator -Db AmTr_v1.0_scaffold00001.fasta -indexname test2.idx -suf -lcp -des -ssp -sds -dna
# Run LTRHarvest on fasta file as:
#  gt ltrharvest -index test2.idx -gff3 test2.out.gff3
# Run LTRDigest on fasta file as:
#  gt ltrdigest -hmms PF10551.hmm PF03108.hmm PF00075.hmm PF00552.hmm PF00077.hmm PF02093.hmm PF02813.hmm PF10536.hmm PF00078.hmm PF00096.hmm PF00872.hmm PF04195.hmm PF06815.hmm PF06817.hmm PF02992.hmm PF07727.hmm PF02994.hmm PF02022.hmm PF01710.hmm PF03017.hmm PF00607.hmm PF03732.hmm PF00385.hmm PF00665.hmm PF01359.hmm PF04827.hmm PF01498.hmm PF08284.hmm PF05699.hmm PF09299.hmm PF03101.hmm PF03004.hmm -- scaffold00001.ltrharvest.gff3 AmTr_v1.0_scaffold00001.idx
# With tRNAS as
#  gt ltrdigest -trnas Athal-tRNAs.fa -hmms PF10551.hmm PF03108.hmm PF00075.hmm PF00552.hmm PF00077.hmm PF02093.hmm PF02813.hmm PF10536.hmm PF00078.hm PF00096.hmm PF00872.hmm PF04195.hmm PF06815.hmm PF06817.hmm PF02992.hmm PF07727.hmm PF02994.hmm PF02022.hmm PF01710.hmm PF03017.hmm PF00607.hmm PF03732.hmm PF00385.hmm PF00665.hmm PF01359.hmm PF04827.hmm PF01498.hm PF08284.hmm PF05699.hmm PF09299.hmm PF03101.hmm PF03004.hmm -- scaffold00001.ltrharvest.gff3 AmTr_v1.0_scaffold00001.idx > scaffold0001_ltrs.gff3
#
# Macbook TEST AS
#  jamie-macbook-pro:test_ltrharvest jestill$ gt suffixerator -db sc1.fasta -indexname sc1.fasta -tis -suf -lcp -des -ssp -sds -dna
#  jamie-macbook-pro:test_ltrharvest jestill$ gt ltrharvest -index sc1.fasta -gff3 sc1.ltrharvest.gff3
#  jamie-macbook-pro:test_ltrharvest jestill$ gt ltrdigest sc1.ltrharvest.gff3 sc1.fasta > sc1.ltrharvest.ltrdigest.gff3
# With tRNAs as
#  jamie-macbook-pro:test_ltrharvest jestill$ gt ltrdigest -trnas tRNA_aragorn_default.fasta sc1.ltrharvest.gff3 sc1.fasta > sc1.ltrharvest.ltrdigest.gff3
# With HMMs as
#  Not available on the older version of gt I was able to get to work on my MacBook Pro
# HMMS can accept wildcards

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
# USER ENV OPTIONS            |
#-----------------------------+
# If variable not defined in the ENV then
# assume it is in the user's path
my $gt_path = $ENV{GENOME_TOOLS} ||
    "gt";

#-----------------------------+
# VARIABLE SCOPE              |
#-----------------------------+
my $indir;
my $outdir;

# Vars needing 
my $name_root;
my $gff_dir;
my $ls_out;
my $fl_gff_outpath;
my $do_gff_convert = 1;

# LTRHarvest varss
my $min_tsd = 3;

# LTRDigest vars
my $hmm_dir;                      # Dir with hmm models for ltrdigest
my $trna_file;                    # fasta file with tRNA sequences

# Counters/Index Vals
my $i = 0;                     # Array index val
my $file_num = 0;              # File number
my $proc_num = 0;              # Process number

# GFF3 output options
my $program = "ltrharvest"; 
#my $param = "default";
#my $param = "old_long";
my $param = "65_5000";
my $delim = "_";

# BOOLEANS
my $quiet = 0;
my $verbose = 0;
my $show_help = 0;
my $show_usage = 0;
my $show_man = 0;
my $show_version = 0;
my $do_test = 0;                  # Run the program in test mode
my $do_force = 0;                # For the program to try to muddle
                                 # through with unexpected errors
#-----------------------------+
# COMMAND LINE OPTIONS        |
#-----------------------------+
my $ok = GetOptions(# REQUIRED OPTIONS
		    "i|indir=s"   => \$indir,
                    "o|outdir=s"  => \$outdir,
		    # ADDITIONAL OPTIONS
		    "param=s"     => \$param,
		    "program=s"   => \$program,
		    "delim=s"     => \$delim,
		    "gt-path=s"   => \$gt_path,
		    "m|hmm=s"     => \$hmm_dir,
		    "t|trna=s"    => \$trna_file,
		    "force"       => \$do_force,
		    "q|quiet"     => \$quiet,
		    "verbose"     => \$verbose,
		    # ADDITIONAL INFORMATION
		    "usage"       => \$show_usage,
		    "test"        => \$do_test,
		    "version"     => \$show_version,
		    "man"         => \$show_man,
		    "h|help"      => \$show_help,);

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
if ( (!$indir) || (!$outdir) ) {
    print "\a";
    print STDERR "\n";
    print STDERR "ERROR: An input directory was not specified at the".
	" command line\n" if (!$indir);
    print STDERR "ERROR: An output directory was specified at the".
	" command line\n" if (!$outdir);
    print_help ("usage", $0 );
}

#-----------------------------+
# CHECK FOR SLASH IN DIR      |
# VARIABLES                   |
#-----------------------------+
# If the indir does not end in a slash then append one
# TO DO: Allow for backslash
unless ($indir =~ /\/$/ ) {
    $indir = $indir."/";
}

unless ($outdir =~ /\/$/ ) {
    $outdir = $outdir."/";
}

if ($hmm_dir) {
    unless ($hmm_dir =~ /\/$/ ) {
	$hmm_dir= $hmm_dir."/";
    }
}

#-----------------------------+
# Get the FASTA files from the|
# directory provided by the   |
# var $indir                  |
#-----------------------------+
opendir( DIR, $indir ) || 
    die "Can't open directory:\n$indir"; 
my @fasta_files = grep /\.fasta$|\.fa$/, readdir DIR ;
closedir( DIR );

my $count_files = @fasta_files;

#-----------------------------+
# SHOW ERROR IF NO FILES      |
# WERE FOUND IN THE INPUT DIR |
#-----------------------------+
if ($count_files == 0) {
    print STDERR "\a";
    print STDERR "\nERROR: No fasta files were found in the input directory\n".
	"$indir\n".
	"Fasta files must have the fasta or fa extension.\n\n";
    exit;
}

print STDERR "NUMBER OF FILES TO PROCESS: $count_files\n" if $verbose;

#-----------------------------------------------------------+
# MAIN PROGRAM BODY                                         |
#-----------------------------------------------------------+

#-----------------------------+
# Check the gt is installed
#-----------------------------+
my $gt_version;

system ("$gt_path -version");

# Check to see if the genome tools program worked
if ($? == -1) {
    print STDERR "\nThe genome program tools program failed to execute: $!\n";
    print STDERR "It is possible that the genome tools program is not".
	" installed as: $gt_path\n";
    if ($do_force) {
	print STDERR "Trying to force run.\n";
    }
    else {
	print STDERR "You can try to force run batch_ltrharvest.pl".
	    " with --force\n";
	exit;
    }
}
else {
    my $gt_version_cmd = `$gt_path -version`;
    if ($gt_version_cmd =~ /gt(.*)/ ) {
	print STDERR "\n\nGenome Tools version: $1\n";
    }
}

print STDERR "\n\n";

for my $ind_file (@fasta_files) {

    
    # Get root file name
    if ($ind_file =~ m/(.*)\.masked\.fasta$/) {
	# file ends in .masked.fasta
	$name_root = "$1";
    }
    elsif ($ind_file =~ m/(.*)\.fasta$/ ) {	    
	# file ends in .fasta
	$name_root = "$1";
    }  
    elsif ($ind_file =~ m/(.*)\.fa$/ ) {	    
	# file ends in .fa
	$name_root = "$1";
    } 
    else {
	$name_root = $ind_file;
    }

    print STDERR "\nPROCESSING: $name_root\n"
	if $verbose;

    #-----------------------------+
    # CREATE ROOT NAME DIR        |
    #-----------------------------+
    my $name_root_dir = $outdir.$name_root."/";
    unless (-e $name_root_dir) {
	mkdir $name_root_dir ||
	    die "Could not create dir:\n$name_root_dir\n"
    }

    #-----------------------------+
    # CREATE LTR_FINDER OUTDIR    |
    #-----------------------------+
    # Dir to hold gene prediction output from local software
    my $ltrharvest_dir = $name_root_dir."ltrharvest/";
    unless (-e $ltrharvest_dir) {
	mkdir $ltrharvest_dir ||
	    die "Could not create ltr_finder out dir:\n$ltrharvest_dir\n";
    }

    #-----------------------------+
    # CREATE GFF OUTDIR           |
    #-----------------------------+
    if ($do_gff_convert) {
	$gff_dir = $name_root_dir."gff/";
	unless (-e $gff_dir) {
	    mkdir $gff_dir ||
		die "Could not create gff out dir:\n$gff_dir\n";
	}
    }
    
    $file_num++;

    #-----------------------------+
    # CREATING INDEX              |
    #-----------------------------+
    # Using $ltrharvest_dir in indexname will write the index
    # to the ltrharvest dir
    my $gt_idx_root = $ind_file.".idx";
    my $gt_sfx_idx = $ltrharvest_dir.$ind_file.".idx";
    my $gt_sfx_cmd = $gt_path." suffixerator -db ".$indir.$ind_file.
	" -indexname ".$gt_sfx_idx.
	" -tis -suf -lcp -des -ssp -sds -dna";
    print STDERR "\tSUFX CMD: ".$gt_sfx_cmd."\n"
	if $verbose;

    system ($gt_sfx_cmd);

    #-----------------------------+
    # RUN LTRHARVEST              |
    #-----------------------------+
    my $gt_harvest_out = $ltrharvest_dir.$ind_file.".ltrharvest.txt";
    my $gt_harvest_gff_root = $ind_file.".ltrharvest.gff3";
    my $gt_harvest_gff = $ltrharvest_dir.$ind_file.".ltrharvest.gff3";
    my $gt_harvest_cmd = $gt_path." ltrharvest".
#	" -mintsd ".$min_tsd.
#	" -longoutput".
	" -similar 65".
	" -maxlenltr 5000".
	" -index ".$gt_sfx_idx.
	" -gff3 ".$gt_harvest_gff.
	" > ".$gt_harvest_out;
    print STDERR "\tHARV CMD: ".$gt_harvest_cmd."\n"
	if $verbose;

    system ($gt_harvest_cmd);

    #-----------------------------+
    # RUN LTRDIGEST               |
    #-----------------------------+
    # Can use wildcards as *hmm for hmms
    my $gt_digest_cmd = $gt_path." ltrdigest";
    my $gt_digest_gff = $ltrharvest_dir.$ind_file.".ltrdigest.gff3";
    if ($trna_file) {
	$gt_digest_cmd = $gt_digest_cmd.
	    " -trnas ".$trna_file;
    }
    if ($hmm_dir) {
	$gt_digest_cmd = $gt_digest_cmd.
	    " -hmms ".$hmm_dir."*.hmm --"
	}
    $gt_digest_cmd = $gt_digest_cmd.
	" ".$gt_harvest_gff.
	" ".$gt_sfx_idx.
	" > ".$gt_digest_gff;

    print STDERR "\tDGST CMT: ".$gt_digest_cmd."\n"
	if $verbose;
    system ($gt_digest_cmd);

    #-----------------------------+
    # CONVERT LTRDIGEST RESULT    |
    # TO SO COMPLIENT GFF FILE    |
    #-----------------------------+
    # Adds appropriate sequence name
    open (DIGESTIN, "<".$gt_digest_gff) ||
	next;


    my $so_gff_out = $gff_dir.$name_root.".".$param.".gff3";
    print STDERR "\tGFF3 OUT:".$so_gff_out."\n"
	if $verbose;
    open (GFFOUT, ">".$so_gff_out) ||
	next;
    
    while (<DIGESTIN>) {

	chomp;
	if ( m/^\#/ ) {
	    # We are in comment and pragma section
	    # skip as print out as is
	    print GFFOUT $_."\n";
	    next;
	}
	
	my @gff_parts = split (/\t/, $_);
	
	my $num_gff_parts = @gff_parts;
	if ( $num_gff_parts != 9) {
	    print STDERR "ERROR: GFF Num parts = ".
		$num_gff_parts."\n";
	    print STDERR "ERROR GFF: ".$_."\n";
	}

	my $type = $gff_parts[2];
	my $start =$gff_parts[3];
	my $end = $gff_parts[4];
	my $score = $gff_parts[5];
	my $strand = $gff_parts[6];
	my $phase = $gff_parts[7];
	my $attribute = $gff_parts[8];
	my $id_prefix = $name_root."_".$param."_";
	
	$attribute =~ s/Parent\=/Parent\=$id_prefix/g;
	$attribute =~ s/ID\=/ID\=$id_prefix/g;

	# Append ID to attributes without an ID
	unless ($attribute =~ m/^ID/) {
	    $attribute = "ID=".$id_prefix.$type."_".$start."_".$end.";".
		$attribute;
	}

	my $gff_out = $name_root."\t".
	    $program.$delim.$param."\t".
	    $type."\t".
	    $start."\t".
	    $end."\t".
	    $score."\t".
	    $strand."\t".
	    $phase."\t".
	    $attribute.
	    "\n";
	
	print STDERR $gff_out
	    if $verbose;
	print GFFOUT $gff_out;
#	print $_."\n" if $verbose;
	
    }
    
    close (DIGESTIN);
    close (GFFOUT);
	
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

=item * Genome Tools

The genome tools program includes the ltrfinder and ltrharvest programs.

=back

=head2 Required HMM Models

=item * PFAM Models

=item * GyDB Models

The ltrdigest proram makes use of hmm models. If you are running ltrdigest
on the

=item * tRNA Fasta file

The ltrfinder program can use a fasta file of tRNAs to identify the putative
primer binding site of the LTR retrotransposon model. Multiple programs
for mining tRNAs in genome sequence data exist. I have had good luck with
running the aragoram program. It has minimal dependencies and compiles 
easily from source code.

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
