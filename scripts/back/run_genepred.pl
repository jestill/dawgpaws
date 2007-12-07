#!/usr/bin/perl -w
#-----------------------------------------------------------+
#                                                           | 
# RunFgenesh.pl                                             |
#                                                           |
#-----------------------------------------------------------+
# AUTHORS: James C. Estill, Renyi Lu                        |
# CONTACT: jestill@sourceforge.org                          |
# STARTED: 11/28/2006                                       |
# UPDATED: 04/12/2007                                       |
# DESCRIPTION:                                              |
#  Modified by Jamie Estill from original script by Renyi   |
#  Lu. Jamie added varaibles and made numerous other changes|
#  to the original script.                                  |
#                                                           |
# LICENSE:                                                  |
#  GNU General Public License, Version 3                    |
#  http://www.gnu.org/licenses/gpl.html                     |  
#                                                           |
#-----------------------------------------------------------+
# WORKING:
# Modify to run multiple gene prediction programs

#-----------------------------+
# INCLUDES                    |
#-----------------------------+
use strict;                   # Keeps thing running clean
use LWP::UserAgent;           # Download data from web sources
use Getopt::Std;              # Get options from command line

#-----------------------------+
# VARIABLES                   |
#-----------------------------+
my $debug = 1;
my $usage = "$0 <input dir> <output dir>\n";
#die $usage if(@ARGV < 2);
my $inputDir = "/home/jestill/projects/wheat_annotation/Hex/Bac2Proc";
my $outputDir = "/home/jestill/projects/wheat_annotation/Hex/fgenesh_out";

#my ($inputDir, $outputDir);
#if(@ARGV >= 1){
#    $inputDir = $ARGV[0];
#}
#if(@ARGV >= 2){
#    $outputDir = $ARGV[1];
#}

#-----------------------------+
# READ FILES FROM INPUT       |
# DIRECTORY                   |
#-----------------------------+ 
opendir(IN, $inputDir) || 
    die "Cannot open $inputDir: $!";
my @files = grep /fasta/, readdir IN;


#-----------------------------+
# DIE IF TOO MANY FILES       |
# SELECTED                    |
#-----------------------------+
if(@files > 15)
{
    die "Number of files exceeds 15 ".
	"(server's limit for noncommercial usage per day)";
}

# Temp exit to see how this is working
my $NumFiles = @files;
print "$NumFiles files were found in the input directory\n";
foreach my $indFile(@files)
{
    print "\t".$indFile."\n";
}

#exit;

my $browser = LWP::UserAgent->new();

foreach my $file(@files)
{
    print STDERR "Doing $file\n";
    my $inputFile = "$inputDir/$file";
    my $outputFile = "$outputDir/$file".".fgenesh";

    my $is_success = &fgenesh($inputFile, $outputFile);
    #my $is_success = 1;
    if(!$is_success)
    {
	if($debug)
	{
	    print STDERR "Not successful for file: $file\n";
	}
    }
    sleep 2;
}


#-----------------------------------------------------------+
#                                                           |
# SUBFUNCTIONS                                              |
#                                                           |
#-----------------------------------------------------------+

#-----------------------------+
# fgenesh subfunction as      |
# written by Renyi            |
#-----------------------------+
# Do FGENESH on one input file
# Args: 
#   (1) input file (full path)
#   (2) output file (full path)
# Return: 
#   boolean value indicating if the run is successful
sub fgenesh{
    my ($inFile, $outFile) = @_;

    print "FGENESH Processing ".$inFile."\n";

    open(OUT, ">$outFile") or die "Cannot open $outFile: $!";

    my $url = "http://www.softberry.com/cgi-bin/programs/gfind/fgenesh.pl";

    my @form_data = ('org'=>'z', 
		     'FILE'=>[$inFile]
		     );

    my $resp = $browser->post($url, \@form_data, 
			      'Content-Type'=>'form-data'
			      #,':content_file'=>$outFile
			      );
    print OUT $resp->content();
    #if($debug){
    print "Response status: ", $resp->status_line(), "\n";
    #print STDERR "Response status: ", $resp->status_line(), "\n";
    #}
    return $resp->is_success();
}


#-----------------------------------------------------------+
# NOTES                                                     |
#-----------------------------------------------------------+
# Expected parameters

# org is organism.  Should be one of:
# h => human
# d => drosophila
# n => elegans
# y => yeast
# a => plant

# Looking at the FGENESH source web interface the names of
# the organisms may have changes
# Monocot.dat    - Monocots (Corn,Rice,Wheat,Barley)
# 


# program_name
# fgenesh => FGENESH / HMM based Human/Drosophila/Nematode/Plan (multiple genes, both chains)
# Fgenesh => FGENESH (with Donor GC)/ HMM based Human/Drosophila/Nematode/Plant (multiple genes, both chains)
# fgenes => FGENES / Pattern based Human Gene structure prediction (multiple genes, both chains)
# fgenesm => FGENES-M / Pattern based Human Multiple variants of Gene structure prediction)
# est    => BESTORF - Finding potential coding fragment in EST/mRNA (Human/Plant)
# fgene  => FGENE - Finding gene by exon search and assembling
# fex    => FEX - Finding potential 5'-, internal and 3'-coding exons
# spl    => SPL / search for potential splice sites
# cdsb   => CDSB / search for  E.coli protein coding genes
# nsite  => NSITE / Recognition of Regulatory motifs with statistics
# polyah => POLYAH / Recognition of 3'-end cleavage and polyadenilation region
# tssg   => TSSG / Recognition of human PolII promoter region and start of transcription
# tssw   => TSSW / Recognition of human PolII promoter region and start of transcription</b>
# rnaspl => RNASPL / search for exon-exon junction positions in cDNA
# hbr    => HBR / recognition of Human and E.coli sequences
