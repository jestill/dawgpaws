#!/usr/bin/perl -w
# This is a program from renyi
# This will need to be modified to produce GFF output or to make a 
# *.genbank file that can be opened in Apollo
# Currently the output of this program has some problems with Apollo
# 1) Does not open due to the lack of a Source line with sequenced
#    span information in the genebank formatted output file
# 2) The sequences are not 
# Parse annotation result and genscan output with original seq in genbank format
# TO DO:
# - Fix GFF output
# - Make as function that accepts input output filepaths
#   OR Wirte program to accept command line input output

use strict;
use Bio::Tools::Genscan;
use Bio::SeqIO;

#my $baseDir = "/home/renyi/projects/maize_annot";
my $baseDir = "/home/jestill/projects/wheat_annotation/200702";
my $inputDir = "$baseDir/DQ537335";
my $genscanDir = "$baseDir/genscan_out";
my $outputDir = "$baseDir/genscan_parsed";
my $pepDir = "$baseDir/genscan_pep";

opendir(IN, $inputDir) or die "Cannot open input directory $inputDir: $!";

my @files = grep /fasta$/, readdir IN;

chdir($inputDir);

#-----------------------------+
# FOR EACH OF THE FASTA FILES |
#-----------------------------+
foreach (@files){
    
    my $seqin = Bio::SeqIO->new(-file=>"$inputDir/$_", -format=>'fasta');

    my $outFile = $outputDir. "/" . $_ . ".genscanparsed";
    
    my $seq = $seqin->next_seq();

    my $genscanFile = $genscanDir . "/" . $_ . ".genscan";
    my $genscanIn = Bio::Tools::Genscan->new(-file=>$genscanFile);

    my $pepFile = $pepDir . "/" . $_ . ".pep";
    my $pepOut = Bio::SeqIO->new(-file=>">$pepFile", -format=>'fasta');
    
    #-----------------------------+
    # PARSE GENSCAN FILE HERE     |
    #-----------------------------+
    # This currently uses the genscan SeqFeature format

    while(my $gene = $genscanIn->next_feature())
    {
	$seq->add_SeqFeature($gene->features_ordered());
	my $proteinSeq = $gene->predicted_protein();
	
	if(defined $proteinSeq)
	{
	    my @exons = $gene->exons_ordered();
	    my @exonRange;
	    push @exonRange, ($exons[0]->start(), $exons[-1]->end());
	    foreach my $exon(@exons)
	    {
		push @exonRange, ($exon->start(), $exon->end());
	    }
	    my $newDesc = join(" ", @exonRange);
	    $proteinSeq->desc($newDesc);
	    $pepOut->write_seq($proteinSeq);
	}else{
	    #print STDERR "protein sequences is undefined in $genscanFile\n";
	}

    } # End of while genscanIn Next feature

    #-----------------------------+
    # WRITE SEQ DATA OUT TO       |
    # GENBANK FILE FORMAT         |
    #-----------------------------+
    my $seqout = Bio::SeqIO->new(-file=>">$outFile", -format=>'genbank');
    $seqout->write_seq($seq);

}

if(0){
    opendir(OUT, $genscanDir) || 
	die "Cannot open $genscanDir: $!";

    my @genscanFiles = grep /genscan$/, readdir OUT;
    chdir $genscanDir;

    #-----------------------------+
    # FOR EACH OF THE GENSCAN     |
    # FILES                       |
    #-----------------------------+
    foreach (@genscanFiles){

	open(INF, $_) ||
	    die "Cannot open $_: $!";

	my $pepFile = $pepDir . "/" . $_ . ".pep";
	open(OF, ">$pepFile") ||
	    die "Cannot open $pepFile: $!";

	my $passHeader = 0;
	my $startSeq = 0;

	while(<INF>){
	    if(/^Predicted peptide sequence/){
		$passHeader = 1;
	    }
	    if($passHeader && /^>gi/){
		$startSeq = 1;
	    }
	    if($startSeq){
		print OF $_;
	    }
	}
    }
} # END OF IF 0 whatever the hell that is

