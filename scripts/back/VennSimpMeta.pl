#!/usr/bin/perl -w

print "The Venn meta program has started\n";

# Just a bunch of system commands to do some VennSeq runs
my $prog = "/home/jestill/projects/wheat_annotation/scripts/VennSeq.pl";
my $BaseDir = "/home/jestill/projects/wheat_annotation/VennTest/";
#my $BaseDir = "/home/jestill/projects/wheat_annotation/VennEst/";

print "processing ".$BaseDir."\n";

#opendir(INDIR, $BaseDir) ||
#    die "Can not open input directory:\n$FeatDir: $!";
@Files2Proc = (
	       "HEX1057J19",
	       "HEX1057D14",
	       "HEX1057D07",
	       "HEX1057C05",
	       "HEX1057B03",
	       "HEX0961H02",
	       "HEX0961G12",
	       "HEX0961F17",
	       "HEX0961F07",
	       "HEX0961E20",
	       "HEX0593O07",
	       "HEX0075G21"
	       );
    
my $Num2Proc = @Files2Proc;
my $Count=0;
foreach $FileRoot (@Files2Proc)
{
    $Count++;
    my $InFile = $BaseDir.$FileRoot.".fasta.masked";
    unless (-e $InFile) {print "Infile does not exist:\n$InFile\n";exit;}
    my $OutFile = $BaseDir.$FileRoot.".venn.txt";
    my $DatDir = $BaseDir.$FileRoot."/";
    my $ProgOut = $BaseDir.$FileRoot.".prog.out";
#    my $SysCmd = "$prog -i $InFile -o $OutFile -d $DatDir -q -f tab";
    print "processing $Count of $Num2Proc $FileRoot\n";
    my $SysCmd = "$prog -i $InFile -o $OutFile -d $DatDir -q -f tab > $ProgOut";
#my $SysCmd = "$prog";
#    print "The system command is:\b$SysCmd\n";
    system ($SysCmd);
}
