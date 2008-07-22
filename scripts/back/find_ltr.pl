#!/usr/bin/perl
#-----------------------------------------------------------+
#                                                           |
# find_ltr.pl - Identify LTR Retrotransposons               |
#                                                           |
#-----------------------------------------------------------+
#                                                           |
#  AUTHOR: Mina Rho                                         |
# STARTED: Unknown                                          |
# UPDATED: 09/12/2007                                       |
# CONTRIBUTORS:                                             |
#   James C. Estill - Modifications to make variables       |
#                     easier to use and cleaned up paths.   |
#                                                           |
# DESCRIPTION:                                              |
#   This is J. Estill's modification of the find_ltr.pl     |
#   program.                                                |
#-----------------------------------------------------------+

# outdir not currently working correctly for long file paths
# Will not worry about for now - JCE 09/13/2007

use strict;
use Getopt::Long;

#-----------------------------+
# SET PATH VARIABLES          |
#-----------------------------+
# USE ENVIRONMETNAL VARIABLE FOR PATH ROOT
my $path_tool = $ENV{FIND_LTR_ROOT} || "";

if ($path_tool) {
    unless ($path_tool =~ /\/$/){
	$path_tool = $path_tool."/";
    }
}

my $errmsg;                    # Var to hold error messages

# PATH VARS
my $outdir;                    # Optional path to an outdir root
my $name_root;                 # Used to refer to seqfile base name
my $mask_file;                 # Seq file masked by TRF
my $mem_file;                  # File containing MEMs
my $dist_file;
my $bin_file;
my $ltrpair_file;
my $ltrbigpair_file;
my $ltrseq_file;
my $ltrpos_file;

#-----------------------------+
# MODEL VARIABLES             |
#-----------------------------+
# Jamie Modified command line to accept
# modifications of the following variables
my $min_len_mem = 40;          # minimal length of the exact match
                               # MEM -- Maximal Exact Match
                               # Var sent to the ltr program
my $min_dist = 1100;           # Minimum distance between MEMs 
                               # Var sent to the ltr program
my $max_dist = 16000;          # Max distance between MEMs
                               # Var sent to the ltr program
my $max_gap_prev_next = 40;    # gap between exact matches to combine them, 
my $min_sub_len_ltr = 0;       # VARIABLE DOES NOT APPEAR TO BE USED
my $max_sub_len_ltr = 1000;    # VARIABLE DOES NOT APPEAR TO BE USED
my $min_len_ltr = 100;         # Minimum length of LTR  
                               # Used in the make_pair subfunction    
my $max_len_ltr = 1000;        # Maximum length of LTR
                               # Used in the make_pair subfunction
my $range_bin = 500;           # range in the bin

my $min_identity = 80;         # Minimum Percent identify between merged
                               # MEMs
# ORF VARS
my $min_len_orf = 700;         # length of ORF
my $e_val = "0.000001";        # E value for hmmsearch

#-----------------------------+
# PFAM HMM MODELS             |
#-----------------------------+
# This is from PFAM 19
# These are used in the find_doman subfunction
my @pf=(# REVERSE TRANSCRIPTASE
        "PF00078_fs.hmm",
	"PF07727_fs.hmm",
	# PAO RETRTROTRANSPOSON PEPTIDASE
	"PF05380_fs.hmm",
	# INTEGRASE CORE DOMAIN
	"PF00665_fs.hmm",
	# INTEGRASE DNA BINDING DOMAIN
	"PF00552_fs.hmm", 
	# INTEGRASE ZINC BINDING DOMAIN
	"PF02022_fs.hmm",
	# EUKARYOTIC ASPARTYL PROTEASE
	"PF00026_fs.hmm",
	# RETROVIRAL ASPARTYL PROTEASE
	"PF00077_fs.hmm",
	# RNASEH H
	"PF00075_fs.hmm",);

        #-----------------------------+
        # THE FOLLOWING COMMENTED OUT IN THE 
        # ORIGINAL DOWNLOAD OF find_ltr
        #-----------------------------+
        # RETROTRANSPSON GAG PROTEIN
	#"PF03732_fs.hmm",
        # ENDONUCLEASE/EXONUCLEASE PHOSPHATASE
	#"PF03372_fs.hmm",
        # ENDONUCLEASE I 
        #"PF04231_fs.hmm",
        # ENDONUCLEASE V
        #"PF04493_fs.hmm",
        # ENDONUCLEASE VII
        #"PF02945_fs.hmm",
        # DNA/RNA NON-SPECIFIC ENDONUCLEASE
        # "PF01223_fs.hmm");

# BOOLEANS
my $verbose = 0;

# PATH VARIABLES

my $tool_trf="trf400.linux.exe"; # Tandem Repeat Finder (TRF) name
#my $path_tool;                  # Base path of the external software tools
my $genome_file;                 # Set in get_input subfunction
my $genome_seq;

#-----------------------------+
# GET INPUT FROM COMMAND LINE |
#-----------------------------+
get_input();

#-----------------------------+
# SET OUTPUT FILE PATHS       |
#-----------------------------+
if ($outdir) {

    # Append / to end of outdir if needed
    unless ($outdir =~ /\/$/ ) {
	$outdir = $outdir."/";
    }

    #-----------------------------+
    # MAKE OUTPUT DIR             |
    #-----------------------------+
    unless (-e $outdir) {
	mkdir ($outdir, 0777) ||
	    die "Can not make output directory $outdir\n"; 
    }

    #-----------------------------+
    # GET BASE FILE NAME          |
    #-----------------------------+
    # Supported fasta extensions are:
    #  fasta | fast | fa | tfa | fsa

    # If full path address
    if ($genome_file =~ m/.*\/(.*)\.fasta$/
	|| $genome_file =~ m/.*\/(.*)\.fast$/
	|| $genome_file =~ m/.*\/(.*)\.fa$/  
	|| $genome_file =~ m/.*\/(.*)\.tfa$/ 
	|| $genome_file =~ m/.*\/(.*)\.fsa$/ ) {
	$name_root = "$1";
    }
    # If short path address
    elsif ($genome_file =~ m/(.*)\.fasta$/ 
	   || $genome_file =~ m/(.*)\.fast$/
	   || $genome_file =~ m/(.*)\.fa$/ 
	   || $genome_file =~ m/(.*)\.tfa$/
	   || $genome_file =~ m/(.*)\.fasta$/) {
	$name_root = "$1";
    }
    else {
	die "ERROR: Can not extract base file name from:\n$name_root:\n";
    }
    
    #-----------------------------+
    # SET FILE PATHS              |
    #-----------------------------+ 
    $mask_file       = $outdir.$name_root.".mask";
    $mem_file        = $outdir.$name_root.".mem";
    $dist_file       = $outdir.$name_root.".dist";
    $bin_file        = $outdir.$name_root.".bin";
    $ltrpair_file    = $outdir.$name_root.".ltrpair";
    $ltrbigpair_file = $outdir.$name_root.".ltrbigpair";
    $ltrseq_file     = $outdir.$name_root.".ltrseq";
    $ltrpos_file     = $outdir.$name_root.".ltrpos";
}
else {
     # All of the following are extensions added to the
     # path of the genome file as set by seq at the command line --seq

    $mask_file       = $genome_file.".mask";
    $mem_file        = $genome_file.".mem";
    $dist_file       = $genome_file.".dist";
    $bin_file        = $genome_file.".bin";
    $ltrpair_file    = $genome_file.".ltrpair";
    $ltrbigpair_file = $genome_file.".ltrbigpair";
    $ltrseq_file     = $genome_file.".ltrseq";
    $ltrpos_file     = $genome_file.".ltrpos";
}



# Show variables if running in verbose mode
if ($verbose) {
    print STDOUT "\n";
    print STDOUT "-----------------------------------\n";
    print STDOUT " PROGRAM VARIABLES\n";
    print STDOUT "-----------------------------------\n";
    print STDOUT "EVAL:\t\t$e_val\n";
    print STDOUT "MIN-ORF:\t$min_len_orf\n";
    print STDOUT "MIN-ID:\t\t$min_identity\n";
    print STDOUT "\n";
    print STDOUT "-----------------------------------\n";
    print STDOUT "\n";

    print STDOUT "-----------------------------------\n";
    print STDOUT " FILE PATHS\n";
    print STDOUT "-----------------------------------\n";
    print STDOUT "TOOL DIR:\t$path_tool\n\n";
    print STDOUT "NAME ROOT:\t$name_root\n" if $outdir;
    print STDOUT "MASK:\t\t$mask_file\n";
    print STDOUT "MEM:\t\t$mem_file\n";
    print STDOUT "DIST:\t\t$dist_file\n";
    print STDOUT "BIN:\t\t$bin_file\n";
    print STDOUT "LTRPAIR:\t$ltrpair_file\n";
    print STDOUT "LTRBIGPAIR:\t$ltrbigpair_file\n";
    print STDOUT "LTRSEQ:\t\t$ltrseq_file\n";
    print STDOUT "LTRPOS:\t\t$ltrpos_file\n";
    print STDOUT "-----------------------------------\n";
    print STDOUT "\n";

}

# Temp exit to check output file paths
#exit;

#-----------------------------+
# GET SEQUENCE FILE AND LOAD  |
# TO genome_seq STRING        |
#-----------------------------+
print STDERR "Fetching sequence ..\n" if $verbose;
get_sequence();

#-----------------------------+
# MASK SEQS WITH TRF          |
#-----------------------------+
print STDERR "Masking with Tandem Repeat Finder ..\n" if $verbose;
system($path_tool.$tool_trf." ".$genome_file.
       " 2 7 7 80 10 50 500 -m -h  > /dev/null 2>&1");

my @temp = split(/\//, $genome_file.".2.7.7.80.10.50.500");

system("mv ".$temp[$#temp].".mask ".$mask_file);

system("rm ".$temp[$#temp].".dat");

#-----------------------------+    
# RUN LTR TO GET MEM          |
#-----------------------------+
# LTR is the C program to find the Maximal Exact Matches
# LTR is the program taht uses a suffix array data structure
# LTR is a modified GAME
print STDERR "Finding MEMs with ltr ..\n" if $verbose;
system($path_tool."ltr -i ".$mask_file." -o ".$mem_file." -s ".
       $min_len_mem." -d ".$min_dist." -D ".$max_dist."  > /dev/null 2>&1");

system("rm ".$mask_file);

# THROW NO MEM ERROR if no output found from the LTR program
if(!(-e $mem_file)){	
    die	"ERROR: No MEM\n";
}
elsif(-s $mem_file==0){
    die "ERROR: No MEM within the distance \n";
}

#-----------------------------+
# MAKE BINS                   |
#-----------------------------+
print STDERR "Joining MEMs ...\n" if $verbose;
make_bin();

#-----------------------------+
# PAIR MEMs                   |
#-----------------------------+
print STDERR "Pairing MEM ...\n" if $verbose;
make_pair();

#-----------------------------+
# CHAIN PAIRS                 |
#-----------------------------+
print STDERR "Chaining pairs ...\n" if $verbose;
chain_pairs(); 

#-----------------------------+
# CHECK IDENTITY              |
#-----------------------------+
print STDERR "Checking identity ...\n" if $verbose;
check_identity();


exit;

#-----------------------------------------------------------+
# SUBFUNCTIONS                                              |
#-----------------------------------------------------------+

sub get_input{

    my ($path_root, $seq_name);
    GetOptions(# REQUIRED ARGUMENTS
	       'seq=s'       => \$seq_name,
	       # ADDITIONAL OPTIONS
	       'outdir=s'    => \$outdir,            # Base outdir
	       'root=s'      => \$path_root,	     # Dir of the BINARIES
	       # MEM Options
	       'min-mem=i'   => \$min_len_mem,       # Default is 40
               'min-dist=i'  => \$min_dist,          # Default is 1100
               'max-dist=i'  => \$max_dist,          # Default is 16000
               'max-gap=i'   => \$max_gap_prev_next, # Default is 40
	       # LTR Options
               'min-ltr=i'   => \$min_len_ltr,       # Defulat is 100
               'max-ltr=i'   => \$max_len_ltr,       # Default is 1000
               'range-bin=i' => \$range_bin,         # Default is 500
	       # ORF options
               'min-orf=i'   => \$min_len_orf,       # Default is 700
	       'e-val=s'     => \$e_val,             # Default is 0.000001
	       # BOOLEANS
	       'verbose'     => \$verbose,
	       );
    
    # Jamie changing the $genome_file variable
    # to not be required to be in test/ dir
    # I also changed the the find_ltr binaries to a
    # dir other then a tool subdir 
    if ($path_root) {
	if ($path_root =~ /\/$/){
	    $path_tool = $path_root;
	}
	else
	{
	    $path_tool = $path_root."/";
	}
    } # End of if path_root

    # Genome file is just the seq_name path as passed at the command line
    $genome_file = $seq_name;

    if (! -e $genome_file){
	print "The file $genome_file does not exist.\n";
	usage();
    }

    #-----------------------------+
    # CHECK FOR EXISTENCE OF THE  |
    # REQUIRED EXTERNAL SOFTWARE  |
    # TOOLS                       |
    #-----------------------------+
    
    if (!-e $path_tool."trf400.linux.exe"){
	$errmsg =  "ERROR: The program TRF is not correctly installed ".
	    "in $path_tool directory.";
	die "$errmsg\n";
    }
    if (!-e $path_tool."stretcher"){
	$errmsg = "ERROR: The program stretcher is not correctly installed ".
	    "in $path_tool directory.";
	die "$errmsg\n";
    }
    if (!-e $path_tool."transeq"){
	$errmsg = "ERROR: The program transeq is not correctly installed in ".
	    "$path_tool directory.";
	die "$errmsg\n";
    }
    if (!-e $path_tool."hmmsearch"){
	$errmsg =  "ERROR: The program hmmsearch is not correctly installed ".
	    "in $path_tool directory.";
	die "$errmsg\n";
    }
} # End of get_input subfunction


sub usage {
    die "Usage: find_ltr.pl --root=<main_dir_path> --seq=<seq_name>" ;
}


sub get_sequence{

    open(INPUT, $genome_file)|| die("Couldn't open $genome_file!\n");
    foreach my $each_line (<INPUT>)  {

        if ($each_line =~ m/>/){
            $genome_seq ="";
        }else{
            chomp($each_line);
            $genome_seq .= $each_line;
        }
    }
    close(INPUT);
}

sub make_bin {

    system("sort -k 4 -n $mem_file > $dist_file");
    system("rm $mem_file");
    open(DAT, $dist_file)|| die("ERROR: Couldn't open $dist_file!\n");
    open OUTPUT, ">$bin_file";

    my (@line_array, @line_array_bin);
    my ($prev1, $prev2);
    my $flag = 0;    
    my $temp_file1 = $genome_file."_temp1";
    my $temp_file2 = $genome_file."_temp2";
    my $start_dist = $min_dist;

    my $line_count=0;
    print OUTPUT "------------------------------------------".
	"------------------------------\n";
    foreach my $each_line (<DAT>) {

	@line_array = split(/\t/, $each_line);
	
	if ($flag==0)   {        # to print out the first line
	    open TEMP, ">$temp_file1";
	    print TEMP $each_line;
	    $flag=1;
	    
	}elsif ($line_array[3] < ($start_dist+$range_bin) && $flag==1) {  # for each bin, add lines in the bin into a temp file(temp1)
	    
	    print TEMP $each_line;
	    
	} else  { #at the last line of the bin, sort them by the start position of mem and put them into a temp file(temp2)

	    $prev1 = 0;   #the ending position of previous MEM
	    $prev2 = 0;

            #after sorting each file of a bin, add it into bin_file
	    close(TEMP);
	    system ("sort +0 -1n +1 -2n $temp_file1 > $temp_file2");       
	    #save the sorted lines in a bin into bin file
	    open (TEMP_SORTED, $temp_file2) 
		|| die("Couldn't open sorted bin file!\n");
	    foreach my $each_line_bin (<TEMP_SORTED>) {

		chop($each_line_bin);
		@line_array_bin = split(/\s+/, $each_line_bin);
		
		$line_count++;

		
		print OUTPUT $line_array_bin[0]."\t".
		    $line_array_bin[1]."\t".$line_array_bin[2].
		    "\t".$line_array_bin[3]."\t";
		print OUTPUT ($line_array_bin[0]-$prev1)."\t".
		    ($line_array_bin[1]-$prev2)."\n";
		$prev1 = $line_array_bin[0]+$line_array_bin[2];
		$prev2 = $line_array_bin[1]+$line_array_bin[2];

	    }
	    close(TEMP_SORTED);

            #open a temp file again for the next bin
	    print OUTPUT "-------------------------------------------".
		"-----------------------------\n";
	    	    
	    open TEMP, ">$temp_file1";
	    print TEMP $each_line;
	    
	    @line_array =  split(/\s+/,$each_line);
	    $start_dist = $start_dist + $range_bin;
	    
	    while ($start_dist+$range_bin <= $line_array[3])  {

		$start_dist = $start_dist + $range_bin;
		
	    }
	    $flag =1;
	}	
    }    
    close(DAT);
    close(TEMP);
    close(TEMP_SORTED);
    close(OUTPUT);

    system("rm $temp_file1");
    system("rm $temp_file2");
    system("rm ".$dist_file);
    if ($line_count == 0){
	system("rm ".$bin_file);
	die ("ERROR: no bin");
    }
}

sub make_pair {

    open(DAT, $bin_file)|| die("Couldn't open $bin_file!\n");
    my $ltrpair_file_temp = $genome_file."_temp";
    open OUTPUT_POS, ">$ltrpair_file_temp";

    # Transform input
    my $count=0;
    my $real_count=0;
    my $start_flag=0;
    my $line_count=0;
    my ($start1, $end1, $start2, $end2);
    my @line_array;

    foreach my $each_line (<DAT>) {
	
	if ($each_line =~ m/^--/) {
	    
	    $start_flag=1;
	    
	    if ($count!=0 && ($end1-$start1)>$min_len_ltr && ($end1-$start1)<$max_len_ltr  && 
		($end2-$start2)>$min_len_ltr  && ($end2-$start2)<$max_len_ltr && 
		(($start2-$end1)>$min_dist) && (($start2-$end1)<$max_dist) ) {
		
		$real_count++;
		print OUTPUT_POS $real_count."\t".
		    $start1."\t".$end1."\t".$start2."\t".$end2."\t";
		print OUTPUT_POS eval($end1-$start1)."\t".
		    eval($end2-$start2)."\t".eval($start2-$end1)."\n";
		$line_count++;
	    }
	}else {
	    chop($each_line);
	    @line_array = split(/\t/, $each_line);
	    
	    if ($start_flag==1) { #the first line of each bin
		
		$count++;
		$start1 = $line_array[0];
		$start2 = $line_array[1];
		$end1 = $line_array[0]+$line_array[2];
		$end2 = $line_array[1]+$line_array[2];
		$start_flag=0;
		
	    }elsif (($line_array[4]<$max_gap_prev_next) && ($line_array[4]>-5) && ($line_array[5]<$max_gap_prev_next) && ($line_array[5]>-5) ){
		
		#inside the cluster
		$end1 = $line_array[0]+$line_array[2];
		$end2 = $line_array[1]+$line_array[2];
		
	    }else {
		
		#start of the cluster
 		if ($count!=0 && ($end1-$start1)>$min_len_ltr && ($end1-$start1)<$max_len_ltr  && 
		    ($end2-$start2)>$min_len_ltr  && ($end2-$start2)<$max_len_ltr && 
		    (($start2-$end1)>$min_dist) && (($start2-$end1)<$max_dist) ) {
		
		    $real_count++;
		    print OUTPUT_POS $real_count."\t".$start1."\t".$end1."\t".$start2."\t".$end2."\t";
		    print OUTPUT_POS eval($end1-$start1)."\t".eval($end2-$start2)."\t".eval($start2-$end1)."\n";
		    $line_count++;
		}		
		$count++;
		$start1 = $line_array[0];
		$start2 = $line_array[1];
		$end1 = $line_array[0]+$line_array[2];
		$end2 = $line_array[1]+$line_array[2];
	    }
	}
    }
    close(DAT);
    close(OUTPUT_POS);
    system("rm ".$bin_file);
    if ($line_count==0){
      	system("rm ".$ltrpair_file_temp);
	die ("ERROR: no bin");
    }else {
	system("sort +1 -2n +3 -4n ".$ltrpair_file_temp." > ".$ltrpair_file);
	system("rm ".$ltrpair_file_temp);
    }
}


sub chain_pairs {

    open(DAT, $ltrpair_file)|| die("Couldn't open  file!\n"); #sub pair file
    open OUTPUT, ">$ltrbigpair_file";
    
    my $line_count=1;
    my $first = 1;
    my @temp=();
    my ($start1, $start2, $end1, $end2);

    foreach my $each_line (<DAT>){

	chop($each_line);

	@temp = split(/\t/, $each_line);

	if ($first == 1) {

	    $start1 = $temp[1];
	    $start2 = $temp[3];
	    $end1 = $temp[2];
	    $end2 = $temp[4];
	    $first = 0;

	} else {

	    if (($temp[1] - $end1 < $max_gap_prev_next) && ($temp[3] - $end2 < $max_gap_prev_next)) {
		
		$end1 = $temp[2];
		$end2 = $temp[4];

	    }else {
		
		if (($end1-$start1)>$min_len_ltr && ($end1-$start1)<$max_len_ltr  && 
		    ($end2-$start2)>$min_len_ltr  && ($end2-$start2)<$max_len_ltr && 
		    (($start2-$end1)>$min_dist) && (($start2-$end1)<$max_dist) ) {

		    print OUTPUT $line_count."\t".$start1."\t".$end1."\t".
			$start2."\t".$end2."\n";
		    $line_count++;
		}

		$start1 = $temp[1];
		$start2 = $temp[3];
		$end1 = $temp[2];
		$end2 = $temp[4];	
	    }
	}
    }
    if (($end1-$start1)>$min_len_ltr && ($end1-$start1)<$max_len_ltr  && 
	($end2-$start2)>$min_len_ltr  && ($end2-$start2)<$max_len_ltr && 
	(($start2-$end1)>$min_dist) && (($start2-$end1)<$max_dist) ) {

	print OUTPUT $line_count."\t".$start1."\t".$end1."\t".
	    $start2."\t".$end2."\n";
	$line_count++;
    }
    close(DAT);
    close(OUTPUT);

    if ($line_count==0){
      	system("rm ".$ltrpair_file);
      	system("rm ".$ltrbigpair_file);
	die ("ERROR: no pair");
    }else{
      	system("rm ".$ltrpair_file);
    } 
}

sub check_identity {


    # If the seqs pass this test, they will
    # be passed on to the output files
    my $input1 = $genome_file.".temp1";
    my $input2 = $genome_file.".temp2";
    my $input3 = $genome_file.".temp3";
    my $input4 = $genome_file.".temp4";

    my $count=0;

    my ($seq1, $seq2, $seq3, $similarity, $direction, $is_ltr );
    my (@plus_minus, @long_orf, @pre, @temp);

    if (-e $ltrpos_file){
	system("rm $ltrpos_file");
    }
    if (-e $ltrseq_file){
	system("rm $ltrseq_file");
    }

    open OUTPUT, ">$ltrpos_file";   
    open OUTPUT_SEQ, ">$ltrseq_file";

    open(DAT, $ltrbigpair_file)
	|| die("Couldn't open $ltrbigpair_file!\n"); 
    @pre = (0,0,0,0,0);
    foreach my $each_line (<DAT>){

	chop($each_line);
	@temp = split(/\t/, $each_line);
	
	$seq1 = substr($genome_seq, $temp[1], $temp[2]-$temp[1]);
	$seq2 = substr($genome_seq, $temp[3], $temp[4]-$temp[3]);
	$seq3 = substr($genome_seq, $temp[2]+1, $temp[3]-$temp[2]);

	open SEQ1, ">$input1";
	open SEQ2, ">$input2";
	open SEQ3, ">$input3";
	
	print SEQ1 ">".$temp[0].":1:".$temp[1].":".$temp[2]."\n";
	print SEQ1 $seq1."\n";
	print SEQ2 ">".$temp[0].":2:".$temp[3].":".$temp[4]."\n";
	print SEQ2 $seq2."\n";
	print SEQ3 ">".$temp[0].":rt:".$temp[2].":".$temp[3]."\n";
	print SEQ3 $seq3."\n";
	close(SEQ1);
	close(SEQ2);
	close(SEQ3);
	
	$similarity = find_sim($input1, $input2);

	# Similarity of 80 percent here
	
	#if ($similarity >80)
	# Value was hard coded to 80, JCE changed to below
	if ($similarity >$min_identity)
	{
	    @long_orf=(0,"*");
	    $is_ltr=0;
	    $direction="";

	    @plus_minus = find_domain($input3, $input4);
	    # find_domain[0] is plus strand e value
	    # find_domain[0] is minus strand e value

	    # ANOTHER E VALUE HERE
	    #if ($plus_minus[0] < $plus_minus[1] && $plus_minus[0] < 1e-10){ 
            # JCE CHANGED hard coded EVALUE TO Variable
	    if ($plus_minus[0] < $plus_minus[1] && $plus_minus[0] < $e_val){ 
		# forward strand
		$direction="+";
		$is_ltr=1;

	    }
	    #elsif($plus_minus[1] < $plus_minus[0] && $plus_minus[1] < 1e-10){
	    # JCE Chagned hard coded EVALUE to Variable
	    elsif($plus_minus[1] < $plus_minus[0] && $plus_minus[1] < $e_val){
		# Backward strand
		$direction="-";
		$is_ltr=1;
		
	    }
	    else{
		# $long_orf[0]:length, $long_orf[1]:strand
		@long_orf = check_long_orf($input4);
		
		if ($long_orf[0]>$min_len_orf ){
		    
		    $direction=$long_orf[1];
		    $is_ltr=1;
		}
	    }

	    # is_ltr is boolean based on long orf or evalu threshold
	    # temp[1] is data split from each line of the
            # #ltrbigpair_file is $genome_file.ltrbigpair
            # $temp[1] is start of the of the 5' LTR ?
	    if ($is_ltr==1 && $pre[4] <$temp[1]){

		$count++;
		# PRINT SUMMARY OUTPUT FOR LTR POSITIONS,SIMILARITY
		print OUTPUT $count."\t".$temp[1]."\t".$temp[2].
		    "\t".$temp[3]."\t".$temp[4]."\t".$direction."\t";
		print OUTPUT  eval($temp[2]-$temp[1]).
		    "\t".eval($temp[4]-$temp[3])."\t".
		    eval($temp[3]-$temp[2])."\t".$similarity."\n";

		# PRINT OUT SEQUENCES TO SEQ FILE
		print OUTPUT_SEQ ">".$count."_1\n";
		print OUTPUT_SEQ $seq1."\n";
		print OUTPUT_SEQ ">".$count."_2\n";
		print OUTPUT_SEQ $seq2."\n";

		@pre = @temp;

	    }
	}
    }
    
    close(DAT);
    close(OUTPUT);

    print "\nLTR Retrotransposons discovered: $count\n" if $verbose;
    
    # REMOVE ALL TEMPORARY FILES
    # This could be made a var to just move these somewhere
    # instead of deleting the data
    system("rm ".$input1);
    system("rm ".$input2);
    system("rm ".$input3);
    system("rm ".$input4);
    system("rm ".$ltrbigpair_file);

}

sub find_sim{

    # Find similarity between two sequences passed to the subfunction
    # using the stretcher program

    my $result;
    my $temp_matcher = $path_tool."stretcher -asequence=".$_[0]." -bsequence=".$_[1]." -outfile=stdout -awidth3=4000 2>/dev/null";
    my $str_result = `$temp_matcher`;

    if ($str_result =~ /(Identity:\s+\d+\/\d+\s*\((\d+.\d+)%\))/){
	$result=$2;
    }
    return $result;
}

sub find_domain{

    system($path_tool."transeq -sequence ".$_[0]." -outseq ".
	   $_[1]." -frame=6  2>/dev/null");
    
    my $plus = 1;
    my $minus = 1;
  
    # FOR EACH MODEL IN THE ARRAY 
    # RUN HMMSEARCH
    for (my $j=0; $j<=$#pf; $j++){
		
	# E-value is 1 x e-6
	# AS ORIGINALLY CODED
	#my $temp_tool=$path_tool."hmmsearch -E 0.000001 ".
	my $temp_tool=$path_tool."hmmsearch -E "."$e_val ".
	    $path_tool.$pf[$j]." ".$_[1];
	my $str = `$temp_tool`;
	
	if ($str =~ /\s---\n(\d+\_\d\s+\d+\.\d+\s+((\d|\-|\.|e)+))\s/){
	    my @temp_plus = split(/\s+/, $1);
	    
	    if ($temp_plus[0]=~ /(\_1|\_2|\_3)/){
		$plus = $plus * $temp_plus[2];
	    }else{
		$minus = $minus * $temp_plus[2];
	    }
	}
    }
    return ($plus, $minus);
}


sub check_long_orf{
    
    open(DAT, $_[0])||die("ERROR: Couldn't open $_[0]\n");
    my @temp = <DAT>;
    chomp(@temp);

    my $temp_frame;
    my @max_len;
    my $max_frame;
    $max_len[0]=0;
    my $i=0;
    my $long = join("",@temp);
    while ($long=~ m/\>\d+_\d+((\w|\*)+)/g){
	$temp_frame = $1;
	$i++;
	while($temp_frame =~ m/(\w+)/g){
	    if (length($1) > $max_len[0]){
		$max_len[0] = length($1);
		$max_frame = $i;
	    }
	}
    }
    if ($max_frame == 1 || $max_frame ==2 || $max_frame ==3){
	$max_len[1] = "+";
    }else {
	$max_len[1] = "-";
    }

    return @max_len;
}

=head1 NAME

find_ltr.pl - LTR Finding program

=head1 SYNOPSIS

  USAGE: find_ltr.pl --seq InSeq.fasta 

    REQUIRED ARGUMENTS
    --seq             # Fasta formatted sequence file containg a single record
    ADDITIONAON OPTIONS
    --root            # Directory containing the required external tools
    MEM OPTIONS
    --min-mem         # Minimum length of the Maximal Exact Match (MEM) 
    --min-dist        # Minimum distance between MEMs
    --max-dist        # Maximum distance between MEMs
    --max-gap         # Maximum gap between MEMs
    LTR OPTIONS
    --min-ltr         # Minimum length of the LTR
    --max-ltr         # Maximum length of the LTR
    --range-bin       # 
    # ORF OPTIONS
    --min-orf         # Minimum length of ORF
    --e-val           # Max E value of HMMER result
    # BOOLEANS
    --verbose         # Run program in verbose mode

=head1 AUTHOR

Mina Rho

=head1 CONTRIBUTORS

James Estill

=cut

#-----------------------------+
# CHANGELOG                   |
#-----------------------------+
# 05/24/2007
# - Got rid of line wraps where possible to maintain 
#   the standard 80 character column width
# 09/12/2007
# - Set model varaibles at the command line
# - Made Evalue of hmm a variable
# - Made identity a variable, default is 80%
# - Added ENV{FIND_LTR_ROOT} for the path root
# - Output is now in the same directory as the 
#   input file. Not in tools.
# 09/13/2007
# - Adding some POD documentation 
