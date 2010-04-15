#!/usr/bin/perl -w
# Replace revision information from svn with release id

package DAWGPAWS;

#-----------------------------+
# INCLUDES                    |
#-----------------------------+
use strict;
use File::Copy;
use Getopt::Long;

my $name_root;
#my $indir = "/Users/jestill/code/dawgpaws/branches/release-1.0/docs/html/man";
my $indir = "/Users/jestill/sandbox/doc";
my $new_version = "Release 1.0";

#-----------------------------+
# CHECK FOR SLASH IN DIR      |
# VARIABLES                   |
#-----------------------------+
# If the indir does not end in a slash then append one
# TO DO: Allow for backslash
unless ($indir =~ /\/$/ ) {
    $indir = $indir."/";
}

#-----------------------------+
# Get the HTML files from the |
# directory provided by the   |
# var $indir                  |
#-----------------------------+
opendir( DIR, $indir ) || 
    die "Can't open directory:\n$indir"; 
my @html_files = grep /\.html$|\.htm$/, readdir DIR ;
closedir( DIR );

#-----------------------------+
# RUN REPEAT MASKER AND PARSE |
# RESULTS FOR EACH SEQ IN THE |
# fasta_files ARRAY FOR EACH  |
# REPEAT LIBRARY IN THE       |
# rep_libs ARRAY              |
#-----------------------------+
my $file_num=0;
for my $ind_file (@html_files) {
    
    $file_num++;
    print STDERR "Processing $file_num\n";

    # GET ROOT NAME
    if ($ind_file =~ m/(.*)\.html$/ ) {	    
	$name_root = "$1";
    }  
    elsif ($ind_file =~ m/(.*)\.htm$/ ) {	    
	$name_root = "$1";
    } 
    else {
	$name_root = "UNDEFINED";
    }

    my $cmd = "sed s//";

    my $html_file = $ind_file;
    if (-e $html_file) {
	open (INFILE, $html_file) ||
	    warn "Count not open HTML file for editing\n";
	my @html_content = <INFILE>;
	close INFILE;
	
	# ITERATE THROUGH THE RESULTING HTML TO ADD THE ANALYTICS CODE
	open (FILEOUT, ">$html_file") ||
	    die "Can not open html file:\n$html_file\n";
	foreach my $line (@html_content) {
	    $line =~ s/\$Rev.*\$/$new_version/;
	    print FILEOUT $line;
	}
	close FILEOUT;
    }

    
} # End of for each file in the input folder

exit;
