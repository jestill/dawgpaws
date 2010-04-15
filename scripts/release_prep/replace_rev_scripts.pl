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
my @perl_files = grep /\.pl$|\.perl$/, readdir DIR ;
closedir( DIR );

#-----------------------------+
# RUN REPEAT MASKER AND PARSE |
# RESULTS FOR EACH SEQ IN THE |
# fasta_files ARRAY FOR EACH  |
# REPEAT LIBRARY IN THE       |
# rep_libs ARRAY              |
#-----------------------------+
my $file_num=0;
for my $ind_file (@perl_files) {
    
    $file_num++;

    # GET ROOT NAME
    if ($ind_file =~ m/(.*)\.pl$/ ) {	    
	$name_root = "$1";
    }  
    elsif ($ind_file =~ m/(.*)\.perl$/ ) {	    
	$name_root = "$1";
    } 
    else {
	$name_root = "UNDEFINED";
    }

    print STDERR "Processing $file_num\n";

    my $cmd = "sed s//";
    
    # REPLACE REV
    my $perl_file = $indir.$ind_file;
    if (-e $perl_file) {
	open (INFILE, $perl_file) ||
	    warn "Count not open script for editing\n";
	my @html_content = <INFILE>;
	close INFILE;
	
	# ITERATE THROUGH THE RESULTING HTML TO ADD THE ANALYTICS CODE
	open (FILEOUT, ">$perl_file") ||
	    die "Can not open html file:\n$perl_file\n";
	foreach my $line (@html_content) {
	    $line =~ s/\$Rev.*\$/$new_version/;
	    print FILEOUT $line;
	}
	close FILEOUT;
    }
    
    # REPLACE VERSION
    my $perl_file = $indir.$ind_file;
    if (-e $perl_file) {
	open (INFILE, $perl_file) ||
	    warn "Count not open script for editing\n";
	my @html_content = <INFILE>;
	close INFILE;
	
	my $new_perl_version = "my \(\$VERSION\) \= \"Release 1.0\"\;";
#	my $new_perl_version = 'my ($VERSION) = Release 1.0';
	# ITERATE THROUGH THE RESULTING HTML TO ADD THE ANALYTICS CODE
	open (FILEOUT, ">$perl_file") ||
	    die "Can not open html file:\n$perl_file\n";
	foreach my $line (@html_content) {
	    $line =~ s/\my \(\$VERSION\).*\;/$new_perl_version/;
	    print FILEOUT $line;
	}
	close FILEOUT;
    }

    
} # End of for each file in the input folder

exit;
