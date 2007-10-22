#!/usr/bin/perl -w
#-----------------------------------------------------------+
#                                                           |
# dir_merge.pl - Merge directories                          |
#                                                           |
#-----------------------------------------------------------+
#                                                           |
#  AUTHOR: James C. Estill                                  |
# CONTACT: JamesEstill_@_gmail.com                          |
# STARTED: 07/27/2007                                       |
# UPDATED: 07/27/2007                                       |
#                                                           |
# DESCRIPTION:                                              |
#  Merge a two directories containing subdirs with the same |
#  name into a single dir. One dir can serve as the template|
#  or an entirely new directory set could be created.       |
#                                                           |
# VERSION:                                                  |
# $Id::                                                  $: |
#                                                           |
#                                                           |
# LICENSE:                                                  |
#  GNU General Public License, Version 3                    |
#  http://www.gnu.org/licenses/gpl.html                     |  
#                                                           |
#-----------------------------------------------------------+

=head1 NAME

dir_merge.pl - Merge directories

=head1 VERSION

This documentation refers to program version 0.1

=head1 SYNOPSIS

 Usage:
  dir_merge.pl -i 'dir01,dir02' -o new_dir

=head1 DESCRIPTION

This is what the program does

=head1 REQUIRED ARGUMENTS

=over 2

=item -i,--indir

The directories to merge.

=item -o,--outdir

The directory the merged dirs will be moved to.

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

=head1 DIAGNOSTICS

The list of error messages that can be generated,
explanation of the problem
one or more causes
suggested remedies
list exit status associated with each error

=head1 CONFIGURATION AND ENVIRONMENT

Names and locations of config files
environmental variables
or properties that can be set.

=head1 DEPENDENCIES

Other modules or software that the program is dependent on.

=head1 BUGS AND LIMITATIONS

Any known bugs and limitations will be listed here.

=head1 LICENSE

GNU LESSER GENERAL PUBLIC LICENSE

http://www.gnu.org/licenses/lgpl.html

=head1 AUTHOR

James C. Estill E<lt>JamesEstill at gmail.comE<gt>

=cut

package DAWGPAWS;

#-----------------------------+
# INCLUDES                    |
#-----------------------------+
use strict;
use Getopt::Long;
use File::Copy;

#-----------------------------+
# PROGRAM VARIABLES           |
#-----------------------------+
my $ver = "0.1";

#-----------------------------+
# VARIABLE SCOPE              |
#-----------------------------+
my $indir_list;                # List of dirs passed from the cmd line
my $outdir;                    # Directory to serve as parent for output
my @indirs;                    # List of directories
my $num_indirs;                # Total number of indirs

# Booleans
my $do_overwrite = 0;          # Overwrite existing files
my $help = 0;
my $quiet = 0;
my $show_help = 0;
my $show_usage = 0;
my $show_man = 0;
my $show_version = 0;
my $err = 0;                  # Boolean to indicate program encountered error
my $verbose = 0;

#-----------------------------+
# COMMAND LINE OPTIONS        |
#-----------------------------+
my $ok = GetOptions("i|indir=s"   => \$indir_list,
                    "o|outdir=s"  => \$outdir,
		    # Booleans
		    "overwrite"   => \$do_overwrite,
		    "usage"       => \$show_usage,
		    "version"     => \$show_version,
		    "man"         => \$show_man,
		    "h|help"      => \$show_help,
		    "$verbose"    => \$verbose,
		    "q|quiet"     => \$quiet,);


#-----------------------------+
# SHOW REQUESTED HELP         |
#-----------------------------+
if ($show_usage) {
    print_help("");
}

if ($show_help || (!$ok) ) {
    print_help("full");
}

if ($show_version) {
    print "\n$0:\nVersion: $ver\n\n";
    exit;
}

if ($show_man) {
    # User perldoc to generate the man documentation.
    system("perldoc $0");
    exit($ok ? 0 : 2);
}

# Exit if indir_list and outdir not specified at cmdline
if ( (!$indir_list) || (!$outdir) ) {
    print "\a";
    print "ERROR: An input list of dirs and and output dir is required.\n";
    print "For help use:\n$0 --help\n";
    exit;
}

#-----------------------------------------------------------+
# MAIN PROGRAM BODY                                         |
#-----------------------------------------------------------+

#-----------------------------+
# GET THE INPUT DIRS TO MERGE |
#-----------------------------+
@indirs = split(/\,/, $indir_list);

# List dirs to split
$num_indirs = @indirs;

# Got rid of this
# Throw error if multiple dirs not indicated
#if ($num_indirs < 2) {
#    print "\a";
#    print "ERROR: Multiple dirs shold be indicated by --indir\n";
#    print ""; 
#}

print "\nDirs to merge ($num_indirs):\n";
for (my $i=0; $i < $num_indirs; $i++) {

    # ADD SLASH TO END OF DIR PATH
    unless ($indirs[$i] =~ /\/$/ ) {
	$indirs[$i] = $indirs[$i]."/";
    }

    # CHECK FOR EXISTENCE OF INPUT DIRECTORIES
    if (-e $indirs[$i]) {
	print "\t-".$indirs[$i]."\n";
    }
    else {
	$err = 1;
	print "\a";
	print "ERROR: Directory does not exist\n\t".$indirs[$i]."\n";
    }
}

#-----------------------------+
# GET THE OUTPUT DIR          |
#-----------------------------+
unless ($outdir =~ /\/$/ ) {
    $outdir = $outdir."/";
}
print "New dir: $outdir\n\n";

# CREATE OUTDIR IF IT DOES NOT EXIST
unless (-e $outdir) {
    print "Creating output dir ...\n" unless $quiet;
    mkdir $outdir ||
	die "Could not create the output directory:\n$outdir";
}


#-----------------------------+
# FOR EACH DIR TO COPY VIA THE
# merge SUBFUNCTION
#-----------------------------+
for my $ind_dir (@indirs) {
    merge ($ind_dir, $outdir);
}

exit;

#-----------------------------------------------------------+ 
# SUBFUNCTIONS                                              |
#-----------------------------------------------------------+

sub print_help {

    # Print requested help or exit.
    # Options are to just print the full 
    my ($opt) = @_;

    my $usage = "USAGE:\n". 
	"dir_merge.pl -i DirsToMerge -o OutDir";
    my $args = "REQUIRED ARGUMENTS:\n".
	"  --indir        # Path to the directories to merge\n".
	"                 # Dirs separated by comma\n".
	"  --outdir       # Path to the output directory\n".
	"\n".
	"OPTIONS::\n".
	"  --version      # Show the program version\n".     
	"  --usage        # Show program usage\n".
	"  --help         # Show this help message\n".
	"  --man          # Open full program manual\n".
	"  --quiet        # Run program with minimal output\n";
	
    if ($opt =~ "full") {
	print "\n$usage\n\n";
	print "$args\n\n";
    }
    else {
	print "\n$usage\n\n";
    }
    
    exit;
}


sub merge {
# SUBFUCNTION SRC:
# http://software.hixie.ch/utilities/unix/merge/merge.pl

    my($a, $b) = @_;
    #   SOURCE: not -e    -f    -d   empty -d else
    # TARGET:
    #  not -e    err      mv    mv    mv      err
    #
    #      -f    err      err   err   err     err
    #
    #      -d    err      err   loop  rmdir   err
    #
    #    else    err      err   err   err     err
    print "merge: merging $a and $b\n";
    if (not -e $a) {
	#-----------------------------+
	# $a DOES NOT EXIST           |
	#-----------------------------+
        print STDERR "merge:$a: doesn't exist\n";

    } 
    elsif (not (-f $a or -d $a)) {
	#-----------------------------+
	# $a IS NOT A FILE OR DIR     |
	#-----------------------------+
        print STDERR "merge:$a: not a normal file\n";

    } 
    elsif (not -e $b) {
	#-----------------------------+
	# $b DOES NOT EXIST           |
	#-----------------------------+
        print "merge: moving $a to $b\n" if $verbose;
        rename($a, $b) || 
	    print STDERR "merge:$a: could not rename to $b, $!\n";;
	# JAMIE MOD BELOW -- curently not working
        #print "dir_merge: moving $a to $b\n" if $verbose;
	#copy ($a, $b) ||
	#    print STDERR "merge:$a: could not rename to $b, $!\n";
    } 
    elsif (-d $b) {
	#-----------------------------+
	# $b IS A DIRECTORY           |
	#-----------------------------+
        if (-d $a) {
            my @entries = getdir($a);
            if (@entries) {
                # not empty
                # recurse through it to give us a chance to make it empty
                print "merge: going through contents of $a\n" if $verbose;
                foreach my $entry (@entries) {
                    my $c = "$a/$entry";
                    $c =~ s|//|/|gos;
                    my $d = "$b/$entry";
                    $d =~ s|//|/|gos;
                    &merge($c, $d);
                }
            }
            # empty now?
            @entries = getdir($a);
            if (not @entries) {
                print "merge: deleting empty directory $a\n" if $verbose;
                rmdir($a) or print STDERR "merge:$a: could not delete ".
		    "directory, $!\n";
            } 
	    else {
                print STDERR "merge:$a: could not delete directory,".
		    " directory is not empty\n";
            }
        } 
	else {
            print STDERR "merge:$a: conflicts with directory $b\n";
        }
    } 
    else {
	#-----------------------------+
	# THE FILE ALREADY EXIST WITH |
	# A FILE IN B                 |
	#-----------------------------+
	print  "merge:$a: conflicts with non-directory $b\n";
	
	if ($do_overwrite) {
	    # DELETE EXISTING FILE AND REPLACE WITH NEW FILE
	    #print STDERR "File $b exists and will be overwritten\n";
	    #my $rmcmd = "rm $b";
	    #system ($rmcmd);
	    #rename($a, $b) || 
	    #	print STDERR "merge:$a: could not rename to $b, $!\n";
	    move($a, $b) || 
		print STDERR "\n\nmerge:$a: could not rename to $b, $!\n";
	}
	else {
	    print STDERR "merge:$a: conflicts with non-directory $b\n";
	}
	
    }
}


sub getdir {
# SUBFUCNTION SRC:
# http://software.hixie.ch/utilities/unix/merge/merge.pl
    my($a) = @_;
    local *DIR;

    unless (opendir(DIR, $a)) {
        print STDERR "dir_merge:$a: can't open directory\n";
        return;
    }

    my @entries;
    while (my $entry = readdir(DIR)) {
        if ($entry !~ m/^\.\.?$/o) {
            push(@entries, $entry);
        }
    }
    
    closedir(DIR) || 
	print STDERR "dir_merge:$a: could not close directory,$!\n";
    
    return @entries;
}

=head1 HISTORY

STARTED: 07/27/2007

UPDATED: 07/27/2007

VERSION: $Id:$

=cut

#-----------------------------------------------------------+
# HISTORY                                                   |
#-----------------------------------------------------------+
#
# 07/27/2007
# - Program started
# - Added getdir from external source
