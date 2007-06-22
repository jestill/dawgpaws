#!/usr/bin/perl -w
#-----------------------------------------------------------+
#                                                           |
# ParseFGenesh                                              |
#                                                           |
#-----------------------------------------------------------+
#                                                           |
# UPDATED: 11/28/2006                                       |
# Jamie modified this from something that Renyi had.        |
#                                                           |
# A parser of the Fgenesh contig-by-contig results;

use strict;

my $insertDB = 1;
my $printScreen = 1;
my $debug = 0;
#my $sql = SQL->new('root','riofuy','briggs_elegans_synteny');
#my $output_dir = 'output/';
my $usage = "$0 [<fgenesh output dir>]\n";
my $inputDir = "/home/renyi/projects/maize_annot/fgenesh_txt";
if(@ARGV > 0){
	$inputDir = $ARGV[0];
}
chdir $inputDir;
opendir (DIR, $inputDir) or die "Cannot open $inputDir: $!";
my @files = grep /fasta/, readdir DIR;

my $dbh = getConnection("maize_annot", "renyi", "suizhou");
die "Cannot connect to database.\n" if (!defined $dbh);
my $bac_sth = $dbh->prepare("INSERT INTO bac_info_fgenesh (file, clone, length, ".
                         "total_genes, total_genes_plus, total_genes_minus, ".
						 "total_exons, total_exons_plus, total_exons_minus) VALUES ".
						 "(?, ?, ?, ?, ?, ?, ?, ?, ?)");
my $protein_sth = $dbh->prepare("INSERT INTO proteins_fgenesh (file, clone, fgenesh_id, ".
                       "method, strand, exons, plength, start, stop, sequence".
					   ") VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)");

my $feature_sth = $dbh->prepare("INSERT INTO features_fgenesh (predicted_id, clone, fgenesh_id, ".
                         "strand, exon_id, feature, start, stop, score, orf_start, orf_stop, ".
						 "flength) VALUES (?,?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)");

my $counter = 0;
foreach  my $file (@files) {
  $counter++;
  print STDERR "$counter\t$file\n";
  my $features = {};
  my $proteins = {};
  my $d = {};
  $d->{parsing_predicted} = 0;
  $d->{parsing_protein} = 0;
  $d->{parsing_mRNA} = 0;
  # Open the file for reading
  open (IN,"$file") or die "Cannot open $file: $!";
  while (my $line = <IN>) {
    chomp $line;
    next if ($line eq ''); # skip the empties
    $d->{clone}     = ($line =~ /Seq\sname:\s(\S+)\s(.*)/) ? $1 : $d->{clone};
    
    if ($line =~ /Length\sof\ssequence:\s(\d+)/) {
      $d->{length} = $1;
      #$d->{gc}     = 0;
    }
    
    if ($line=~/.*predicted\sgenes\s(\d+)\sin\s\+chain\s(\d+)\sin\s-chain\s(\d+).*/) {
		if($debug){
		#	print STDERR "Feature line: $line\n";
		}
      $d->{total_genes}  = $1;
      $d->{total_genes_plus}  = $2;
      $d->{total_genes_minus} = $3;
    }
    
    if ($line=~/.*predicted\sexons\s(\d+)\sin\s\+chain\s(\d+)\sin\s-chain\s(\d+).*/) {
		if($debug){
			print STDERR "Feature line: $line\n";
		}
      $d->{total_exons}  = $1;
      $d->{total_exons_plus}  = $2;
      $d->{total_exons_minus} = $3;
    }
    
    # Have we already encoutnered the header line that signals
    # start of the data?
    if ($d->{flag}) {
      # Parse each line - store multiple information together by gene
      # Store each predicted protein in a hash of arrays of arrays
      
      if ($line =~ /\s*(\d+)\s+([+-])\s+(\d+)\s+(\w+)\s+(\d+)\s+-\s+(\d+)\s+(\d+\.\d+)\s+(\d+)\s+-\s+(\d+)\s+(\d+).*/) {
     if($debug){
		# print STDERR "Feature list: $line\n";
	 }
	# This parses correctly for the most part
	my $id = $1;
	# Each feature for this gene should be in a hash
	# order: clone,strand,exon_id, feature,start,end,score,orf_start,orf-end,length)
	push (@{$features->{$id}},[ $d->{clone},$2,$3,$4,$5,$6,$7,$8,$9, $10 ]);
      } elsif ($line =~ /\s*(\d+)\s+([\+-])\s+(\w+)\s+(\d+)\s+(\d+\.\d+)/) {
		  if($debug){
			  print STDERR "Feature list: $line\n";
		  }
	# some records look like this
	# 1 -     PolA     617                1.05  
	my $id = $1;
	# Each feature for this gene should be in a hash
	# order: clone,strand,exon_id, feature,start,end,score,orf_start,orf-end,length)
	push (@{$features->{$id}},[ $d->{clone},$2,'', $3,$4,'',$5,'','','' ]);

      } else {}
    }

    # Am I parsing the predcited proteins now?
    if ($d->{parsing_predicted}) {
      # first line is the identifier line
	  if($line =~ /.*FGENESH:\s+(\d+)\s+(\d+)\s+exon\s+\(s\)\s+(\d+)\s+-\s+(\d+)\s+(\d+).*chain\s+([\+-])/){
		  if($debug){
			  print STDERR "protein line $line\n";
		  }
      $d->{parsing_protein} = 1;
	  $d->{parsing_mRNA} = 0;
	# Parse out each of the fields
	# method id exons position aa strand
	# still need to account for positions, aa, and strand
	# >FGENESH   1   8 exon (s)    203  -   3137    601 aa, chain +
	
	# This regex isn't working
	#$line =~ />FGENESH\s+(\d+)\s+(\d+)\sexon.*(\d+)\s+-\s+(\d+)\s+(\d+)\saa,\s+chain\s([+-])/;
	
	my $id = $1;
	$proteins->{$id}->{clone}      = $d->{clone};
	$proteins->{$id}->{method}     = 'FGENESH';
	$proteins->{$id}->{fgenesh_id} = $1;
	$proteins->{$id}->{exons}      = $2;
	$proteins->{$id}->{start}      = $3;
	$proteins->{$id}->{stop}       = $4;
	$proteins->{$id}->{plength}    = $5;
	$proteins->{$id}->{strand}     = $6;
	
	# Need to store the ID for later - storing sequences
	$d->{current_id} = $id;
      
   } 
	  
	  if($d->{parsing_protein} && (!$d->{parsing_mRNA}) && ($line !~ /FGENESH/)){  # These are sequence lines then
	# keep concatenating sequence lines
	      $proteins->{$d->{current_id}}->{sequence} .= $line;
	  }
    

	if($line =~ /FGENESH:\s*\[mRNA\]/){
		$d->{parsing_protein} = 0;
		$d->{parsing_mRNA} = 1;
	}
  }
    
    # remaining lines are the sequence
    # Am I done with the data set?
    if ($line =~ /Predicted\s+protein/) {
      $d->{parsing_predicted} = 1;
      delete $d->{flag};
    }
    
    # Extracting the salient information is going to be a pain
    # Need to find out what exactly all of the fields are.
    if ($line =~ /G\s+Str\s+Feature\s+/) {
      # We've arrived at the features.  Set a flag.
      $d->{flag}++;
	  if($debug){
		  print STDERR "Entering feature list section.\n";
	  }
    }
  } # next line
  
  # Now let's print everything:
  # First, the easy stuff.  Clone info.  


  #delete $d->{parsing_predicted};
  #delete $d->{current_id};
  #delete $d->{parsing_protein};
  #delete $d->{parsing_mRNA};
  
  # Insert the contig info
 my @bac_info = ($file, $d->{clone}, $d->{length}, $d->{total_genes}, 
                 $d->{total_genes_plus}, $d->{total_genes_minus}, $d->{total_exons},
				 $d->{total_exons_plus}, $d->{total_exons_minus});
 if($printScreen){
	 print "Insert bac_info_fgenesh: @bac_info\n";
 }
 if($insertDB){
	 $bac_sth->execute(@bac_info);
 }
 
  foreach my $id (sort keys %$proteins) {
    my $ref = $proteins->{$id};
    #foreach my $f (keys %{$proteins->{$id}}) {
    #  $vals{$f} = $proteins->{$id}->{$f};
    #}
    
    # Do the insert
#    $sql->insert(-table=>'predicted_proteins',
#		 -field_vals=>\%vals);
    my @vals = ($file, $ref->{clone}, $ref->{fgenesh_id}, $ref->{method},
	            $ref->{strand}, $ref->{exons}, $ref->{plength}, $ref->{start},
				$ref->{stop}, $ref->{sequence});
	if($printScreen){
		print "Insert protein_fgenesh @vals\n";
	}
	my $pred_id = 0;
	if($insertDB){
		$protein_sth->execute(@vals);
		$pred_id = $dbh->{'mysql_insertid'};
	}
    
    # Processing features
    #my @cols = (qw/clone strand feature start stop 
	#	score orfstart orfstop flength/ );
    foreach my $array (@{$features->{$id}}) {
      #my $count = 0;
      #my %featvals;
      #$featvals{predicted_id} = $pred_id;
      #$featvals{fgenesh_id} = $id;
      #foreach my $col (@cols) {
	#$featvals{$col} = $array->[$count];
	#$count++;
     # }
      # Do the insert
      #$sql->insert(-table=>'features',
	#	   -field_vals=>\%featvals);
	  my @fvals = @$array;
	  my $clone = shift @fvals;
	  unshift @fvals, ($pred_id, $clone, $id);
	  if($printScreen){
		  print "Insert features_fgenesh @fvals\n";
	  }
	  if($insertDB){
		  $feature_sth->execute(@fvals);
	  }
    }
  }
  if($debug && $counter >= 1){
	  die "That's it for test!\n";
  }
} # next file

# Three tables:
# General clone info
# Predicted proteins
# unique id features
