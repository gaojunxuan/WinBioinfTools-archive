#!/usr/bin/perl
#create infofile
# 1st line: number of genomes (k)
# 2nd to 2+kth line: length of genomes

#use strict;



$debug = 1;                  # if you're inclined to have debug levels...
$| = 1;                      # flush output to stdout immediately

my $argc = @ARGV;
if ($argc < 1)
   {
     print "usage: \n perl create_info_file.pl memspex_matchfile \n";
     exit -1;
   }
my $infilename = $ARGV[0];
print "starting create_info_file.pl $infilename\n";

my $dir_name = qx/dirname $infilename/;
chop $dir_name;
print "dirname is: $dir_name \n";



open(matchfile,"<".$infilename)             # filehandles don't use $name...
    or die "Oops! can't open input file ".$infilename."\n";


$number = 0; #number of genomes
@genome_length;

while ($infileline = <matchfile>){           # retrieve file, line by line


if($infileline=~/\#/){# Kommentarzeile
      #Laenge der 3 Genome rausparsen:
	if ($infileline=~/\#\s+file=.+\s+(\d+)\s+\d+/){
		my $length = $1;
		print "create_info_file: length = $length\n";
		push(@genome_length, $length);
	}

	elsif($infileline=~/width(\d+)/){
		$number = $1 + 1; #number of genomes

	}
	elsif($infileline=~/queryfile=.+length\s+(\d+)/){
		my $length = $1;
		print "create_info_file: length = $length\n";
		push(@genome_length, $length);
	}
	elsif($infileline=~/output the matches/){
		last; #matches beginnen, fertig mit kommentarzeilen
	}
	
}
}
close(matchfile);

$infofile=$infilename.".info";
#1. Zeile: Anzahl der Genome (k); dann k Zeilen mit der LÃ¤nge des 1. bis k-ten genoms; dann alle outputfiles (2^(k.1))
open(infof,">$infofile")
  or die "failed to open $infofile \n";

print "printing information to infofile $infofile : number = $number\n";
print infof "$number\n";

# Vorsicht: im memspex output ist genomlänge nicht enthalten. wenn sie später gebraucht wird,
# anders besorgen!
foreach (@genome_length){
	print infof "$_\n";
}
#foreach (@outfilenames){
#	print infof "$_\n";
#}

close(infof);

