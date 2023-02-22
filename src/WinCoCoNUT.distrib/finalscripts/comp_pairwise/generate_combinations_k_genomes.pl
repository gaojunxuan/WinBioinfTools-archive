#!/usr/bin/perl
#

$debug = 1;                  # if you're inclined to have debug levels...
$| = 1;                      # flush output to stdout immediately

my $argc = @ARGV;
my $outfilename;
my $k;
if ($argc == 1)
   {
	 $k = $ARGV[0];
	 #print to stdout
	 }
if ($argc == 2)
   {
	 $k = $ARGV[0];
	 $outfilename = $ARGV[1];
	 open(outfile,">".$outfilename)
			or die "Oops! can't open output file ".$outfilename."\n";
		close(outfile);
	 }
if (($argc > 2) || ($argc == 0))
	{
	print "usage: \n perl generate_combinations_k_genomes.pl k [outfile] \n";
	exit -1;
	}
	
	
for ($i = 1; $i < $k; $i++)
{
	for ($j = $i+1; $j <= $k ; $j++ ){
		my $combination = "$i"."x$j";
		if($argc == 1)
			{print "$combination\t";}
		if($argc == 2)
		{
			open(outfile,">>".$outfilename)
			or die "Oops! can't open output file ".$outfilename."\n";
			print outfile "$combination\t";
			close(outfile);
		}
	}
}
