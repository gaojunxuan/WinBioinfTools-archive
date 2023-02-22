#!/usr/bin/perl
#

sub dec2bin{
#decimal zahl wird zu binärzahl mit führenden nullen konvertiert
#usage $bin = dec2bin($dec, $size)
	my $dec = $_[0];
	my $size = $_[1];
	my $bin = "";
	if ($dec == 0) {
		$bin = "0".$bin;
	}
	else{
		while ($dec > 0){
			my $mod = $dec % 2; #Rest
			$bin = $mod.$bin;
			$dec = $dec / 2;
			if($mod == 1){$dec = $dec -0.5}
		}
	}
	#fuehrende nullen anhangen, falls $bin weniger stellen hat als $size
	while (length($bin)< $size){
		$bin = "0".$bin;
	}
	$bin;
}


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
    print "usage: \n perl generate_suffixes_k_genomes.pl k [outfile] \n";
    exit -1;
}


my $number_of_suffixes = 2**($k-1); # matches in the first genome are always on the leading strand

for ($i = 0; $i < $number_of_suffixes; $i++)
{
	my $bin = dec2bin($i, $k-1);
	# 0 durch - ersetzen; 1 durch +
	$bin  =~ s/0/m/g;
	$bin =~ s/1/p/g;
	#print bin
	if($argc == 1){
	    print "$bin\t";
	    
	}
	if($argc == 2){
	    open(outfile,">>".$outfilename)
		or die "Oops! can't open output file ".$outfilename."\n";
	    print outfile "$bin\t";	    
	    close(outfile);
	}
}

if($argc == 1){
    print "\n";
    
}
if($argc == 2){
    open(outfile,">>".$outfilename)
	or die "Oops! can't open output file ".$outfilename."\n";
    print outfile "\n";	    
    close(outfile);
}
