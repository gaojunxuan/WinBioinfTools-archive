use strict;

my $filename = shift(@ARGV);
open(IN,"<$filename");
open(OUT,">$filename.sorted");

my @lines, my @sorted;
my $line,my $index = 0;
while(chomp($line = <IN>)){
    @lines[$index] = $line;
    $index++;
}

@sorted = sort compare_kd @lines;

sub compare_kd {
    my $x, my $y;
    $a =~ /\d+\s(\d+)\s/;
    $x = $1;
    $b =~ /\d+\s(\d+)\s/;
    $y = $1;
    if    ($x == $y) { return  0 }
    elsif ($x >  $y) { return -1 }
    else             { return  1 }
}

select(OUT);
while($index > 0)
{
    $index--;
    print @sorted[$index];
    print "\n";
}
