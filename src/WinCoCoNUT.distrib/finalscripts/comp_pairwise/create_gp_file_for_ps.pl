#!/usr/bin/perl


#perl /home/hockel/bin/create_gp_file.pl ${FRAGMENTFILE}.dat ${FRAGMENTFILE}.gp

$argc = @ARGV;

if ($argc != 2)
   {
     print "usage: \n perl create_gp_file_for_ps.pl infile.dat outfile.gp\n";
   }

my $datfile = $ARGV[0];
my $gpfile = $ARGV[1];

open(OUT, ">$gpfile");

print OUT "set terminal postscript landscape color\n";
print OUT "set out '$gpfile.ps'\n";
print OUT "set nokey\n";
print OUT "set nogrid\n";
print OUT "set xlabel \"\"\n";
print OUT "set ylabel \"\"\n";
print OUT "set xrange [0:*]\n";
print OUT "set yrange [0:*]\n";
print OUT "plot \"$datfile\" index 0:3 title \"Fwd\" w lp lt 1 lw 2 pt 1 ps .75, \"$datfile\" index 4:8 title \"Rev\" w lp lt 2 lw 2 pt 1 ps .75\n";

close(OUT);
