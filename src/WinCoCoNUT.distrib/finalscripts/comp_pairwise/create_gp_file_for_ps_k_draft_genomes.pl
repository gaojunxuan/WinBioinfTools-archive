#!/usr/bin/perl


$argc = @ARGV;

if ($argc < 4){
   print "usage: \n create_gp_file_for_ps_k_draft_genomes.pl infile.dat number_of_genomes projection desfile outfile.gp forward_flag xlabel ylabel \n";
   exit;
}

my $datfile = $ARGV[0];
my $desfile = $ARGV[3];
my $gpfile = $ARGV[4];
my $forward_flag= $ARGV[5];
my $x_label= $ARGV[6];
my $y_label= $ARGV[7];
my $k=$ARGV[1];
my $projection= $ARGV[2];
my $num_suffixes=2**($k-1);
my $split_index = 2**($k-2);
my $split_index_minus_1 =$split_index-1;

open(OUT, ">$gpfile");

print OUT "set terminal postscript landscape color\n";
print OUT "set out '$gpfile.ps'\n";
print OUT "set nokey\n";
print OUT "set nogrid\n";
print OUT "set xlabel \"$x_label\"\n";
print OUT "set ylabel \"$y_label\"\n";

#print OUT "set xrange [0:*]\n";
#print OUT "set yrange [0:*]\n";


open(DESFILE, "<$desfile");

$draft_genome=1;
@projection_array=split(/x/,$projection);
my $x_axis=int($projection_array[0]);
#print "X axis ".$x_axis."\n";
my $y_axis=int($projection_array[1]);
$counter=1;
$sep=0;
$xsep=0;
$ysep=0;
while($desfilelines=<DESFILE>){

     if($desfilelines=~/\#---/){
	 $sep=0;
	 $draft_genome++;	
	 if($open_falg==1){
	     print OUT ")\n";
	     $open_falg=0;
	 }
	$start_flag=0;
	 next;
     }
     else{
	 if($start_flag==1){
	     print OUT ", ";
	 }
     }
    
     @array = split(/\ /,$desfilelines);
     $sep=$array[1]+$sep+1;
     if($draft_genome==$x_axis){
	 if($start_flag==0){
	     $start_flag=1;
	     print OUT "set x2tics (";
	     $open_falg=1;
	 }
	 print OUT " \" \" $sep";
	 #print OUT "set arrow $counter from $sep, graph 0 to $sep, graph 1  nohead  lt 9 lw 1 \n";
	 #print OUT "set arrow $counter from $sep, graph 0 to $sep, graph 1  nohead  lt 3  \n";
	 #print OUT "set arrow from $sep, graph 0 to $sep, graph 1  head \n";
	 $counter++;
	 $xsep=$sep;
     }

     if($draft_genome==$y_axis){
	 if($start_flag==0){
	     $start_flag=1;
	     print OUT "set y2tics (";
	     $open_falg=1;
	 }
	 print OUT " \" \" $sep";
	 #print OUT "set arrow $counter from $sep, graph 0 to $sep, graph 1  nohead  lt 9 lw 1 \n";
	 #print OUT "set arrow $counter from $sep, graph 0 to $sep, graph 1  nohead  lt 3  \n";
	 #print OUT "set arrow from $sep, graph 0 to $sep, graph 1  head \n";
	 $counter++;
	 $ysep=$sep;
     }
     
}
if($open_falg==1){
    print OUT ")\n";
    $open_falg=0;
}

print OUT "set xrange [0:$xsep]\n";
print OUT "set yrange [0:$ysep]\n";

print OUT "\n set grid  noxtics noytics  x2tics y2tics \n";
close(DESFILE);

my $commanc;
my @suffixes;

if(($argc > 4) and ($forward_flag==1)){
    $result="p";
    for($x=2;$x<$k;$x++){
	$result=$result."p";
    }
    @suffixes = split(/\s/, $result);
}
else{
#compute suffixes and combinations for k genomes:
    $command = "generate_suffixes_k_genomes_multimat.pl ".$k ;
    $result = qx($command);
#print $result;
    @suffixes = split(/\s/, $result);
}
$x_axis=int(substr($projection,0,1));
#print "X axis ".$x_axis."\n";
$y_axis=int(substr($projection,1,2));
#print "Y axis ".$y_axis."\n";
my $x_axis_direction;
my $y_axis_direction;

#@suffixes = (pppp,pppm,ppmp,ppmm,pmpp,pmpm,pmmp,pmmm);
$number = @suffixes;
print OUT "plot ";
for($i = 0; $i < $number; $i++){
	$suffix = "p".$suffixes[$i];
	#print $suffix." \n";
	$x_axis_direction=substr($suffix,($x_axis-1),1);
	$y_axis_direction=substr($suffix,($y_axis-1),1);
	#print $x_axis_direction." xdir\n"; 
	#print $y_axis_direction." ydir\n"; 
	#print substr($suffix, 0,5)." NN\n";
	#$color = $i+1;
	if(($x_axis_direction eq "p")and ($y_axis_direction eq "p") ){
	    $color = 1;
	}
	elsif(($x_axis_direction eq "p")and ($y_axis_direction eq "m") ){
	    $color= 2;
	}
	elsif(($x_axis_direction eq "m")and ($y_axis_direction eq "p") ){
	    $color= 2;
	}
	elsif(($x_axis_direction eq "m")and ($y_axis_direction eq "m") ){
	    $color= 1;
	}
	print OUT "\"$datfile\" index $i title \"$suffix\" w lp lt $color lw 2 pt 1 ps .75";
	if ($i < $number-1) {
		print OUT ", ";
	}
	else{ #letzter eintrag
		print OUT " \n";
	}

}

#print OUT "plot \"$datfile\" index 0:7:1 w lp lt 1 lw 2 pt 1 ps .75\n";

#old one for projection 12
#print OUT "plot \"$datfile\" index 0:$split_index_minus_1 w lp lt 2 lw 2 pt 1 ps .75, ";
#print OUT "\"$datfile\" index $split_index:$num_suffixes w lp lt 1 lw 2 pt 1 ps .75 \n";


close(OUT);

