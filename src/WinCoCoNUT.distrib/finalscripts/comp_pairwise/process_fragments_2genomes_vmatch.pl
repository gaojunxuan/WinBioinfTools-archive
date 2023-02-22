#!/usr/bin/perl
#

#$debug = 1;                  # if you're inclined to have debug levels...
#$| = 1;                      # flush output to stdout immediately

my $argc = @ARGV;

my $fragmentfile;
my $k; #number of involved sequences / genomes
my $verbose_flag=0;
my $plot_flag=0;
$draft_flag=0;

if ($argc < 4)
{
    print "Usage: \n"."process_fragments_kgenomes_vmatch.pl fragment_file number_of_genomes chaining_prog verbose_flag parameter_file plot_flag forwardflag usematch_flag usechain_flag in_seq_string\n";
    exit -1;
}


$fragmentfile = $ARGV[0];
$k = $ARGV[1];
#$chainer = $ARGV[2];
#$chainer ="C:/cygwin/home/H/MicroSoftCoCoNUT/bin/chainer/chainer.exe";
$chainer ="Home/MicroSoftCoCoNUT/bin/chainer/chainer";

$verbose_flag=$ARGV[3];
$plot_flag=$ARGV[5];
$forward_flag=$ARGV[6];
$usematch_flag=$ARGV[7];
$usechain_flag=$ARGV[8];

@input_files;
if($plot_flag==1){
    for(my $i=9; $i<$argc; $i++){ 	    		
	push(@input_files,$ARGV[$i]);	
    }
     
}
#print "\nHIIIIIIIIIIIII: ".$ARGV[2]."\n";
#print "HIIIIIIIIIIIII: ".$ARGV[9]."\n";

#change directories and parameters according to your needs:
#bindir: directory, where create_info_file.pl, generate_suffixes_k_genomes.pl, 
# generate_combinations_k_genomes.pl, memspex2chainer.pl and ccn2gnuplot.pl can be found
# (with a slash at the end)
#$bindir = "./";
#chainer: chainer executable


# default parameters for first and second (recursivly) time chaining
#my $chainerparameter1 = " -l -chainerformat -s 0 ";
#my $chainerparameter2 = " -l -chainerformat -gc 10000 -s 50000 -lw 100 ";


#fragmentfile in directory und filename zerlegen
my $dirname = qx/dirname $fragmentfile/;
chop $dirname;
if($verbose_flag==1){
    print "process_fragments_2genomes_vmatch.pl: dirname is: $dirname \n";
}
my $basename = qx/basename $fragmentfile/;
chop $basename;
if($verbose_flag==1){
    print "process_fragments_kgenomes_multimat.pl: basename = $basename \n";
}


open(fragfile,"<".$fragmentfile."+")             # filehandles don't use $name...
	or die "Oops! can't open fragment file ".$fragmentfile."\n";    
close(fragfile);


open(fragfile,"<".$fragmentfile."-")             # filehandles don't use $name...
	or die "Oops! can't open fragment file ".$fragmentfile."\n";    
close(fragfile);



my $dir_match2chaainer = qx/dirname $chainer/;
chop $dir_match2chaainer;
$dir_match2chaainer="home/Microsoftcoconut/bin/chainer/";


# LOADING CHAINING PARAMETERS IF A CHAIN PARAMETER IS GIVEN
my @chaining_param;    
my $chaining_number=0;


if($argc > 3 ){   
    $param_file=$ARGV[4];
    my $infileline;
    open(param_file,"<".$param_file)             # filehandles don't use $name...
	or die "Oops! can't open parameter file \n";     
    
    while ($infileline = <param_file>)           # retrieve file, line by line
    {   #chop $infileline;
	#print $infileline;
	# check for commented lines starting with #
	$tmp_infileline=$infileline;
	@tmp_tokens=split(" ",$tmp_infileline);
	@tmp_tokens=split('',$tmp_tokens[0]);
	if($tmp_tokens[0] eq "#"){next;}

	#print $infileline;
	
	if($infileline=~/CHAINING=(.+)/){ # get path without endline
	    $chaining_param[$chaining_number]=$1;
	    $chaining_number++;

	}

    }
    close(param_file);
}
if($chaining_number==0){
 #   $chaining_param[0]=$chainerparameter1;
 #   $chaining_param[1]=$chainerparameter2;
 #   $chaining_number=2;    
    print "Failuer: Chaining parameters are not given in parameter file\n";
    exit;
}

my $command;
my $argstring;


# remove previous files
#$command ="rm  ".$fragmentfile.".*p*";
#$result = qx($command);

#create info file
#if (!(-e $fragmentfile.".info")){
#	$command = $myperl.$bindir."create_info_file.pl ".$fragmentfile;
#       print "command = $command\n";
#	qx($command);
#}


$contig_num1=get_number_of_contigs("$input_files[0].seqinfo");
push(@contig_num,$contig_num1);
$contig_num=get_number_of_contigs("$input_files[1].seqinfo");
push(@contig_num,$contig_num1);
print "Contig numbers: ".$contig_num[0]." ".$contig_num[1]."\n";

    
if(($contig_num[0]>1) or ($contig_num[0]>1)){
    #modify the parameters
    $draft_flag=1;
    #$chaining_number=1; # make only one chaining step
    for(my $x=0;$x<$chaining_number;$x++){
	$tmp_str=$chaining_param[$x];
	
	if($tmp_str=~/-l\s+/){	    
	    $tmp_str=~ s/-l\s+/ -draft /;	    
	    $chaining_param[$x]=$tmp_str;
	}
    }
        	   
}

for(my $x=0;$x<$chaining_number;$x++){
    if($verbose_flag==1){
	print "Parameters of chaining call no. ".$x." :".$chaining_param[$x]."\n";
    }
}

$forward_flag=0;
if(($forward_flag==0) or ($argc < 7) ){
#compute suffixes and combinations for k genomes:
    $command = "perl Home/MicroSoftCoCoNUT/finalscripts/comp_pairwise/generate_suffixes_k_genomes_multimat.pl ".$k ;
    $result = qx($command);
 #   print "HALOO22222OOOOOOOO: ".$result."\n";exit 0;
    @suffixes = split(/\s/, $result);

    
}
else{
    $result="p";
    for($x=2;$x<$k;$x++){
	$result=$result."p";
    }
    #@suffixes = split(/\s/, $result);
    $suffixes[0]=$result;
}

#$command = "generate_combinations_k_genomes.pl ".$k ;
#$result = qx($command);
$result ="1x2";
@combinations = split(/\s/,$result);

#format conversion from multimat to chainer input:
if($verbose_flag==1){
    print "Start fragment processing $fragmentfile \n";
}


if($usematch_flag==0){
    $argstring="home/Microsoftcoconut/bin/chainer/gt2chainer -i ".$fragmentfile."+"." -o ".$fragmentfile.".pp";

    if($verbose_flag==1){
	print "Running: $argstring\n";
    }

    my $retcode =system($argstring);
    if($retcode ne 0){
	print STDERR "failure: $argstring\n";
	exit 1;
    }
    # reverse complement
    $argstring=$dir_match2chaainer."/gt2chainer -i ".$fragmentfile."-"." -o ".$fragmentfile.".pm";

    if($verbose_flag==1){
	print "Running: $argstring\n";
    }

    my $retcode =system($argstring);
    if($retcode ne 0){
	print STDERR "failure: $argstring\n";
	exit 1;
    }
 
}
    


#chaining:
my $file_extension="ccn";
if ($usechain_flag==0)
{

if($verbose_flag==1){
    print "Chaining the first time\n";
}



foreach $suffix (@suffixes){

#    print "Suffixes: ".$suffix."\n";
     if(-e $fragmentfile.".p".$suffix){
	if($draft_flag){
	    $command = $chainer." ".$chaining_param[0]."  -des ".$fragmentfile.".seqinfo.p$suffix  ".$fragmentfile.".p".$suffix;
	}
	else{	
	    $command = $chainer." ".$chaining_param[0]."  ".$fragmentfile.".p".$suffix;
	}
	if($verbose_flag==1){
	    print "command = $command\n";
	}
        system($command);
    }
    
}



generate_info_file($input_files[0].".seqinfo",$input_files[1].".seqinfo",$fragmentfile);
#exit 0;

###############################################################
#DISPLAYING THE FIRST CHAIN 

if($plot_flag==1){
#prepare plot (format transformations etc):
if($verbose_flag==1){
    print "Running ccn2gnuplot_multimat.pl\n";
}


#cleaning files up:
$tmp_name=$fragmentfile.".".$file_extension.".dat";

if(-e $tmp_name){
    if($verbose_flag==1){
	print "Deleting".$tmp_name."\n";
    }
    unlink($tmp_name);    
}
foreach $combination (@combinations){
    if (-e $fragmentfile.".".$file_extension.".".$combination.".dat"){
	if($verbose_flag==1){
	    print "Removing: ".$fragmentfile.".".$file_extension.".".$combination.".dat"."\n";
	}
	unlink($fragmentfile.".".$file_extension.".".$combination.".dat");
    }
    foreach $suffix (@suffixes){
		
	if (-e $fragmentfile.".p".$suffix.".".$file_extension.".".$combination){
	    my $tmp_name=$fragmentfile.".p".$suffix.".".$file_extension.".".$combination;
	    if($verbose_flag==1){
		print  "Removing ".$tmp_name." \n";
	    }
	    unlink($fragmentfile.".p".$suffix.".".$file_extension.".".$combination);
	}
	if (-e ($fragmentfile.".p".$suffix.".".$file_extension.".dat")){	    
	    if($verbose_flag==1){
		print  "Removing ".$fragmentfile.".p".$suffix.".".$file_extension.".dat \n";
	    }
	    unlink($fragmentfile.".p".$suffix.".".$file_extension.".dat");
	}
    }
}
##########end cleaning up

if(-e $fragmentfile.".info"){
	$command = "perl Home/MicroSoftCoCoNUT/finalscripts/comp_pairwise/ccn2gnuplot_multimat.pl ".$fragmentfile.".info ".$file_extension;
	if($verbose_flag==1){
	    print "command = $command\n";

	}
	system($command);
}
else {
	print "Error: can't find infofile ".$fragmentfile.".info\n";
	exit;
}

#debug
#foreach $suffix (@suffixes){
#    print "HIIIIIIIIIIII $suffix"."\n";
#}




#if (! (-e $tmp_name)){
foreach $suffix (@suffixes){	    
    $command = "cat ".$fragmentfile.".p".$suffix.".".$file_extension.".dat >>".$fragmentfile.".".$file_extension.".dat";
    if($verbose_flag==1){
	print $command."\n";
    }
    system($command);
    #file is needed as input for chainer2permutation, which creates permutations from chainer output data
}
#}



foreach $combination (@combinations){
    if (-e $fragmentfile.".".$file_extension.".".$combination.".dat"){
	if($verbose_flag==1){
	    print "Removing: ".$fragmentfile.".".$file_extension.".".$combination.".dat"."\n";
	}
	unlink($fragmentfile.".".$file_extension.".".$combination.".dat");
    }
    foreach $suffix (@suffixes){
	$command = "cat ".$fragmentfile.".p".$suffix.".ccn.".$combination." >>".$fragmentfile.".".$file_extension.".".$combination.".dat";
	if($verbose_flag==1){
	    print $command."\n";
	}
	system($command);
    }
}


## create extension
$des_ext=".seqinfo.p";
for($y=1;$y<$k;$y++){
    $des_ext=$des_ext."p";
}
foreach $combination (@combinations){
    @projection_array=split(//,$combination);
    $x_label=$input_files[$projection_array[0]-1];
    $y_label=$input_files[$projection_array[1]-1];
    
    if($draft_flag==0){
    $argstring="perl home/Microsoftcoconut/finalscripts/comp_pairwise/create_gp_file_for_ps_kgenomes.pl ".$fragmentfile.".".$file_extension.".".$combination.".dat ".$k." ".$combination." ".$fragmentfile.".".$file_extension.".".$combination.".gp ".$forward_flag." ".$x_label." ".$y_label;    
}
    else{
	$argstring="home/Microsoftcoconut/finalscripts/comp_pairwise/create_gp_file_for_ps_k_draft_genomes.pl ".$fragmentfile.".".$file_extension.".".$combination.".dat ".$k." ".$combination." ".$fragmentfile.$des_ext." ".$fragmentfile.".".$file_extension.".".$combination.".gp ".$forward_flag." ".$x_label." ".$y_label; 
    }
    if($verbose_flag==1){
	print $argstring."\n";
    }
    system($argstring);
    $argstring="gnuplot ".$fragmentfile.".".$file_extension.".".$combination.".gp";
    if($verbose_flag==1){
	print $argstring."\n";

    }
    system($argstring); 





#open (GP, "|/usr/bin/gnuplot -persist") or die "no gnuplot";

print "\n \nHelooooooooooo H-Process fragments \n\n";

# force buffer to flush after each write
#use FileHandle;
#GP->autoflush(1);
#print GP "set term x11;plot 'testdata/chlamd/fragment.mm.pp.ccn.dat' \n";


#print GP "set term x11;plot '".$fragmentfile.".".$file_extension.".".$combination.".dat'\n";

print GP "set terminal postscript landscape color;plot 'home/Microsoftcoconut/testdata/chlamd/fragment.mm.pp.ccn.dat' index 0 title 'pm' w lp lt 2 lw 2 pt 1 ps .75, 'home/Microsoftcoconut/testdata/chlamd/fragment.mm.pm.ccn.dat' index 1 title 'pp' w lp lt 1 lw 2 pt 1 ps .75' "





















   
}

if($verbose_flag==1){
    print "Cleaning Up: Removing intermediate files\n";
}
#cleaning up:
foreach $combination (@combinations){
    foreach $suffix (@suffixes){
		
	if (-e $fragmentfile.".p".$suffix.".".$file_extension.".".$combination){
	    my $tmp_name=$fragmentfile.".p".$suffix.".".$file_extension.".".$combination;
	    if($verbose_flag==1){
		print  "Removing ".$tmp_name." \n";
	    }
	    unlink($fragmentfile.".p".$suffix.".".$file_extension.".".$combination);
	}
    }
}
}# end plotting check flag

} # end usechain_flag

exit;
######################################################
#CHAINING THE SECOND TIME

if($chaining_number!=2){
    #print "\n Done $chaining_number \n";
    exit;
}
if($verbose_flag==1){
    print "Chaining the second time\n";
}
#Note that $file_extension="ccn";


foreach $suffix (@suffixes){
#if (!(-e $fragmentfile.$suffix.".ccn.ccn")){
 $command = $chainer." ".$chaining_param[1]."  ".$fragmentfile.".p".$suffix.".".$file_extension;
 if($verbose_flag==1){
     print "$command\n";
 }
 system($command);
#}
}


$file_extension=$file_extension.".ccn";

if($plot_flag!=1){
    exit;
}

#cleaning files up:
$tmp_name=$fragmentfile.".".$file_extension.".dat";

if(-e $tmp_name){
    if($verbose_flag==1){
	print "Deleting".$tmp_name."\n";
    }
    unlink($tmp_name);    
}
foreach $combination (@combinations){
    if (-e $fragmentfile.".".$file_extension.".".$combination.".dat"){
	if($verbose_flag==1){
	    print "Removing: ".$fragmentfile.".".$file_extension.".".$combination.".dat"."\n";
	}
	unlink($fragmentfile.".".$file_extension.".".$combination.".dat");
    }
    foreach $suffix (@suffixes){
		
	if (-e $fragmentfile.".p".$suffix.".".$file_extension.".".$combination){
	    my $tmp_name=$fragmentfile.".p".$suffix.".".$file_extension.".".$combination;
	    if($verbose_flag==1){
		print  "Removing ".$tmp_name." \n";
	    }
	    unlink($fragmentfile.".p".$suffix.".".$file_extension.".".$combination);
	}
	if (-e ($fragmentfile.".p".$suffix.".".$file_extension.".dat")){	    
	    if($verbose_flag==1){
		print  "Removing ".$fragmentfile.".p".$suffix.".".$file_extension.".dat \n";
	    }
	    unlink($fragmentfile.".p".$suffix.".".$file_extension.".dat");
	}
    }
}
##########end cleaning up

#prepare plot (format transformations etc):
if($verbose_flag==1){
    print "Running ccn2gnuplot.pl\n";
}
if(-e $fragmentfile.".info"){
	$command = "perl home/Microsoftcoconut/finalscripts/ccn2gnuplot_multimat.pl ".$fragmentfile.".info"." ".$file_extension;
	if($verbose_flag==1){
	    print "$command\n";
	}
	system($command);
}
else {
	print "cant find infofile ".$fragmentfile.".info\n";
}




if(-e ($fragmentfile.".".$file_extension.".dat")){
	unlink($fragmentfile.".".$file_extension.".dat");
	print "Removing: ".$fragmentfile.".".$file_extension.".dat"."\n" ;
}

if (! (-e $fragmentfile.".".$file_extension.".dat")){
	foreach $suffix (@suffixes){
		$command = "cat ".$fragmentfile.".p".$suffix.".".$file_extension.".dat >>".$fragmentfile.".".$file_extension.".dat";
		if($verbose_flag==1){
		    print "$command\n";
		}
		system($command);
	#file is needed as input for chainer2permutation, which creates permutations from chainer output data
	}
}



foreach $combination (@combinations){
	if (-e $fragmentfile.".".$file_extension.".".$combination.".dat"){
		unlink($fragmentfile.".".$file_extension.".".$combination.".dat");
	}
	foreach $suffix (@suffixes){
	    $command = "cat ".$fragmentfile.".p".$suffix.".".$file_extension.".".$combination." >> ".$fragmentfile.".".$file_extension.".".$combination.".dat";
	    if($verbose_flag==1){
		print "$command\n";
	    }
	    system($command);
	}
}



## create extension
$des_ext="seqinfo";
for($y=0;$y<$k;$y++){
    $des_ext=$des_ext."p";
}

foreach $combination (@combinations){
    @projection_array=split(//,$combination);
    $x_label=$input_files[$projection_array[0]-1];
    $y_label=$input_files[$projection_array[1]-1];

    if(-e $fragmentfile.".".$file_extension.".".$combination.".gp"){
	unlink($fragmentfile.".".$file_extension.".".$combination.".gp");
    }
    if($draft_flag==0){
    $argstring="perl home/Microsoftcoconut/finalscripts/comp_pairwise/create_gp_file_for_ps_kgenomes.pl ".$fragmentfile.".".$file_extension.".".$combination.".dat ".$k." ".$combination." ".$fragmentfile.".".$file_extension.".".$combination.".gp"." ".$forward_flag." ".$x_label." ".$y_label;    
}
    else{
	$argstring="perl home/Microsoftcoconut/finalscripts/comp_pairwise/create_gp_file_for_ps_kgenomes.pl ".$fragmentfile.".".$file_extension.".".$combination.".dat ".$k." ".$combination." ".$fragmentfile.$des_ext." ".$fragmentfile.".".$file_extension.".".$combination.".gp"." ".$forward_flag." ".$x_label." ".$y_label;  
    }
    if($verbose_flag==1){
	print $argstring."\n";
    }
    system($argstring);
    $argstring="gnuplot ".$fragmentfile.".".$file_extension.".".$combination.".gp";

    if($verbose_flag==1){
	print $argstring."\n";
      
    }
    system($argstring);
    
}

if($verbose_flag==1){
    print "Cleaning up some intermediate files\n";
}
exit;
#cleaning up:
foreach $combination (@combination){
    foreach $suffix (@suffixes){
	if (-e $fragmentfile.".p".$suffix.".ccn.ccn.".$combination){
	    unlink($fragmentfile.".p".$suffix.".ccn.ccn.".$combination);
	}
    }
}

exit;


############### End program ###################################


sub get_number_of_contigs
{
my @arg   =@_;
my $in_fragment_file=$arg[0];
my $linecounter=0;

open(matchfile,"<".$in_fragment_file)             # filehandles don't use $name...
    or die "Oops! can't open input1 file ".$in_fragment_file."\n";

my @contig_num=();

#print $in_fragment_file."\n";

while ($infileline = <matchfile>)           # retrieve file, line by line
{
    # Kommentarzeile
    if($infileline=~/\#/){
	next;
    }
    else{
	$linecounter++;
    }
    
}
push(@contig_num, $linecounter); 
close(matchfile);
return @contig_num;
}

sub generate_info_file
{
#($input_files[0].".seqinfo",$input_files[1].".seqinfo")
my @arg   =@_;
my $in_file1=$arg[0];
my $in_file2=$arg[1];
my $basename=$arg[2];
my $linecounter=0;
my $outfile=$basename.".info";

open(seqinfofile1,"<".$in_file1)             # filehandles don't use $name...
    or die "Oops! can't open input2 file ".$in_file1."\n";



#print $in_fragment_file."\n";
my $genome1_size=0;
my $genome2_size=0;
while ($infileline = <seqinfofile1>)           # retrieve file, line by line
{
    # Kommentarzeile
    if($infileline=~/\#/){	
	next;
    }
    else{
	$infileline=~/\d+\s+(\d+)/;
	$length=$1;
	$genome1_size=$genome1_size+$length;
	$linecounter++;
    }
    
}
close(seqinfofile1);


open(seqinfofile2,"<".$in_file2)             # filehandles don't use $name...
    or die "Oops! can't open input3 file ".$in_file2."\n";

while ($infileline = <seqinfofile2>)           # retrieve file, line by line
{
    # Kommentarzeile
    if($infileline=~/\#/){
	next;
    }
    else{
	$infileline=~/\d+\s+(\d+)/;
	$length=$1;
	$genome2_size=$genome2_size+$length;
	$linecounter++;
    }
    
}


open(outfile,">".$outfile)             # filehandles don't use $name...
    or die "Oops! can't open input file ".$outfile."\n";
print outfile "2 $genome1_size $genome2_size\n";
print outfile $basename.".pp\n";
print outfile $basename.".pm\n";

close(outfile);


close(seqinfofile2);
}
