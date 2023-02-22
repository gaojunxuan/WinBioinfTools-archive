#!/usr/bin/perl
#

my $argc = @ARGV;

my $fragmentfile;
my $k; #number of involved sequences / genomes
my $verbose_flag=0;
my $plot_flag=0;

if ($argc < 4)
{
    print "Usage: \n"."second_chaining.pl fragment_file number_of_genomes chaining_prog verbose_flag parameter_file plot_flag forward_flag usematch_flag usechain_flag file_extension in_seq_string\n";
    exit -1;
}


$fragmentfile = $ARGV[0];
$k = $ARGV[1];
$chainer = $ARGV[2];
$verbose_flag=$ARGV[3];
$plot_flag=$ARGV[5];
$forward_flag=$ARGV[6];
$usematch_flag=$ARGV[7];
$usechain_flag=$ARGV[8];
$in_file_extension=$ARGV[9];

my @input_files;

if($plot_flag==1){
    for(my $i=10; $i<$argc; $i++){ 	    		
	push(@input_files,$ARGV[$i]);	
    }
     
}

#fragmentfile in directory und filename zerlegen
my $dirname = qx/dirname $fragmentfile/;
chop $dirname;
if($verbose_flag==1){
    print "iterate_chaining.pl: dirname is: $dirname \n";
}
my $basename = qx/basename $fragmentfile/;
chop $basename;
if($verbose_flag==1){
    print "iterate_chaining.pl: basename = $basename \n";
}


open(fragfile,"<".$fragmentfile)             # filehandles don't use $name...
	or die "Oops! can't open fragment file ".$fragmentfile."\n";    
close(fragfile);


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
	    if($verbose_flag==1){
		print "Parameters of chaining call no. ".$chaining_number." :".$chaining_param[($chaining_number-1)]."\n";
	    }
	}

    }
    close(param_file);
}


@contig_num=get_number_of_contigs($fragmentfile);
    
if(($contig_num[0]>1) or ($contig_num[0]>1)){
    #modify the parameters
    $draft_flag=1;
    for(my $x=0;$x<$chaining_number;$x++){
	$tmp_str=$chaining_param[$x];
	
	if($tmp_str=~/-l\s+/){	    
	    $tmp_str=~ s/-l\s+/ -draft /;	    
	    $chaining_param[$x]=$tmp_str;
	}
    }	   
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

if(($forward_flag==0) or ($argc < 7) ){
#compute suffixes and combinations for k genomes:
    $command = "generate_suffixes_k_genomes_multimat.pl ".$k ;
    $result = qx($command);
    @suffixes = split(/\s/, $result);
     
}
else{
    $result="p";
    for($x=2;$x<$k;$x++){
	$result=$result."p";
    }
    @suffixes = split(/\s/, $result);
}

$command = "generate_combinations_k_genomes.pl ".$k ;
$result = qx($command);
@combinations = split(/\s/,$result);



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

$file_extension=$in_file_extension;

foreach $suffix (@suffixes){
#if (!(-e $fragmentfile.$suffix.".ccn.ccn")){
    if($draft_flag){
	$command = $chainer." ".$chaining_param[1]."  -des ".$fragmentfile.".seqinfo.p$suffix  ".$fragmentfile.".p".$suffix.".".$file_extension;
    }
    else{
	$command = $chainer." ".$chaining_param[1]."  ".$fragmentfile.".p".$suffix.".".$file_extension;
    }
    if($verbose_flag==1){
	print "$command\n";
    }
    $retcode=system($command);
    if($retcode ne 0){
	print STDERR "\nFailure: $command\n";
	exit -1;
    }
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

if(-e $fragmentfile.".info"){
	$command = "ccn2gnuplot_multimat.pl ".$fragmentfile.".info"." ".$file_extension;
	if($verbose_flag==1){
	    print "$command\n";
	}
	$retcode=system($command);
	if($retcode ne 0){
	    print STDERR "\nFailure: $command\n";
	    exit -1;
	}
}
else {
	print "Error: can't find infofile ".$fragmentfile.".info\n";
}



if(-e $fragmentfile.".".$file_extension.".dat"){
	unlink($fragmentfile.".".$file_extension.".dat");
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
$des_ext=".seqinfo.p";
for($y=1;$y<$k;$y++){
    $des_ext=$des_ext."p";
}


foreach $combination (@combinations){
    @projection_array=split(/x/,$combination);
    $x_label=$input_files[$projection_array[0]-1];
    $y_label=$input_files[$projection_array[1]-1]; 

    if($draft_flag==0){	
	$argstring="create_gp_file_for_ps_kgenomes.pl ".$fragmentfile.".".$file_extension.".".$combination.".dat ".$k." ".$combination." ".$fragmentfile.".".$file_extension.".".$combination.".gp ".$forward_flag." ".$x_label." ".$y_label;    
    }
    else{
	$argstring="create_gp_file_for_ps_k_draft_genomes.pl ".$fragmentfile.".".$file_extension.".".$combination.".dat ".$k." ".$combination." ".$fragmentfile.$des_ext." ".$fragmentfile.".".$file_extension.".".$combination.".gp ".$forward_flag." ".$x_label." ".$y_label; 
    }

    if($verbose_flag==1){
	print $argstring."\n";
    }
    $retcode=system($argstring);
    if($retcode ne 0){
	print STDERR "\nFailure: $argstring\n";
	exit -1;
    }

    $argstring="gnuplot ".$fragmentfile.".".$file_extension.".".$combination.".gp";

    if($verbose_flag==1){


	print $argstring."\n";
    }
    $retcode=system($argstring);
    if($retcode ne 0){
	print STDERR "\nFailure: $argstring\n";
	exit -1;
    }
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

################## End program


sub get_number_of_contigs
{
my @arg   =@_;
my $in_fragment_file=$arg[0];
my $linecounter=0;

open(matchfile,"<".$in_fragment_file)             # filehandles don't use $name...
    or die "Oops! can't open input file ".$infilename."\n";

my @contig_num=();

#print $in_fragment_file."\n";

while ($infileline = <matchfile>)           # retrieve file, line by line
{
    $linecounter++;
    # Kommentarzeile
    if($infileline=~/\#/){	
	if ($infileline=~/databaselength=\s*(\d+)/){
	    print $infileline;
	    my $length = $1;
	    #print "length = $length\n";
	    
	    if($infileline=~/including\s+(\d+)/){
		push(@contig_num, $1+1); 
	    }
	    else{
		push(@contig_num, 1); 
		#print "contig num: ".$contig_num[0]."\n";
	    }
	}
	if ($infileline=~/querylength=\s*(\d+)/){
	    my $length = $1;
	    #print $infileline;
	    #print "length = $length\n";
	    
	    if($infileline=~/including\s+(\d+)/){
		push(@contig_num, $1+1); 
	    }
	    else{
		push(@contig_num, 1); 
		#print "contig num: ".$contig_num[0]."\n";
	    }
	}
    }
    else{	
	last;
    }
    
}
#print "contig num: ".$contig_num[0]."\n";
close(matchfile);
return @contig_num;
}
