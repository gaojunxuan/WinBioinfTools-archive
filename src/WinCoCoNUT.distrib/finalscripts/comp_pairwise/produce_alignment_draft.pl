#!/usr/bin/perl
#

#$debug = 1;                  # if you're inclined to have debug levels...
#$| = 1;                      # flush output to stdout immediately

my $argc = @ARGV;

my $fragmentfile;
my $k=0; #number of involved sequences / genomes
my $verbose_flag=0;
my $plot_flag=0;
my $numofargs = scalar @ARGV;
my $draft_flag=0;
my $relative_flag=0;
my $align_param;    
my $align_param_flag=0;
my $usealign_flag=0;
my $forward_flag=0;

if ($argc < 4)
{
    print "Usage: \n"."produce_alignment_draft.pl align_prog clustalw_prog fragment_file_prefix verbose_flag parameter_file filter_ratio usealign_flag $forward_flag sequences\n";
    exit -1;
}

$align_prog=$ARGV[0];
$clustalw_prog=$ARGV[1];
$fragmentfile = $ARGV[2];
#$k = $ARGV[1];
$verbose_flag=$ARGV[3];
$filter_ration=$ARGV[5]; # filter ration determines if to plot or not
$usealign_flag=$ARGV[6]; # to skip re-aligning
$forward_flag==$ARGV[7];
if($filter_ration>0){
    $plot_flag=1;
}
else{
    $plot_flag=0;
}

my @queries = ();

my $sequence_string=" ";

for(my $i=8; $i<$numofargs; $i++)
{ 
    open(query_file,"<".$ARGV[$i])             # filehandles don't use $name...
	or die "Oops! can't open input file ".$ARGV[$i]."\n";    
    push(@queries,$ARGV[$i]);
    push(@input_files,$ARGV[$i]);
    $sequence_string=$sequence_string." ".$ARGV[$i];
    $k++;
}





# get directory of  align_program to set the invert_match_prog
my $alig_prog_dir= qx/dirname $align_prog/;
chop $alig_prog_dir;
my $inver_match_prog=$alig_prog_dir."/invert_matches";

#fragmentfile in directory und filename zerlegen
my $dirname = qx/dirname $fragmentfile/;
chop $dirname;
if($verbose_flag==1){
    print "produce_alignment_draft.pl: dirname is: $dirname \n";
}
my $basename = qx/basename $fragmentfile/;
chop $basename;
if($verbose_flag==1){
    print "produce_alignment_draft.pl: basename = $basename \n";
}


open(fragfile,"<".$fragmentfile)             # filehandles don't use $name...
	or die "Oops! can't open fragment file ".$fragmentfile."\n";    
close(fragfile);


# LOADING ALIGNING PARAMETERS IF A PARAMETER FILE IS GIVEN



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

    if($infileline=~/ALIGN=(.+)/){ # get path without endline
	$align_param=$1;
	$align_param_flag=1;

	last;
    }
    
}
close(param_file);

#print "HALLO Parameters of alignment: $align_param \n";
#exit;

if($align_param_flag==0){
    
    $align_param="";
    if($verbose_flag==1){
	print "Warning: Alignment parameters not given default will be taken \n";
    }
}


my $command;

my $argstring;


########## Detecting darft genomes##################################
@contig_num=get_number_of_contigs($fragmentfile);
    
if(($contig_num[0]>1) or ($contig_num[0]>1)){
    #modify the parameters
    $draft_flag=1;
    

    if($align_param=~/-draft/){
    }
    else{
	
	$align_param=" -draft ".$align_param;
	
    }
        
    #print "HALOO: ".$tmp_str."\n";    
	
    if($align_param=~/-relative/){
	$align_param=~s/-relative/ /;	
	$relative_flag=1;		
    }
    
}
else{
    
    
    if($align_param=~/-draft/){
	$align_param=~ s/-draft\s+/  /;	    
    }
    if($align_param=~/-relative/){	    
	$align_param=~s/-relative/ /;
	
    }
}

if($verbose_flag==1){
    print "Parameters of alignment: $align_param \n";
}

#exit;


 




#compute suffixes and combinations for k genomes:
$command = "generate_suffixes_k_genomes_multimat.pl ".$k ;
$result = qx($command);
#print $result;
@suffixes = split(/\s/, $result);

#$command = "generate_combinations_k_genomes.pl ".$k ;
#$result = qx($command);
$result = "1x2";
@combinations = split(/\s/,$result);


if($verbose_flag==1){
    print "Start alignment: calling alichainer\n";
}


my $file_extension="chn";
my $cntg_file_extension="chn.ctg";


if($usealign_flag==0){
# remove ordered files
foreach $suffix (@suffixes){
    if(-e $fragmentfile.".p$suffix.$file_extension.ordered"){	
	unlink($fragmentfile.".p$suffix.$file_extension.ordered");
	if($verbose_flag==1){
	    print "Removing: $fragmentfile.p$suffix.$file_extension.ordered  \n";
	}
	
    }
    if(-e $fragmentfile.".p$suffix.$cntg_file_extension.ordered"){	
	unlink($fragmentfile.".p$suffix.$cntg_file_extension.ordered");
	if($verbose_flag==1){
	    print "Removing: $fragmentfile.p$suffix.$cntg_file_extension.ordered  \n";
	}
	
    }
}

my $ready=0;

foreach $suffix (@suffixes){
    #print $suffix."\n";
    if (-e $fragmentfile.".p".$suffix.".".$file_extension){
	$command=$inver_match_prog."  ".$fragmentfile.".p".$suffix.".".$file_extension;
	if($verbose_flag==1){
	    print "Running: $command\n";
	}	
	system($command);
	if ($relative_flag){
	    $command=$inver_match_prog."  ".$fragmentfile.".p".$suffix.".".$cntg_file_extension;
	    if($verbose_flag==1){
		print "Running: $command\n";
	    }	
	    system($command);
	}
	if($relative_flag){
	    if($ready==0){
		$command = "$align_prog -clustalw $clustalw_prog ".$align_param."  -chain ".$fragmentfile.".p".$suffix.".".$file_extension.".ordered "." -relative ".$fragmentfile.".p".$suffix.".".$cntg_file_extension.".ordered  ". " -dirc p".$suffix."  -info ".$fragmentfile.".info". "  -f ".$sequence_string." > $fragmentfile.p$suffix.$file_extension.align";
		if($draft_flag==1){$ready=1;}
	    }
	    else{
		$command = "$align_prog -clustalw $clustalw_prog ".$align_param."  -ready -chain ".$fragmentfile.".p".$suffix.".".$file_extension.".ordered "." -relative ".$fragmentfile.".p".$suffix.".".$cntg_file_extension.".ordered  ". " -dirc p".$suffix."  -info ".$fragmentfile.".info". "  -f ".$sequence_string." > $fragmentfile.p$suffix.$file_extension.align";
	
	    }

	}
	else{
	    if($ready==0){
		$command = "$align_prog -clustalw $clustalw_prog ".$align_param."  -chain ".$fragmentfile.".p".$suffix.".".$file_extension.".ordered "." -dirc p".$suffix."  -info ".$fragmentfile.".info". "  -f ".$sequence_string." > $fragmentfile.p$suffix.$file_extension.align";
		if($draft_flag==1){$ready=1;}
	    }
	    else{
		$command = "$align_prog -clustalw $clustalw_prog ".$align_param."  -ready  -chain ".$fragmentfile.".p".$suffix.".".$file_extension.".ordered "." -dirc p".$suffix."  -info ".$fragmentfile.".info". "  -f ".$sequence_string." > $fragmentfile.p$suffix.$file_extension.align";
		
	    }

	}
	if($verbose_flag==1){
	    print "Running $command\n";
	}
        system($command);
    }
}

if($verbose_flag==1){
    print "Removing temporary files: tmpffXg, tmpffXg.gde, tmpffXg.dnd \n";
    if(-e "tmpffXg"){
	$command = unlink("tmpffXg"); 
    }
    if(-e "tmpffXg.gde"){
	$command = unlink("tmpffXg.gde"); 
    }
    if(-e "tmpffXg.dnd"){
	$command = unlink("tmpffXg.dnd"); 
    }

}

# remove ordered files
foreach $suffix (@suffixes){
    if(-e $fragmentfile.".p$suffix.$file_extension.ordered"){	
	unlink($fragmentfile.".p$suffix.$file_extension.ordered");
	if($verbose_flag==1){
	    print "Removing: $fragmentfile.p$suffix.$file_extension.ordered  \n";
	}
	
    }
    if(-e $fragmentfile.".p$suffix.$cntg_file_extension.ordered"){	
	unlink($fragmentfile.".p$suffix.$cntg_file_extension.ordered");
	if($verbose_flag==1){
	    print "Removing: $fragmentfile.p$suffix.$cntg_file_extension.ordered  \n";
	}
	
    }
}

}
#######################################################################################
#### Plotting the aligned sequences
my $file_extension="ccn.filtered";
$filter_prog=$alig_prog_dir."/filteralign.x";

#$plot_flag=1;
if($plot_flag==1){


#$filter_ration=0.5;

# remove filtered files
foreach $suffix (@suffixes){
    if(-e $fragmentfile.".p$suffix.$file_extension"){	
	unlink($fragmentfile.".p$suffix.$file_extension");
	if($verbose_flag==1){
	    print "Removing: $fragmentfile.p$suffix.$file_extension \n";
	}
	
    }
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


foreach $suffix (@suffixes){
    
    if (-e $fragmentfile.".p".$suffix){
	
	$command = "$filter_prog  ".$fragmentfile.".p".$suffix."  ".$filter_ration;
	
	if($verbose_flag==1){
	    print "Running $command\n";
	}
        system($command);
    }
    
}

#prepare plot (format transformations etc):
    if($verbose_flag==1){
	print "Running ccn2gnuplot_multimat.pl\n";
    }


    if(-e $fragmentfile.".info"){
	$command = "ccn2gnuplot_multimat.pl ".$fragmentfile.".info ".$file_extension;
	if($verbose_flag==1){
	    print "command = $command\n";
	}
	system($command);
    }
    else {
	print "cant find infofile ".$fragmentfile.".info\n";
	exit;
    }

#debug
#foreach $suffix (@suffixes){
#    print "HIIIIIIIIIIII $suffix"."\n";
#}


$tmp_name=$fragmentfile.".".$file_extension.".dat";
if(-e $tmp_name){
    if($verbose_flag==1){
	print "Deleting".$tmp_name."\n";
    }
    unlink($tmp_name);
}

if (! (-e $tmp_name)){
	foreach $suffix (@suffixes){	    
	    $command = "cat ".$fragmentfile.".p".$suffix.".".$file_extension.".dat >>".$fragmentfile.".".$file_extension.".dat";
	    if($verbose_flag==1){
		print $command."\n";
	    }
            system($command);
	#file is needed as input for chainer2permutation, which creates permutations from chainer output data
	}
}



foreach $combination (@combinations){
	if (-e $fragmentfile.".".$file_extension.".".$combination.".dat"){
	    if($verbose_flag==1){
		print "Unlink".$fragmentfile.".".$file_extension.".".$combination.".dat"."\n";
	    }
	    unlink($fragmentfile.".".$file_extension.".".$combination.".dat");
	}
	foreach $suffix (@suffixes){
	    $command = "cat ".$fragmentfile.".p".$suffix.".".$file_extension.".".$combination." >>".$fragmentfile.".".$file_extension.".".$combination.".dat";
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
    system($argstring);
    $argstring="gnuplot ".$fragmentfile.".".$file_extension.".".$combination.".gp";
    if($verbose_flag==1){
	print $argstring."\n";
    }
    system($argstring);
    
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





exit; ########### End program


#################### Extra Code for handling Draft ####################

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
	    #print $infileline;
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

