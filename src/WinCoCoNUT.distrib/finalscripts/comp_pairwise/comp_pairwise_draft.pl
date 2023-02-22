#!/usr/bin/perl
#use strict;
# use warnings;

my $numofargs = scalar @ARGV;

my $syntenic_flag=0;
my $synteny_param;

if($numofargs < 4)
{
  print STDERR "Usage: perl coconut.pl   <Options> seq_1 seq_2 \n";

#  print STDERR "Arguments:\n";
#  print STDERR "   -cf        : confiuration file, containing paths for required modules\n";
  print STDERR "\nOptions:\n";
  print STDERR "  -pr        : parameter file (optional), if not given defaults are computed\n";
  print STDERR "  -v         : verbose mode, i.e., deisplay of the program steps \n";
  print STDERR "  -forward   : run the comparison for forward strands only\n";
 # print STDERR "  -align     : compute alignment using clustalw\n";
  print STDERR "  -plot      : produce Postscript 2D plots of the chains\n";
 # print STDERR "  -plotali X : filter out alignments with idenetity < X (0< X <= 1) and produce 2D plots \n";
  print STDERR "  -indexname : specify the index, if constructed\n";
  print STDERR "  -useindex  : do not construct index again\n";
  print STDERR "  -usematch  : do not compute matches again, this construct no index\n";
 # print STDERR "  -usechain  : use the computed chains and complete processing\n";
 # print STDERR "  -usealign  : use the computed alignments and complete processing\n";
  print STDERR "  -prefix    : specify a prefix name for the output files\n";
 # print STDERR "  -syntenic  : computes sysntenic regions by removing repeats and applying\n";
  #print STDERR "	       1D chaining over all dimensions.\n";
  print STDERR "\nExample 1: \n";
  print STDERR "  >coconut.pl -pairwise -v seq1.fna seq2.fna -plot\n";
#  print STDERR "\nExample 2 (without coconut interface): \n";
#  print STDERR "   FPATH/comp_pairewise_draft.pl -cf Home/MicroSoftCoCoNUT/config -v seq1.fna seq2.fna -plot -align -syntenic\n";
#  print STDERR "   FPATH: is the path to the script finalscripts/comp_pairwise\n";
  exit -1;
}

my $fragment_prog;
my $fragment_prog_flag=0;
my $config_file_flag=0;
my $config_file;
my $param_file;
my $param_file_flag=0;
my $verbose_flag=0;
my $plot_flag=0;
my $clustalw_flag=0;
my $forward_flag=0;
my $draft_flag=0;
my $useindex_flag=0;
my $usematch_flag=0;
my $usechain_flag=0;
my $usealign_flag=0;

my $plotali_flag=0;

my $prefix_flag=0;
my $prefix="./fragment.mm";

my $in_seq_string=" ";

my @queries = ();

my $filter_ratio=0;

my $index_name_flag=0;
my $index_name="";
#print $numofargs."\n";


   
$fragment_prog="gt";
$fragment_prog_flag=1;
	

for(my $i=0; $i<$numofargs; $i++)
{    
  
 
    if($ARGV[$i]  eq "-cf"){
	$i++;
	$config_file=$ARGV[$i];
	$config_file_flag=1;
	next;
    }
    if($ARGV[$i] eq "-pr"){
	$i++;
	$param_file=$ARGV[$i];
	$param_file_flag=1;
	next;
    }
    if($ARGV[$i] eq "-indexname"){
	$i++;
	$index_name=$ARGV[$i];
	$index_name_flag=1;
	$useindex_flag=1;
	next;
    }
    if($ARGV[$i] eq "-plotali"){
	$i++;
	$filter_ratio=$ARGV[$i];
	
	
	$x=$filter_ratio;
	$numres=is_number($x);
	
	if(($numres !=1)||($x>1)||($x<=0)){
	    print STDERR "The option -plotali requires a float number > 0 and <= 1\n";
	    exit -1;
	}	
	$plotali_flag=1;
	next;
    }
    
    if($ARGV[$i] eq "-v"){
	$verbose_flag=1;
	next;
    }
    if($ARGV[$i] eq "-plot"){
	$plot_flag=1;
	next;
    }
    if($ARGV[$i] eq "-align"){
	$clustalw_flag=1;
	next;
    }
    if($ARGV[$i] eq "-useindex"){
	$useindex_flag=1;
	next;
    }
    if($ARGV[$i] eq "-usematch"){
	$usematch_flag=1;
	$useindex_flag=1;
	next;
    }
    if($ARGV[$i] eq "-usechain"){
	$usematch_flag=1;
	$useindex_flag=1;
	$usechain_flag=1;
	next;
    }
    if($ARGV[$i] eq "-usealign"){
	$usematch_flag=1;
	$useindex_flag=1;
	$usechain_flag=1;
	$usealign_flag=1;
	next;
    }
    if($ARGV[$i] eq "-syntenic"){
	$syntenic_flag=1;
	next;
    }
    if($ARGV[$i] eq "-prefix"){
	$i++;
	$prefix=$ARGV[$i];
	$prefix_flag=1;	
	next;
    }
    if($ARGV[$i] eq "-forward"){
	$forward_flag=1;
	next;
    }
    if((scalar @queries)<2){
	# check if the query file is valid
	open(query_file,"<".$ARGV[$i]) 
	    or die "Oops! can't open input file ".$ARGV[$i]."\n";    
	push(@queries,$ARGV[$i]);
	$in_seq_string= $in_seq_string." ".$ARGV[$i];
    }
}


my $dir_name = qx/dirname $queries[0]/;
chop $dir_name;
#print "dirname is: $dir_name \n";

# CHECK ARGUMENTS

if($fragment_prog_flag==0){
    print STDERR "Error: Fragment generation progam should be specified\n";
    exit -1;
}
if($config_file_flag==0){
    print STDERR "Error: config file should be specified\n";
    exit -1;
}
if($param_file_flag==0){
    print "Warning: Parameter file not specified. Default will be taken \n";
    $param_file=generate_parameter_file(@queries);    
    #exit -1;
}
else{
    $clustalw_flag=0;
    $syntenic_flag=0;  
    $forward_flag=0;
}




# LOAD CONFIG
# get paths for fragment generation and chaining
#my$config_file="Home/MicroSoftCoCoNUT/finalscripts/config";

open(config_file,"<".$config_file)             # filehandles don't use $name...
    or die "Error: can't open config file \n";    

my $infileline="";
my $fragment_generation_path;
my $chaining_path;
my $align_path;
while ($infileline = <config_file>)           # retrieve file, line by line
{    	#chop $infileline;
	#print $infileline;

    	$tmp_infileline=$infileline;
	@tmp_tokens=split(" ",$tmp_infileline);
	@tmp_tokens=split('',$tmp_tokens[0]);
	if($tmp_tokens[0] eq "#"){next;}

	if($infileline=~/FRAGMENT=(.+)/){ # get path withoit endline
	    #print $infileline;
	    #$fragment_generation_path=$1;
           $fragment_generation_path="home/Microsoftcoconut/bin/genometools-current/bin";

	    #print $fragment_generation_path."HHHH\n";
	}
	
	elsif($infileline=~/CHAINING=(.+)/){ # get path without endline
	    #print $infileline;
	    #$chaining_path=$1;
	    $chaining_path="home/Microsoftcoconut/bin/chainer";		
	    #print $1."HHHHH\n";
	}
	elsif($infileline=~/ALIGN=(.+)/){ # get path without endline
	    $align_path=$1;	    	    
	}
	
}
close(config_file);


$prog_frag_gen=$fragment_generation_path."/".$fragment_prog;
if($verbose_flag==1){
    print "Fragment generation program: ".$prog_frag_gen."\n";
}
$prog_chainer=$chaining_path."/chainer";
if($verbose_flag==1){
    print "Chaining program: ".$prog_chainer."\n";
}


my $prog_synteny=$chaining_path."/chainer2permutation.x";

# LOAD PARAMETERS
#my $param_file="/home/mibrahim/Myprojects/CoCoNUT/finalscripts/parameters";

open(param_file,"<".$param_file)             # filehandles don't use $name...
    or die "Error: can't open parameter file \n";    

my $fragment_param;
my @chaining_param;
my $chaining_number=0;

while ($infileline = <param_file>)           # retrieve file, line by line
{    	#chop $infileline;
	#print $infileline;
    	# check for commented lines starting with #

	$tmp_infileline=$infileline;
	@tmp_tokens=split(" ",$tmp_infileline);
	@tmp_tokens=split('',$tmp_tokens[0]);
	if($tmp_tokens[0] eq "#"){next;}

	if($infileline=~/FRAGMENT=(.+)/){ # get path withoit endline
	    $fragment_param=$1;

	}	
	elsif($infileline=~/CHAINING=(.+)/){ # get path without endline
	    $chaining_param[$chaining_number]=$1;
	    $chaining_number++;
	}
	elsif($infileline=~/ALIGN=(.+)/){ # get path without endline	    
	    $clustalw_flag=1;
	}
	elsif($infileline=~/SYNTENY=(.+)/){ # get path without endline
	    $syntenic_flag=1;	
    	    $synteny_param=$1;
	}

}
close(param_file);



if($clustalw_flag==1){
    $align_prog=$align_path."/alichainer";
    $clustalw_prog=$align_path."/clustalw";
    if(!(-e $align_prog)){
	print STDERR "Error: can't find $align_prog \n";
	exit -1;
    }
    if(!(-e $clustalw_prog)){
	print STDERR "Error: can't find $clustalw_prog \n";
	exit -1;
    }
    if($verbose_flag==1){
	print "Aigninging program: $align_prog \n";
	print "Clustalw program: $clustalw_prog \n";
    }
}
else{
    $plotali_flag=0;
}

if($verbose_flag==1){
    print "Align program: ".$align_prog."\n";
}


my $numofqueries = 2;

if($verbose_flag==1){
    for(my $i=0; $i<$numofqueries; $i++){
	print "Query file no. ".$i.": ".$queries[$i]."\n";
    }

    print "Fragment generation parameters: ".$fragment_param."\n";
    print "Number of recursive chaining calls: ".$chaining_number."\n";

    for(my $i=0; $i<$chaining_number; $i++){
	print "Parameters of chaining call no. ".($i+1) .": ".$chaining_param[$i]."\n";
    }
}



#*********************************************************###
#          The Fragment Generation Phase                  ###
#*********************************************************###

my $argstring;
#my $index_name=$dir_name."/Index";

if($fragment_param=~/indexname(\w)/){
    $index_name=$1;
    $index_name_flag=1;
}
else{
    if($index_name_flag==0){
	$index_name=$dir_name."/Index";
	if($prefix_flag==1){
	    $index_name=$prefix.".index";
	}    
    } 
}





my $fragment_file;
if($prefix_flag==0){
    $fragment_file= $dir_name."/fragment.mm";
}
else{
   $fragment_file=$prefix;
}


if($fragment_prog eq "gt"){
    $prog_indexing=$fragment_generation_path."/gt suffixerator";
    my $query_str="";
 #   for(my $i=0; $i<$numofqueries; $i++){
#	$query_str=$query_str." ".$queries[$i];
 #   }
    
    
    if($useindex_flag==0){
	#$argstring=$chaining_path."/concatmultifasta $queries[0] > $queries[0].seqinfo";	
       $argstring="Home/MicroSoftCoCoNUT/bin/chainer/concatmultifasta $queries[0] > $queries[0].seqinfo";

	if($verbose_flag==1){
	    print "Running 1 ".$argstring."\n";
	}
	my $retcode = system($argstring);
	if($retcode ne 0){
	    print STDERR "\nFailure: $argstring\n";
	    exit -1;
	}
	if($index_name_flag==0){	    
	    $argstring=$prog_indexing." -pl -suf -bwt -lcp -tis -bck  -dna  -v -indexname ".$index_name." -db ".$queries[0].".ready";
	}
	else{	
	    $argstring=$prog_indexing." -pl -suf -bwt -lcp -tis -bck  -dna  -v -indexname ".$index_name." -db ".$queries[0].".ready";
	}
	if($verbose_flag==1){
	    print "Running ".$argstring."\n";
	}
	my $retcode = system($argstring);
	if($retcode ne 0){
	    print STDERR "\nFailure: $argstring\n";
	    exit -1;
	}
    }
    
    if($usematch_flag==0){

	$argstring=$chaining_path."/concatmultifasta $queries[1] > $queries[1].seqinfo";
	if($verbose_flag==1){
	    print "Running ".$argstring."\n";
	}
	my $retcode = system($argstring);
	if($retcode ne 0){
	    print STDERR "\nFailure: $argstring\n";
	    exit -1;
	}
	# forward fragments
	$argstring=$prog_frag_gen." dev maxpairs ".$fragment_param." -q $queries[1] -ii ".$index_name."  > ".$fragment_file."+";
	if($verbose_flag==1){
	    print "Running ".$argstring."\n";
	}
	my $retcode = system($argstring);
	if($retcode ne 0){
	    print STDERR "Failure: $argstring\n";
	    exit -1;
	}
	# reverse complement fragments
	$argstring=$chaining_path."/reversecomp $queries[1].ready > $queries[1].ready.rc";
	if($verbose_flag==1){
	    print "Running ".$argstring."\n";
	}
	my $retcode = system($argstring);
	if($retcode ne 0){
	    print STDERR "\nFailure: $argstring\n";
	    exit -1;
	}
	$argstring=$prog_frag_gen." dev maxpairs ".$fragment_param." -q $queries[1].ready.rc -ii ".$index_name."  > ".$fragment_file."-";
	if($verbose_flag==1){
	    print "Running ".$argstring."\n";
	}
	my $retcode = system($argstring);
	if($retcode ne 0){
	    print STDERR "Failure: $argstring\n";
	    exit -1;
	}
	

    }
    
    if($verbose_flag==1){
	print "Fragments are stored in:  ".$fragment_file."\n";
	print "Done :-) \n";
	print "CHAINING PHASE \n";
    }
    
    # for draft: check for multiple contigs in the header of the fragment file
    
    $contig_num1=get_number_of_contigs($queries[0].".seqinfo");
    push(@contig_num,$contig_num1);
    $contig_num=get_number_of_contigs($queries[1].".seqinfo");
    push(@contig_num,$contig_num1);
    print "Contig numbers: ".$contig_num[0]." ".$contig_num[1]."\n";


    if(($contig_num[0]>1) or ($contig_num[0]>1)){
	generate_seqinfo_files(@queries);
	 $draft_flag=1;
    }
    


    # CHAINING PHASE and Plotting PHASE
    
    $argstring="perl Home/MicroSoftCoCoNUT/finalscripts/comp_pairwise/process_fragments_2genomes_vmatch.pl  ".$fragment_file."  ".$numofqueries."  ".$prog_chainer. "  ".$verbose_flag."  ".$param_file. "  ".$plot_flag." ".$forward_flag." ".$usematch_flag." ".$usechain_flag." "." ".$in_seq_string;
    
    
    if($verbose_flag==1){
	print "Running ".$argstring."\n";
    }
    
    my $retcode = system($argstring);
    if($retcode ne 0){
	print STDERR "Failure: $argstring\n";
	exit -1;
    }
    
    if($verbose_flag==1){
	print "Chains are stored in:  ".$dir_name."\n";
	#print "Done :-) \n";
    }


    # ALIGNMENT PHASE
    cleanup();
    
    if($clustalw_flag==1){
	$argstring="produce_alignment_draft.pl $align_prog $clustalw_prog $fragment_file $verbose_flag $param_file  $filter_ratio $usealign_flag $forward_flag $in_seq_string";
	
	if($verbose_flag==1){
	    print "Running ".$argstring."\n";
	}    
	my $retcode = system($argstring);
	if($retcode ne 0){
	    print STDERR "Failure: $argstring\n";
	    exit -1;
	}
    }
    
    ########### SECOND CHAINING

    if($chaining_number==2){
	if($plotali_flag==0){
	    $argstring="iterate_chaining.pl  ".$fragment_file."  ".$numofqueries."  ".$prog_chainer. "  ".$verbose_flag."  ".$param_file. "  ".$plot_flag." ".$forward_flag." ".$usematch_flag." ".$usechain_flag." "." "." ccn ".$in_seq_string;
    }
	else{
	    $argstring="iterate_chaining.pl  ".$fragment_file."  ".$numofqueries."  ".$prog_chainer. "  ".$verbose_flag."  ".$param_file. "  ".$plotali_flag." ".$forward_flag." ".$usematch_flag." ".$usechain_flag." "." "." ccn.filtered ".$in_seq_string;
	    
	}
	if($verbose_flag==1){
	    print "Running ".$argstring."\n";
	}
	#exit;
	my $retcode = system($argstring);
	if($retcode ne 0){
	    print STDERR "failure: $argstring\n";
	    exit 1;
	}
    
	if($verbose_flag==1){
	    print "Recursive Chains are stored in:  ".$dir_name."\n";
	    print "Done :-) \n";
	}
    
    }
    ##########end second chaining
    



}

####################################################################
##                          Handling Synteny                      ##
####################################################################


#print "Handling synteny: $draft_flag, $syntenic_flag  $numofqueries $synteny_param \n";  

if(($draft_flag==1)&&($syntenic_flag==1)){
    $synteny_param=$synteny_param." -draft ";
    for($i=0;$i<$numofqueries; $i++){
	$synteny_param=$synteny_param." ".$queries[$i].".seqinfo ";
    }
        #print "Warning: -syntenic option is not performed because the genomes are draft\n";    
}

#if(($draft_flag==0)&&($syntenic_flag==1)){
if(($syntenic_flag==1)){
    if(($plot_flag!=1) && ($plotali_flag==0)){
	my $file_extension="ccn";

	$tmp_name=$fragment_file.".*".$file_extension.".dat";			
	$command="rm  -f $tmp_name"; 
	if($verbose_flag==1){
	    print $command."\n";
	}
	system($command);	

	
	if(-e $fragment_file.".info"){
	    $command = "ccn2gnuplot_multimat.pl ".$fragment_file.".info ".$file_extension;
	    if($verbose_flag==1){
		print "command = $command\n";
	    }
	    $retcode=system($command);
	    if($retcode ne 0){
		print STDERR "\nFailure: $command\n";
		exit -1;
	    }
	}
	else {
	    print STDERR "Error: can't find infofile ".$fragment_file.".info\n";
	    exit -1;
	}
	$command = "generate_suffixes_k_genomes_multimat.pl ".$numofqueries;
	$result = qx($command);
	@suffixes = split(/\s/, $result);
	foreach $suffix (@suffixes){	    
	    $command = "cat ".$fragment_file.".p".$suffix.".".$file_extension.".dat >>".$fragment_file.".".$file_extension.".dat";
	    if($verbose_flag==1){
		print $command."\n";
	    }
	    $retcode=system($command);
	    if($retcode ne 0){
		print STDERR "\nFailure: $command\n";
		exit -1;
	    }
	    #file is needed as input for chainer2permutation, which creates permutations from chainer output data
	}
	# preparing for case of second chaining phase
	if($chaining_number>1){	    
	    if($plotali_flag==0){
		$file_extension2="ccn.ccn";
	    }
	    else{
		$file_extension2="ccn.filtered.ccn";
	    }
	    
	    $tmp_name=$fragment_file.".*".$file_extension2.".dat";			
	    $command="rm  -f $tmp_name"; 
	    if($verbose_flag==1){
		print $command."\n";
	    }
	    system($command);
	    
	    if(-e $fragment_file.".info"){
		$command = "ccn2gnuplot_multimat.pl ".$fragment_file.".info ".$file_extension2;
		if($verbose_flag==1){
		    print "command = $command\n";
		}
		$retcode=system($command);
		if($retcode ne 0){
		    print STDERR "\nFailure: $command\n";
		    exit -1;
		}
	    }
	    else {
		print "Error: can't find infofile ".$fragment_file.".info\n";
		exit;
	    }
	    foreach $suffix (@suffixes){	    
		$command = "cat ".$fragment_file.".p".$suffix.".".$file_extension2.".dat >>".$fragment_file.".".$file_extension2.".dat";
		if($verbose_flag==1){
		    print $command."\n";
		}
		system($command);
		#file is needed as input for chainer2permutation, which creates permutations from chainer output data
	    }
	}
	# preparing for the case of plot_ali	
	if($plotali_flag==1){	    
	    $file_extension2="ccn.filtered";
	    
	    $tmp_name=$fragment_file.".*".$file_extension2.".dat";
	    $command="rm  -f $tmp_name"; 
	    if($verbose_flag==1){
		print $command."\n";
	    }
	    system($command);
	    
	    if(-e $fragment_file.".info"){
		$command = "ccn2gnuplot_multimat.pl ".$fragment_file.".info ".$file_extension2;
		if($verbose_flag==1){
		    print "command = $command\n";
		}
		$retcode=system($command);
		if($retcode ne 0){
		    print STDERR "\nFailure: $command\n";
		    exit -1;
		}
	    }
	    else {
		print "Error: can't find infofile ".$fragment_file.".info\n";
		exit;
	    }
	    foreach $suffix (@suffixes){	    
		$command = "cat ".$fragment_file.".p".$suffix.".".$file_extension2.".dat >>".$fragment_file.".".$file_extension2.".dat";
		if($verbose_flag==1){
		    print $command."\n";
		}
		system($command);
		#file is needed as input for chainer2permutation, which creates permutations from chainer output data
	    }
	}

    }

   
    if($plotali_flag==0){
	$argstring="$prog_synteny $fragment_file.ccn.dat $numofqueries $synteny_param > $fragment_file.ccn.syn";
	if($verbose_flag==1){
	    print "Running ".$argstring."\n";
	} 
	my $retcode = system($argstring);
	if($retcode ne 0){
	    print STDERR "Failure: $argstring\n";
	    exit -1;
	}
    }
    if($plotali_flag==1){
	$argstring="$prog_synteny $fragment_file.ccn.filtered.dat $numofqueries $synteny_param > $fragment_file.ccn.filtered.syn";
	if($verbose_flag==1){
	    print "Running ".$argstring."\n";
	} 
	my $retcode = system($argstring);
	if($retcode ne 0){
	    print STDERR "Failure: $argstring\n";
	    exit -1;
	}
    }
    if($chaining_number>1){
	if($plotali_flag==0){
	    $file_extension2="ccn.ccn";
	}
	else{
	    $file_extension2="ccn.filtered.ccn";
	}
	$argstring="$prog_synteny $fragment_file.$file_extension2.dat $numofqueries $synteny_param > $fragment_file.$file_extension2.syn";
	if($verbose_flag==1){
	    print "Running ".$argstring."\n";
	} 
	my $retcode = system($argstring);
	if($retcode ne 0){
	    print STDERR "Failure: $argstring\n";
	    exit -1;
	}
    }
 

    #$command = "generate_combinations_k_genomes.pl ".$numofqueries ;
    $result ="1x2";
    @combinations = split(/\s/,$result);
    if($plotali_flag==0){
	$synteny_prefix=$fragment_file.".ccn.dat.syn.";
	foreach $combination (@combinations){
	    @projection_array=split(//,$combination);
	    $tmp_f=$synteny_prefix.$combination.".gp";
	    open (SynFile, "<$tmp_f") or die "Can't open $tmp_f \n";
	    my @fsynlines = <SynFile>;
	    close SynFile;	
	    open (SynFile, ">$tmp_f") or die "Can't open $tmp_f \n";
	    foreach $fline (@fsynlines){
		if($fline=~/xlabel/){
		    $fline="set xlabel \"$queries[$projection_array[0]-1]\"\n";
		}
		if($fline=~/ylabel/){
		    $fline="set ylabel \"$queries[$projection_array[1]-1]\"\n";
		}
		print SynFile $fline;
	    }
	    close SynFile;
	    $argstring="gnuplot ".$tmp_f;
	    if($verbose_flag==1){
		print $argstring."\n";
	    }
	    $retcode=system($argstring);
	    if($retcode ne 0){
		print STDERR "\nFailure: $argstring\n";
		exit -1;
	    }
	}
    }
    if($chaining_number>1){
	if($plotali_flag==0){
	    my $file_extension2="ccn.ccn";
	    $synteny_prefix=$fragment_file.".$file_extension2.dat.syn.";
	}
	else{
	    my $file_extension2="ccn.filtered.ccn";
	    $synteny_prefix=$fragment_file.".$file_extension2.dat.syn.";
	}
	foreach $combination (@combinations){
	    @projection_array=split(/x/,$combination);
	    $tmp_f=$synteny_prefix.$combination.".gp";
	    open (SynFile, "<$tmp_f") or die "Can't open $tmp_f \n";
	    my @fsynlines = <SynFile>;
	    close SynFile;	
	    open (SynFile, ">$tmp_f") or die "Can't open $tmp_f \n";
	    foreach $fline (@fsynlines){
		if($fline=~/xlabel/){
		    $fline="set xlabel \"$queries[$projection_array[0]-1]\"\n";
		}
		if($fline=~/ylabel/){
		    $fline="set ylabel \"$queries[$projection_array[1]-1]\"\n";
		}
		print SynFile $fline;
	    }
	    close SynFile;
	    $argstring="gnuplot ".$tmp_f;
	    if($verbose_flag==1){
		print $argstring."\n";
	    }
	    $retcode=system($argstring);
	    if($retcode ne 0){
		print STDERR "\nFailure: $argstring\n";
		exit -1;
	    }
	}
    }
    if($plotali_flag==1){
	$synteny_prefix=$fragment_file.".ccn.filtered.dat.syn.";
	foreach $combination (@combinations){
	    @projection_array=split(/x/,$combination);
	    $tmp_f=$synteny_prefix.$combination.".gp";
	    open (SynFile, "<$tmp_f") or die "Can't open $tmp_f \n";
	    my @fsynlines = <SynFile>;
	    close SynFile;	
	    open (SynFile, ">$tmp_f") or die "Can't open $tmp_f \n";
	    foreach $fline (@fsynlines){
		if($fline=~/xlabel/){
		    $fline="set xlabel \"$queries[$projection_array[0]-1]\"\n";
		}
		if($fline=~/ylabel/){
		    $fline="set ylabel \"$queries[$projection_array[1]-1]\"\n";
		}
		print SynFile $fline;
	    }
	    close SynFile;
	    $argstring="gnuplot ".$tmp_f;
	    if($verbose_flag==1){
		print $argstring."\n";
	    }
	    $retcode=system($argstring);
	    if($retcode ne 0){
		print STDERR "\nFailure: $argstring\n";
		exit -1;
	    }
	}
    }
}
      

# removing extra files

if($forward_flag==1){    
    $command = "generate_suffixes_k_genomes_multimat.pl ".$numofqueries ;
    $result = qx($command);
    
    @suffixes = split(/\s/, $result);
        
    foreach $suffix (@suffixes){
	#print $fragment_file.".p".$suffix."\n";
	if($suffix=~/m/){	
	    if (-e $fragment_file.".p".$suffix){
		$command ="rm -f $fragment_file".".p".$suffix."*"; 
		if($verbose_flag==1){
		    #print $command."\n";
		     print "Removing: ".$fragment_file.".p".$suffix."*"."\n";
		}
		system($command);
	    }
	}
    }

}

exit 0;
##########################End program #############################
sub generate_parameter_file
{
    my @file_sizes= ();
    my @in_genome_files=@_;
    my $numoffiles = scalar @in_genome_files;
    my $fsize;
    my $average_f_size=0;
    my $length=0;
    my $rare=0;
    my $gc_value=1000000000;
    my $filter_length_value=0;
    my $lw_value=1;
    use File::stat;

    for(my $i=0; $i<$numoffiles; $i++){	
	$fsize= stat( $in_genome_files[$i])->size;	
	$average_f_size=$average_f_size+$fsize;
	push(@file_sizes,$fsize);
	
	
    }
    $average_f_size=$average_f_size/$numoffiles;
    

    
    #-----------------------------------------------------------
    # estimation of the minimum fragment length
    # this is done based on the table of my dissertation
    #
    #-----------------------------------------------------------
    
    if($numoffiles==2){
	if($average_f_size<80000){	
	    $length=13; 	    
	}
	elsif(($average_f_size >= 80000)and($average_f_size<100000)){	
	    $length=16; 	    
	}
	elsif(($average_f_size >= 100000)and($average_f_size<200000)){
	    $length=17;	    
	}
	elsif(($average_f_size >= 200000)and($average_f_size<500000)){
	    $length=18;	    
	}
	elsif(($average_f_size >= 500000)and($average_f_size<1000000)){
	    $length=19;	    
	}
	elsif(($average_f_size >= 1000000)and($average_f_size<2000000)){
	    $length=20;	    
	}
	elsif(($average_f_size >= 2000000)and($average_f_size<5000000)){
	    $length=22;	    
	}
	elsif(($average_f_size >= 5000000)and($average_f_size<7000000)){
	    $length=24;	    
	}
	elsif(($average_f_size >= 7000000)and($average_f_size<10000000)){
	    $length=24;	    
	}
	elsif(($average_f_size >= 10000000)and($average_f_size<30000000)){
	    $length=25;	    
	}
	elsif(($average_f_size >= 30000000)and($average_f_size<60000000)){
	    $length=26;	    
	}
	elsif(($average_f_size >= 60000000)and($average_f_size<80000000)){
	    $length=27;	    
	}
	elsif(($average_f_size >= 80000000)and($average_f_size<100000000)){
	    $length=28;	    
	}   
	else{
	    $length=30;	    
	}
	$gc_value=100;
	
    }
    elsif($numoffiles==3){
	if($average_f_size<80000){	
	    $length=9; 	    
	}
	elsif(($average_f_size >= 80000)and ($average_f_size<100000)){	
	    $length=10;
	}
	elsif(($average_f_size >= 100000)and($average_f_size<200000)){
	    $length=11;
	}
	elsif(($average_f_size >= 200000)and($average_f_size<500000)){
	    $length=12;
	}
	elsif(($average_f_size >= 500000)and($average_f_size<1000000)){
	    $length=13;
	}
	elsif(($average_f_size >= 1000000)and($average_f_size<2000000)){
	    $length=14;
	}
	elsif(($average_f_size >= 2000000)and($average_f_size<5000000)){
	    $length=15;
	}
	elsif(($average_f_size >= 5000000)and($average_f_size<7000000)){
	    $length=16;
	}
	elsif(($average_f_size >= 7000000)and($average_f_size<10000000)){
	    $length=18;
	}
	elsif(($average_f_size >= 10000000)and($average_f_size<30000000)){
	    $length=20;
	}
	elsif(($average_f_size >= 30000000)and($average_f_size<60000000)){
	    $length=21;
	}
	elsif(($average_f_size >= 60000000)and($average_f_size<80000000)){
	    $length=22;
	}
	elsif(($average_f_size >= 80000000)and($average_f_size<100000000)){
	    $length=23;
	}   
	else{
	    $length=24;
	}
	$gc_value=150;
    }   
    else{
	if($average_f_size<80000){
	    $length=9;
	}
	elsif(($average_f_size >= 80000)and ($average_f_size<100000)){	
	    $length=10;
	}
	elsif(($average_f_size >= 100000)and($average_f_size<200000)){
	    $length=11;
	}
	elsif(($average_f_size >= 200000)and($average_f_size<500000)){
	    $length=12;
	}
	elsif(($average_f_size >= 500000)and($average_f_size<1000000)){
	    $length=13;
	}
	elsif(($average_f_size >= 1000000)and($average_f_size<2000000)){
	    $length=14;
	}
	elsif(($average_f_size >= 2000000)and($average_f_size<5000000)){
	    $length=15;
	}
	elsif(($average_f_size >= 5000000)and($average_f_size<7000000)){
	    $length=16;
	}
	elsif(($average_f_size >= 7000000)and($average_f_size<10000000)){
	    $length=16;
	}
	elsif(($average_f_size >= 10000000)and($average_f_size<30000000)){
	    $length=17;
	}
	elsif(($average_f_size >= 30000000)and($average_f_size<60000000)){
	    $length=17;
	}
	elsif(($average_f_size >= 60000000)and($average_f_size<80000000)){
	    $length=18;
	}
	elsif(($average_f_size >= 80000000)and($average_f_size<100000000)){
	    $length=19;
	}   
	else{
	    $length=20;
	}
	$gc_value=40*$numoffiles;
    }
    $gc_value=500;

    $lw_value=int(($numoffiles*$gc_value)/(2*$length)+2+0.5);
    $filter_length_value=(2*$length);

    $parameter_file_name="home/Microsoftcoconut/parameters.auto";
    open (FILE, ">$parameter_file_name") or die "Failure: Couldn't open $parameter_file_name \n";
    
    if($fragment_prog eq "multimat"){
	

	if($forward_flag==1){
	    print FILE "FRAGMENT=  -l $length \n";
	}
	else{
	    print FILE "FRAGMENT=  -l $length \n";
	}
	print FILE "CHAINING= -l -chainerformat -lw $lw_value -gc $gc_value -length $filter_length_value \n";
    }
    elsif($fragment_prog eq "gt"){
	if($forward_flag==1){
	    print FILE "FRAGMENT=   -l $length \n";
	}
	else{
	    print FILE "FRAGMENT=  -l $length \n";
	}
	print FILE "CHAINING= -l -chainerformat -lw $lw_value -gc $gc_value -length $filter_length_value \n";	
    }
    else{
	$length=$length-4;
	$rare=5;
	print FILE "FRAGMENT= $length $rare \n";
	print FILE "CHAINING= -l -chainerformat -lw $lw_value -gc $gc_value -length $filter_length_value \n";
    }
    
    if($clustalw_flag==1){
	print FILE "ALIGN= -palindrome \n";	
    }
    
    if($syntenic_flag==1){
	
	$block_gap;
	if($average_f_size<200000){
	    $block_gap=$gc_value;
	}   
	elsif(($average_f_size >= 200000)and($average_f_size<500000)){
	    $block_gap=1000;
	}   
	elsif(($average_f_size >= 500000)and($average_f_size<1000000)){
	    $block_gap=10000;
	} 
	elsif(($average_f_size >= 1000000)and($average_f_size<5000000)){
	    $block_gap=20000;
	} 
	elsif(($average_f_size >= 5000000)and($average_f_size<10000000)){
	    $block_gap=50000;
	} 
	elsif(($average_f_size >= 10000000)and($average_f_size<20000000)){
	    $block_gap=80000;
	} 
	elsif(($average_f_size >= 20000000)and($average_f_size<30000000)){
	    $block_gap=200000;
	} 
	elsif(($average_f_size >= 30000000)and($average_f_size<50000000)){
	    $block_gap=500000;
	} 
	else{
	    $block_gap=1000000;
	}

	#print FILE "CHAINING= -neighbor -chainerformat  -gc $block_gap  -length $filter_length_value \n";
	print FILE "SYNTENY= -overlap1 0.2 -filterrep 0.7\n";
    }

    

    close FILE;
    
    return $parameter_file_name;
}


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
    # Kommentarzeile
    if($infileline=~/\#/){
	next;
    }
    else{
	$linecounter++;
	$infileline=~/\d+\s+(\d+)/;
	
    }
    
}
push(@contig_num, $linecounter); 
close(matchfile);
return @contig_num;
}

sub generate_seqinfo_files
{

    my @queries=@_;
    my @arg   =@_;
    
#    my $in_index_name=$arg[0];
#    my $db_name=$arg[1];
#    my $query_name=$arg[1];

    $command = "generate_suffixes_k_genomes_multimat.pl ".$numofqueries ;
    $result = qx($command);    
    @suffixes = split(/\s/, $result);
    

    my $prog_seqinfo=$fragment_generation_path."/vseqinfo";
       
    $argstring=$prog_seqinfo." ".$index_name." > ".$queries[0].".seqinfo";
    if($verbose_flag==1){
	print "Running ".$argstring."\n";
    }
    my $retcode = system($argstring);
    if($retcode ne 0){
	print STDERR "\nFailure: $argstring\n";
	exit -1;
    }

#    for($i=1;$i<$numofqueries;$i++){
#	# generate index
#	
#	my $seq_name = qx/basename $queries[$i]/;
#	chop $seq_name;
#	my $in_index_name=$dir_name."/$seq_name".".index";
#	if($prefix_flag==1){
#	    my $prefix_dir = qx/dirname $prefix/;
#	    chop $prefix_dir;
#	    $in_index_name=$prefix_dir."/$seq_name.index";
#	}
#	$argstring=$prog_indexing." -tis -ois -dna  -v -indexname ".$in_index_name." -db ".$queries[$i];
#	if($verbose_flag==1){
#	    print "Running ".$argstring."\n";
#	}
#	my $retcode = system($argstring);
#	if($retcode ne 0){
#	    print STDERR "\nFailure: $argstring\n";
#	    exit -1;
#	}
#	$argstring=$prog_seqinfo." ".$in_index_name." > ".$queries[$i].".seqinfo";	
#	if($verbose_flag==1){
#	    print "Running ".$argstring."\n";
#	}
#	my $retcode = system($argstring);
#	if($retcode ne 0){
#	    print STDERR "\nFailure: $argstring\n";
#	    exit -1;
#	}
#    }
   
    
    foreach $suffix (@suffixes){
	
	# create file
	$outfile=$fragment_file.".seqinfo.p".$suffix;
	$argstring="cp ".$queries[0].".seqinfo"."  $outfile";
	my $retcode = system($argstring);
	if($retcode ne 0){
	    print STDERR "\nFailure: $argstring\n";
	    exit -1;
	}
	#print $argstring." \n";

	open(out_file_handle,">>".$outfile) 
	    or die "Oops! can't open seqinfo input file \n";  
	$tmp_suf="p".$suffix;

	my @dirc=split("",$tmp_suf);
	
	
	for($j=1;$j<$numofqueries;$j++){
	    print out_file_handle "#--- \n";
	    
	    if($dirc[$j]eq "p"){
		
		# add file
		open(seqinfo_file,$queries[$j].".seqinfo")
		    or die "Oops! can't open seq input file "."\n";  
		my @in=<seqinfo_file>;
		print out_file_handle @in;
		close(seqinfo_file);
	    }
	    else{
		# add reverse file
		
		open(seqinfo_file,$queries[$j].".seqinfo")
		    or die "Oops! can't open seq input file "."\n";  
		my @in=<seqinfo_file>;
		my @rev=();
		@rev= reverse(@in);				
		print out_file_handle @rev;
		close(seqinfo_file);
	    }
	}
	close(out_file_handle);
    }
    
}

# removing extra files
sub cleanup{
if($forward_flag==1){    
    $command = "generate_suffixes_k_genomes_multimat.pl ".$numofqueries ;
    $result = qx($command);
    
    @suffixes = split(/\s/, $result);
        
    foreach $suffix (@suffixes){
	#print $fragment_file.".p".$suffix."\n";
	if($suffix=~/m/){
	
	    if (-e $fragment_file.".p".$suffix){
		$command ="rm -f $fragment_file".".p".$suffix."*"; 
		if($verbose_flag==1){
		    #print $command."\n";
		     print "Removing: ".$fragment_file.".p".$suffix."*"."\n";
		}
		system($command);
	    }
	}
    }


}
}


sub is_number
{
    my @input=@_;
    my $x=$input[0];
    #print "First Arg $x \n";
    $x =~ s/\.//;
    #print "After Sub: $x \n";
    if($x =~/(\D)/){
	return 0;
    }
    else{
	#print STDERR "Correct: $x is a number\n";	
	if ($x =~/^-?\d*\.*\d+/){
	    #print STDERR "Correct: $x is a number\n";	    
	   return 1;
	} 
	else{
	    return 0;
	}
    }
}
