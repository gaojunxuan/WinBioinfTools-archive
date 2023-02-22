#!/usr/bin/perl
# use strict;
# use warnings;


my $numofargs = scalar @ARGV;

my $task_comp_multiple_flag=0;
my $task_comp_pairwise_flag=0;
my $task_map_cdna_flag=0;
my $task_repeat_flag=0;
my $verbose_flag=0;

if($numofargs< 1){
    print STDERR "\nUsage: coconut.pl -task_name arguments\n";
    print STDERR "\n";
    print STDERR "           -pairwise --> compare two finished or draft genomes\n";
    print STDERR "To view arguments for each task, run coconut.pl with task_name\n";
    print STDERR "Example: \n > coconut.pl -multiple \n\n";
    exit -1;
}

if($ARGV[0] eq "-multiple"){
    $task_comp_multiple_flag=1;
}
elsif ($ARGV[0] eq "-pairwise"){
    $task_comp_pairwise_flag=1;
}
else{
    print STDERR "Error: unknown task name\n";
    exit -1;
}

my $argstr="  ";

for(my $i=1; $i<$numofargs; $i++)
{    
    $argstr=$argstr." ".$ARGV[$i];
    if($ARGV[$i] eq "-v"){
	$verbose_flag=1;
    }
    #print $ARGV[$i]."\n";
}

$argstr=" -cf home/Microsoftcoconut/config ".$argstr;

print "\n";
if($task_comp_multiple_flag==1){
    $ENV{"PATH"} = "finalscripts/comp_finished:".$ENV{"PATH"};
    $argstr="comp_multiple_genome.pl ".$argstr;
}
elsif($task_comp_pairwise_flag==1){
    $ENV{"PATH"} = "finalscripts/comp_pairwise:".$ENV{"PATH"};
    $argstr="perl Home/MicroSoftCoCoNUT/finalscripts/comp_pairwise/comp_pairwise_draft.pl ".$argstr;
}


if($verbose_flag==0){
  $argstr=$argstr." -v  > logfile ";  
}
my $ret_val=system($argstr);
if($ret_val != 0){
    print STDERR "Error running: $argstr\n";
    exit -1;
}
print "\n";

$argstr="rm -f  *-match";
system($argstr);
$argstr="rm -f *+match";
system($argstr);

exit 0;



#for(my $i=0; $i<$numofargs; $i++)
#{    
  
#    print $ARGV[$i]."\n";
#}
#print "No. of arg. $numofargs\n";
#print $0;
#print "\n";

#$original_path=$ENV{PATH};
#print $original_path;
#print "\n";
#print "new set path\n";
#print $ENV{'PATH'};
#$ENV{PATH}="finalscripts/comp_finished:$ENV{PATH}";
#$argstr="PATH=finalscripts/comp_finished:$PATH; export PATH";
#$argstr="export PATH=$PATH:finalscripts/comp_finished";
#system("printenv PATH");
#qx(export PATH=$PATH:finalscripts/comp_finished);
#exec("export PATH=$PATH:finalscripts/comp_finished");
#system("printenv PATH");
#system($argstr);
#print $argstr."\n";
#$retval=system($argstr);
#print "Ret_Val: $retval \n";

#$new_path=$ENV{PATH};
#print $new_path."\n";
