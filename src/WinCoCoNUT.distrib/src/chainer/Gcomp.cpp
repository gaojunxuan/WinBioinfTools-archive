//////////////////////////////////////////////////////////////////////
// Gcomp.cpp: Implementatiion of the main function.
//
//////////////////////////////////////////////////////////////////////


/*
  Copyright by Mohamed I. Abouelhoda (C) 2003
  =====================================
  You may use, copy and distribute this file freely as long as you
   - do not change the file,
   - leave this copyright notice in the file,
   - do not make any profit with the distribution of this file
   - give credit where credit is due
  You are not allowed to copy or distribute this file otherwise
  The commercial usage and distribution of this file is prohibited
  Please report bugs and suggestions to <mibrahim@informatik.uni-ulm.de>
*/

#include <sys/times.h>
#include <unistd.h> 
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "kdChainer.h"
#include "ProcessFragmentFile.h"
#include "ProcessOptions.h"
#include "ContigSorting.h"
#include "myglobaldef.h"
//necessary for testing
#define SUCCESS 0
#define FAIL 1

int main(long argc, char **argv)
 {


   int res;
   ProcessOptions OptionsObj;
   
   char tempstr[1000];
   
   if(argc<3){
       printf("*//////////////////////////////////////////////////*\n");
       printf("*          This is CHAINER version 3.0             *\n");
       printf("*                 %d-bits version                  *\n",VERSION);
       printf("*   Copyright by Mohamed I. Abouelhoda  (C) 2004   *\n");
       printf("*//////////////////////////////////////////////////*\n");
       printf("Usage:\n");
       printf("chainer  [options]  fragment_file\n");
       printf("\nOptions:\n");
       printf("-l:               generate local chains \n");
       printf("-g:               generate a global chain without gap cost\n");
       printf("-gg:              generate a global chain including gap costs\n");
       printf("-draft:           compare a draft genome to a finished genome\n");
       printf("-est:             map cDNA/EST database to a genomic sequence\n");
       printf("-r:               generate local chains for repeats\n");
       printf("-des:             description file for the contigs or the ESTs\n");
       printf("-s:               minimum score of the reported chains (default =0) \n"); 
       printf("-length:          minimum average chain length \n");
       printf("-perc:            minimum percentage score w.r.t. highest score \n");
       printf("-d:               chain size (min number of fragments in a chain) \n");
       printf("-lw:              multiplies the fragment weight by a given value\n"); 
       printf("-gc:              maximum gap allowed between fragments\n");       
       printf("-gcfile:          file specifying max. gap allowed between fragments \n");
       printf("                  in each genome (only with options -r and -rc)\n");       
       printf("-overlap:         maximum amount of overlap between fragments\n");  
       printf("-side:            sequence(s) in which overlap is allowed\n");   
       printf("-coverage:        coverage ratio of the EST\n");
       printf("-chainerformat:   output to files in chainer-format \n");
       printf("-stdout:          output to standard output in table-format \n");
       printf("-cluster:         output cluster of local chains \n");
       printf("-neighbor:        output cluster of neighboring fragments \n");
       printf("-rc:              output cluster of local chains for repeats \n");
       printf("...               For more information, see the manual\n");
       return(FAIL);
   }
 

   ProcessFragmentFile InfoObj;//=new ProcessFragmentFile;
   
   // processing options

   int Operation=OptionsObj.ParseOptions(argc,argv);

   // generating statistics file
   if(Operation>0){
       res=InfoObj.GetFragmentFileStatisticsFormat2(argv[OptionsObj.fileindexarray[0]]);
       if(res<0)return(FAIL);
   }
   else{      
       printf("Invalid options: No action is given in the option \n");
       return(FAIL);
   }
   
   kdChainer* obj=  new kdChainer(InfoObj.genome_size_array);
   
   // Processing the options be calling the corresponding function
   /* The option number is given according to the time point at which this option
      is added to the program
   */

   if(Operation==1){/*Operation Global Chaining */
       obj->set_multiplying_factor((double)OptionsObj.optionsarray[1]);
       obj->set_stdout_flag((int)OptionsObj.optionsarray[2]);
       obj->set_chainer_format_flag((int)OptionsObj.optionsarray[3]);
       
       if(((int)OptionsObj.optionsarray[4])<(int)((double)InfoObj.min_fragment_length)){
	   obj->set_overlapping_value((int)OptionsObj.optionsarray[4]);
       }
       else{
	   printf("Maximum overlapping value is > the minimum sequence length\n");
	   delete obj;
	   return (FAIL);
       }
       if(((int)OptionsObj.optionsarray[5])<=InfoObj.no_of_genomes){
	   obj->set_overlapping_side((int)OptionsObj.optionsarray[5]);
       }
       else{
	   printf("Overlapping side is > no. of genomes\n");
	   delete obj;
	   return (FAIL);
       }

       res=obj->KdtreebasedModgetchain(&InfoObj,0,0);              
       delete obj;
       
   }
   else if(Operation==2){ /*Operation Global Chaining Include Extra Weight file-obselete option*/
       obj->KdtreebasedModgetchain(&InfoObj,1,0); 
       delete obj;
   }
   else if(Operation==3){ /*Operation Global Chaining With gap costs */
       obj->set_multiplying_factor((double)OptionsObj.optionsarray[1]);
       obj->set_stdout_flag((int)OptionsObj.optionsarray[2]);
       obj->set_chainer_format_flag((int)OptionsObj.optionsarray[3]);       
       if(((int)OptionsObj.optionsarray[4])<(int)((double)InfoObj.min_fragment_length)){
	   obj->set_overlapping_value((int)OptionsObj.optionsarray[4]);
       }
       else{
	   printf("Maximum overlapping value is > the minimum sequence length\n");
	   delete obj;
	   return (FAIL);
       }
       if(((int)OptionsObj.optionsarray[5])<=InfoObj.no_of_genomes){
	   obj->set_overlapping_side((int)OptionsObj.optionsarray[5]);
       }
       else{
	   printf("Overlapping side is > no. of genomes\n");
	   delete obj;
	   return (FAIL);
       }
       res=obj->GappedKdtreebasedModgetchain(&InfoObj,0,0);
       delete obj;
   }
   else if(Operation==5){/*Operation Local Chaining*/
       obj->set_chain_depth((long)OptionsObj.optionsarray[2]);
       obj->set_chain_score((double)OptionsObj.optionsarray[0]);
       obj->set_stdout_flag((int)OptionsObj.optionsarray[7]);
       obj->set_chainer_format_flag((int)OptionsObj.optionsarray[8]);
       obj->set_cluster_file_flag((int)OptionsObj.optionsarray[9]);
       obj->set_chain_average_length(OptionsObj.optionsarray[10]);

       //printf("\n%f\n",OptionsObj.optionsarray[5]);
       obj->set_multiplying_factor(OptionsObj.optionsarray[5]);
       if(OptionsObj.optionsarray[1]>=(double)0){
		 obj->set_percentage_threshold((double)OptionsObj.optionsarray[1]);
       }
       if(OptionsObj.optionsarray[6]>0){
	   obj->set_gap_threshold_flag();
	   obj->set_gap_threshold((long)OptionsObj.optionsarray[6]);
	   
       }

       if(((int)OptionsObj.optionsarray[11])<(int)(InfoObj.min_fragment_length)){
	   obj->set_overlapping_value((int)OptionsObj.optionsarray[11]);
       }
       else{
	   printf("Maximum overlapping value is > the minimum sequence length\n");
	   delete obj;
	   return (FAIL);
       }
       if(((int)OptionsObj.optionsarray[12])<=InfoObj.no_of_genomes){
	   obj->set_overlapping_side((int)OptionsObj.optionsarray[12]);
       }
       else{
	   printf("Overlapping side is > no. of genomes\n");
	   delete obj;
	   return (FAIL);
       }
       res=obj->KdtreeSW(&InfoObj,0,0,0);
       delete obj;
   }
   else if(Operation==9){ // operation for handling repeats
       obj->set_chain_depth((long)OptionsObj.optionsarray[2]);
       obj->set_chain_score((double)OptionsObj.optionsarray[0]);
       obj->set_stdout_flag((int)OptionsObj.optionsarray[7]);
       obj->set_chainer_format_flag((int)OptionsObj.optionsarray[8]);
       obj->set_cluster_file_flag((int)OptionsObj.optionsarray[9]);
       obj->set_chain_average_length(OptionsObj.optionsarray[10]);
       obj->set_detect_repeat_flag();
       //printf("\n%f\n",OptionsObj.optionsarray[5]);
       obj->set_multiplying_factor(OptionsObj.optionsarray[5]);
       if(OptionsObj.optionsarray[1]>(double)0){
		 obj->set_percentage_threshold((double)OptionsObj.optionsarray[1]);
       }
       if(OptionsObj.optionsarray[6]>0){
	   obj->set_gap_threshold_flag();
	   obj->set_gap_threshold((long)OptionsObj.optionsarray[6]);
	   
       }

       if(((int)OptionsObj.optionsarray[11])<(int)(InfoObj.min_fragment_length)){
	   obj->set_overlapping_value((int)OptionsObj.optionsarray[11]);
       }
       else{
	   printf("Maximum overlapping value is > the minimum sequence length\n");
	   delete obj;
	   return (FAIL);
       }
       if(((int)OptionsObj.optionsarray[12])<=InfoObj.no_of_genomes){
	   obj->set_overlapping_side((int)OptionsObj.optionsarray[12]);
       }
       else{
	   printf("Overlapping side is > no. of genomes\n");
	   delete obj;
	   return (FAIL);
       }
       res=obj->KdtreeSW(&InfoObj,0,0,0);
       delete obj;
   }
   
   else if(Operation==11){ // operation for clustering repeats
       obj->set_chain_depth((long)OptionsObj.optionsarray[2]);
       obj->set_chain_score((double)OptionsObj.optionsarray[0]);
       obj->set_stdout_flag((int)OptionsObj.optionsarray[7]);
       obj->set_chainer_format_flag((int)OptionsObj.optionsarray[8]);
       obj->set_cluster_file_flag((int)OptionsObj.optionsarray[9]);
       obj->set_chain_average_length(OptionsObj.optionsarray[10]);
       obj->set_detect_repeat_cluster_flag();
       //printf("\n%f\n",OptionsObj.optionsarray[5]);
       obj->set_multiplying_factor(OptionsObj.optionsarray[5]);
       if(OptionsObj.optionsarray[13]==1){	   
	   obj->set_gcfile(argv[OptionsObj.fileindexarray[2]]);
       }
       
       if(OptionsObj.optionsarray[1]>(double)0){
		 obj->set_percentage_threshold((double)OptionsObj.optionsarray[1]);
       }
       if(OptionsObj.optionsarray[6]>0){
	   obj->set_gap_threshold_flag();
	   obj->set_gap_threshold((long)OptionsObj.optionsarray[6]);
	   
       }

       if(((int)OptionsObj.optionsarray[11])<(int)(InfoObj.min_fragment_length)){
	   obj->set_overlapping_value((int)OptionsObj.optionsarray[11]);
       }
       else{
	   printf("Maximum overlapping value is > the minimum sequence length\n");
	   delete obj;
	   return (FAIL);
       }
       if(((int)OptionsObj.optionsarray[12])<=InfoObj.no_of_genomes){
	   obj->set_overlapping_side((int)OptionsObj.optionsarray[12]);
       }
       else{
	   printf("Overlapping side is > no. of genomes\n");
	   delete obj;
	   return (FAIL);
       }
       //res=obj->KdtreeSW(&InfoObj,0,0,0);
       res=obj->ContigESTChainer(&InfoObj,0,0,0,NULL,1);
       delete obj;
   }
   /////// end of repeat cluster

   // neighboring option
   else if(Operation==10){/*Operation Neighboring Clustering*/       
       
       obj->set_chain_depth((long)OptionsObj.optionsarray[2]);
       obj->set_chain_score((double)OptionsObj.optionsarray[0]);
       obj->set_stdout_flag((int)OptionsObj.optionsarray[7]);
       obj->set_chainer_format_flag((int)OptionsObj.optionsarray[8]);
       obj->set_cluster_file_flag((int)OptionsObj.optionsarray[9]);
       obj->set_chain_average_length(OptionsObj.optionsarray[10]);

       //printf("\n%f\n",OptionsObj.optionsarray[5]);
       obj->set_multiplying_factor(OptionsObj.optionsarray[5]);
       if(OptionsObj.optionsarray[13]==1){	   
	   obj->set_gcfile(argv[OptionsObj.fileindexarray[2]]);
       }

       //obj->set_cluster_flag();
       //obj->set_chain_score((double)0);
       //obj->set_multiplying_factor(1); // not needed for clustering

       if(OptionsObj.optionsarray[1]>(double)0){
	   	 obj->set_percentage_threshold((double)OptionsObj.optionsarray[1]);
       }
       if(OptionsObj.optionsarray[6]>0){
	   obj->set_gap_threshold_flag();
	   obj->set_gap_threshold((long)OptionsObj.optionsarray[6]);
	   
       }
    

       if(((int)OptionsObj.optionsarray[11])<(int)(InfoObj.min_fragment_length)){
	   obj->set_overlapping_value((int)OptionsObj.optionsarray[11]);
       }
       else{
	   printf("Maximum overlapping value is > the minimum sequence length\n");
	   delete obj;
	   return (FAIL);
       }
       if(((int)OptionsObj.optionsarray[12])<=InfoObj.no_of_genomes){
	   obj->set_overlapping_side((int)OptionsObj.optionsarray[12]);
       }
       else{
	   printf("Overlapping side is > no. of genomes\n");
	   delete obj;
	   return (FAIL);
       }
       if(OptionsObj.optionsarray[7]&&OptionsObj.optionsarray[8]){
	   // contradicted options
	   printf("Contradicted options: -stdout and -chainerformat\n");
	   delete obj;
	   return (FAIL);
       }
       
//       res=obj->KdtreeSW(&InfoObj,0,0,0);
       res=obj->ContigESTChainer(&InfoObj,0,0,0,NULL,1);
       delete obj;
   
   }
   // end of neighboring option
   
   else if(Operation==7){/*Operation finished-draft genome comparison*/
       obj->set_chain_depth((long)OptionsObj.optionsarray[2]);
       obj->set_chain_score(OptionsObj.optionsarray[0]);
       obj->set_multiplying_factor((double)OptionsObj.optionsarray[5]);
       obj->set_stdout_flag((int)OptionsObj.optionsarray[1]);
       obj->set_chainer_format_flag((int)OptionsObj.optionsarray[8]);
       obj->set_absolute_flag((int)OptionsObj.optionsarray[9]);
       obj->set_cluster_file_flag((int)OptionsObj.optionsarray[10]);
       obj->set_chain_average_length(OptionsObj.optionsarray[11]);

       
       if(((int)OptionsObj.optionsarray[12])<(int)((double)InfoObj.min_fragment_length)){
	   obj->set_overlapping_value((int)OptionsObj.optionsarray[12]);
       }
       else{
	   printf("Maximum overlapping value is > the minimum sequence length\n");
	   delete obj;
	   return (FAIL);
       }
       if(((int)OptionsObj.optionsarray[13])<=InfoObj.no_of_genomes){
	   obj->set_overlapping_side((int)OptionsObj.optionsarray[13]);
       }
       else{
	   printf("Overlapping side is > no. of genomes\n");
	   delete obj;
	   return (FAIL);
       }

       if(OptionsObj.optionsarray[8]&&OptionsObj.optionsarray[9]){
	   // contradicted options
	   printf("Contradicted options: -chainerformat and -absolute\n");
	   delete obj;
	   return (FAIL);
       }
       
       if(OptionsObj.optionsarray[1]&&OptionsObj.optionsarray[8]){
	   // contradicted options
	   printf("Contradicted options: -stdout and -chainerformat\n");
	   delete obj;
	   return (FAIL);
       }
       
       if(OptionsObj.optionsarray[4]>=(double)0){
		 obj->set_percentage_threshold((double)OptionsObj.optionsarray[4]);
       }
       if(OptionsObj.optionsarray[6]>(double)0){
	   obj->set_gap_threshold_flag();
	   obj->set_gap_threshold((long)OptionsObj.optionsarray[6]);     
       }
       
       if(OptionsObj.optionsarray[3]==0){ // not include weight file  // obsolete
	   if(OptionsObj.fileindexarray[1]!=0){	 
	       res=obj->ContigESTChainer(&InfoObj,0,0,0,argv[OptionsObj.fileindexarray[1]],0);
	   }
	   else{
	       res=obj->ContigESTChainer(&InfoObj,1,0,0,NULL,0);	       
	   }
       }
       else{
	   if(OptionsObj.fileindexarray[1]!=0){	   
	       
	       res=obj->ContigESTChainer(&InfoObj,1,0,0,argv[OptionsObj.fileindexarray[1]],0);
	   }
	   else{
	       
	       res=obj->ContigESTChainer(&InfoObj,1,0,0,NULL,0);
	       
	   }
       }
       delete obj;
       
       
       if(res<0)return(FAIL);
       
       // Reporting chains in reference coordiantes
       
       if((int)OptionsObj.optionsarray[8]){
	   ContigSorting contigobj;
	   if(OptionsObj.fileindexarray[1]!=0){
	       contigobj.read_kd_contigs(argv[OptionsObj.fileindexarray[1]],0);
	   }
	   else{
	       contigobj.read_kd_contigs(NULL,0);
	   }
	   
	   sprintf(tempstr,"%s.ccn",argv[OptionsObj.fileindexarray[0]]);
	   //contigobj.kd_read_chain_file(tempstr,InfoObj.no_of_genomes+1);
	   contigobj.kd_read_chain_file_kdraftgenome(tempstr,InfoObj.no_of_genomes+1);
	   sprintf(tempstr,"%s.chn",argv[OptionsObj.fileindexarray[0]]);
	   //contigobj.kd_read_chain_file(tempstr,InfoObj.no_of_genomes+1);
	   contigobj.kd_read_chain_file_kdraftgenome(tempstr,InfoObj.no_of_genomes+1);
       }
   }
   else if(Operation==8){/*Operation cDNA/EST placement*/
       
       obj->set_chain_depth((long)OptionsObj.optionsarray[2]);
       obj->set_chain_score(OptionsObj.optionsarray[0]);
       obj->set_multiplying_factor((double)OptionsObj.optionsarray[5]);
       obj->set_stdout_flag((int)OptionsObj.optionsarray[1]);
       obj->set_chainer_format_flag((int)OptionsObj.optionsarray[9]);
       obj->set_absolute_flag((int)OptionsObj.optionsarray[10]);
       obj->set_cluster_file_flag((int)OptionsObj.optionsarray[11]);

       if(((int)OptionsObj.optionsarray[12])<(int)((double)InfoObj.min_fragment_length)){
	   obj->set_overlapping_value((int)OptionsObj.optionsarray[12]);
       }
       else{
	   printf("Maximum overlapping value is > the minimum sequence length\n");
	   delete obj;
	   return (FAIL);
       }
       
       if(((int)OptionsObj.optionsarray[13])<=InfoObj.no_of_genomes){
	   obj->set_overlapping_side((int)OptionsObj.optionsarray[13]);
       }
       else{
	   printf("Overlapping side is > no. of genomes\n");
	   delete obj;
	   return (FAIL);
       }

       if(OptionsObj.optionsarray[9]&&OptionsObj.optionsarray[10]){
	   // contradicted options
	   printf("Contradicted options: -chainerformat and -absolute\n");
	   delete obj;
	   return (FAIL);
       }  
       if(OptionsObj.optionsarray[1]&&OptionsObj.optionsarray[9]){
	   // contradicted options
	   printf("Contradicted options: -stdout and -chainerformat\n");
	   delete obj;
	   return (FAIL);
       }
              
       //Gaps for EST
       if(OptionsObj.optionsarray[8]>(double)0){
	   obj->set_gap_threshold_flag();
	   obj->set_gap_threshold((long)OptionsObj.optionsarray[8]);     
	   
       }

       if(OptionsObj.optionsarray[6]>=(double)0){
	   obj->set_percentage_threshold((double)OptionsObj.optionsarray[6]);
       }
       
       if(OptionsObj.optionsarray[3]==0){ // not include weight file 
	   if(OptionsObj.fileindexarray[1]!=0){
	       res=obj->ContigESTChainer(&InfoObj,0,0,0,argv[OptionsObj.fileindexarray[1]],1);
	   }
	   else{
	       // no description file is given, the database is dealt with as one sequence.
	        res=obj->ContigESTChainer(&InfoObj,0,0,0,NULL,1);
		
	   }
       }
       else{
	   
	   //  res=obj->ContigESTChainer(&InfoObj,1,0,0,argv[OptionsObj.fileindexarray[1]],1);
       }
       delete obj;

       if(res<0)return(FAIL);
       // Reporting chains in reference coordinates
       if((int)OptionsObj.optionsarray[9]){
	   ContigSorting contigobj;
	   if(OptionsObj.fileindexarray[1]!=0){
	       contigobj.read_contigs(argv[OptionsObj.fileindexarray[1]],0);
	   }
	   else{
	       contigobj.read_contigs(NULL,0);
	   }
	   sprintf(tempstr,"%s.ccn",argv[OptionsObj.fileindexarray[0]]);
	   contigobj.kd_read_chain_file(tempstr,InfoObj.no_of_genomes+1);
	   sprintf(tempstr,"%s.chn",argv[OptionsObj.fileindexarray[0]]);
	   contigobj.kd_read_chain_file(tempstr,InfoObj.no_of_genomes+1);
	   
       }
   }
   else{
       delete obj;
       printf("Invalid options\n");
       return(FAIL);
   }
   if(res<0)return(FAIL);
   //printf("SUCCESS SUCCESS %d",SUCCESS);
   return (SUCCESS);
      
 }
