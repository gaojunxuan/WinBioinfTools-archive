#include<stdio.h>
#include<string.h>
#include<stdlib.h>
#include "myglobaldef.h"

#define max(a, b)  ((a)>=(b) ? (a) : (b))

#define FAIL 1
#define SUCCESS 0
int main(int argc,char *argv[]){
  
    //char filename[1000];
  char line[MaxFragmentLine+50];
  char tmp_line[MaxFragmentLine+50];
//  char* line2;
  long length1=0;
//  long  seqno;
  long pos1=0;
//  unsigned char teampchar;
  long length2=0;
//  long seqno2;
  long pos2=0;
//  int type; // match type
//  int filetype=0;
  char  count=0;
//  double  evalue;
//  double doubleval;
//  unsigned int intval;
  FILE* fptr=NULL;
  FILE* fptr2=NULL;
  FILE* fptr3=NULL;
  FILE* logfptr=NULL;
//  int i;
  int k=0;
//  double slope;
  char* token;
  char tempstr[5000];
  char tempstr1[5000];
  // long  no_of_genomes=0;
//  int flag=0;
//  long kk=0;
  long offset=0;
  long genomesizearray[128];
  // long maxfragmentlength=0;
  long maxfragmentlength1=0;
  // long maxfragmentlength2=0;
  //int open_forward_file_flag=1;
//  int open_reverse_file_flag=1;
  // int j;
  long db_length=0;
  long query_length=0;
  long search_length=0;

  if(argc!=2)
    {
    printf("Invalid options\n");
    printf("Usage: rep2chainer filename (in memspe format) \n");    
    printf("Output:\n ");
    printf("  -One file: 'filename.pm' in case of multimat input file \n");    
    return(FAIL);
  } 
  
  
  fptr=FOPEN(argv[1],"r");
  if(fptr==NULL){
      printf("FILE not found \n");
      return(FAIL);
  }
  
  sprintf(tempstr,"%s.pm",argv[1]);
  fptr3=FOPEN(tempstr,"w"); 
  if(fptr3==NULL){
      printf("unable to open file \n"); return(FAIL);
  }
  if(fptr3!=NULL)fprintf(fptr3,">CHA %d\n",(int)2); 


  if(1){
      // compute statistics
      for(k=0;k<128;k++){
	  genomesizearray[k]=0;
      }
      while(feof(fptr)==0){
	  strcpy(line,"");
	  line[0]=0;
	  fgets(line,MaxFragmentLine,fptr);
	  strcpy(tmp_line,line);
	  token = strtok (line, " " );  
	  if(token==NULL)continue;
	  if (strcmp(token,"#")==0) {
	      token = strtok ( NULL, " ");	      	      
	      if(strcmp(token,"output")==0){
		  break;
	      }
	      else if(strcmp(token,"searchlength")==0){
		  token = strtok ( NULL, "=");		      
		  search_length=atol(token);
		  printf("Min. search len. %lu\n",search_length);
	      }
	      else{
		  //printf("%s\n",token);
		  token = strtok (token, "=");
		  //printf("%s\n",token);
		  if((strcmp(token,"databaselength")==0)||(strcmp(token,"querylength")==0)){
		      //printf("%s\n",token);
		      if(strcmp(token,"databaselength")==0){
			  //printf("%s\n",token);
			  token = strtok ( NULL, "=");		      
			  db_length=atol(token);
			  query_length=db_length;
			  printf("Length: %lu\n",db_length);
			  //break;
			  token = strtok (NULL, " ");
			  if(token==NULL) continue;
			  if(strcmp(token,"(including")==0){
			      token = strtok ( NULL, " ");
			      //  printf("Sep %s\n",token);  
			      //db_length=db_length-atol(token);
			  }
			  //break;
		      }
		  }
	      }
	      
	  }
	  
	  
	  	      	 	  
      }
      //printf("db_len %lu query_len %lu \n",db_length,query_length);

      if((db_length==0)){
	  printf("Error, unable to locate genome lengths, use verbose mode\n");
	  fclose(fptr);
	  return 0;
      }
      genomesizearray[0]=db_length;
      genomesizearray[1]=query_length;
      offset=genomesizearray[1];

      while(feof(fptr)==0){	  
	  strcpy(line,"");
	  line[0]=0;
	  fgets(line,MaxFragmentLine,fptr);
	  count=0;
	  token = strtok (line, " " );  
	  if(token==NULL)continue;
	  if (strcmp(token,"#")==0) {
	      //printf ( "The tokens first are " );
	      //printf ( "%s+", token );
	      continue;
	  }
	  do{
	      if(count==0){      
		  //fprintf(fptr2, "%d ", atoi(token));
		  length1=atoi(token);
		  maxfragmentlength1=max(maxfragmentlength1,length1);
		  //printf("len %d\n",length1);
	      }
	      token = strtok ( NULL, " ");
	      count++;
	      
	      if(count==1){
		  pos1=atoi(token);
	      }
	      if(count==2){		  
		  pos2=(-1)*atoi(token);   		  
	      }		      
	  }while(token);
	  { //palindrome
	      if(fptr3!=NULL){
		  if(pos1==pos2){	
		      //printf("hallooooooooooooo %d %d\n",pos1,pos2);
		      length1=length1/2;
		      if(length1>=search_length){
			  pos2=pos1+length1;		      
			  fprintf(fptr3,"#%.2f\n",(float)length1);  
			  fprintf(fptr3,"[%lu,%lu] [%lu,%lu]\n",pos1,pos1+length1-1,offset-(pos2+length1-1)-1,offset-pos2-1);
		      }
		  }
		  else if(pos1<pos2){
		      fprintf(fptr3,"#%.2f\n",(float)length1);  
		      fprintf(fptr3,"[%lu,%lu] [%lu,%lu]\n",pos1,pos1+length1-1,offset-(pos2+length1-1)-1,offset-pos2-1);
		  }
	      }
	      genomesizearray[0]=max(genomesizearray[0],pos1+length1);
	      genomesizearray[1]=max(genomesizearray[1],pos2+length2);	   
	  }
      } // endwhile
      fclose(fptr);
      
      
      
      sprintf(tempstr,"%s.info",argv[1]);

      logfptr=FOPEN(tempstr,"w");            

      fprintf(logfptr,"%d %lu %lu\n",2,db_length,query_length);  
      sprintf(tempstr1,"%s.pp",argv[1]);
      fprintf(logfptr,"%s\n",tempstr1); 
      sprintf(tempstr1,"%s.pm",argv[1]);
      fprintf(logfptr,"%s\n",tempstr1); 
      fclose(logfptr);
      
       
      if(fptr3!=NULL)fclose(fptr3);
  
      return(SUCCESS);
      
        
       
      while(feof(fptr)==0){
	  strcpy(line,"");
	  line[0]=0;
	  fgets(line,MaxFragmentLine,fptr);
	  count=0;
	  token = strtok (line, " " );  
	  if(token==NULL)continue;
	  if (strcmp(token,"#")==0) {
	      //printf ( "The tokens first are " );
	      //printf ( "%s+", token );
	      continue;
	  }
	  do{
	      if(count==0){      
		  //fprintf(fptr2, "%d ", atoi(token));
		  length1=atoi(token);
	      }
	      token = strtok ( NULL, " ");
	      count++;
	      //      if(count==2){ // used in case of palindroms
	      if(count==1){
		  pos1=atoi(token);
		  //fprintf(fptr2, "%d ", atoi(token));
		  //printf("%d %d %d\n", length,(pos1),(pos2)); 
	      }
	      if(count==2){
		  
		  pos2=(-1)*atoi(token);
		  //fprintf(fptr2, "%d\n", atoi(token));
		  //printf("%d %d %d\n", length,(pos1),(pos2));    
	      }
	      /*
	      //if(count==3){  // used in case of pal
	      if(count==2){
		  if(strcmp(token,"D")==0)
		      type=1;
		  //fprintf(fptr3, "%d\n", 0);
		  else type=0;
		  //fprintf(fptr3, "%d\n", 1);
		  //printf("%d %d %d\n", length,(pos1),(pos2)); 
	      }
	      // if(count==6){  // used in case of pal
	      if(count==3){	  
		  length2=atoi(token);  
	      }
	      if(count==4){
		  
		  pos2=atoi(token);
		  //fprintf(fptr2, "%d\n", atoi(token));
		  //printf("%d %d %d\n", length,(pos1),(pos2));    
	      }
	      */
	  }while(token);
	  if(0){
	      //fprintf(fptr2,"%u %u %u\n",length1,pos1,pos2);
	      // the CHAINER format
	      if(fptr2!=NULL)fprintf(fptr2,"#%.2f\n",(float)length1);
	      if(fptr2!=NULL)fprintf(fptr2,"[%lu,%lu] [%lu,%lu]\n",pos1,pos1+length1-1,pos2,pos2+length2-1);
	      
	  }
	  else{
	      //fprintf(fptr3,"%u %u %u\n",length1,pos1,pos2);
	      if(fptr3!=NULL)fprintf(fptr3,"#%.2f\n",(float)length1);
	      if(fptr3!=NULL)
		  fprintf(fptr3,"[%lu,%lu] [%lu,%lu]\n",pos1,pos1+length1-1,offset-(pos2+length2-1)-1,offset-pos2-1);
	      //fprintf(fptr3,"[%d,%d] [%d,%d]\n",pos1,pos1+length1-1,pos2,pos2+length2-1);
	  }
      }
      if(fptr2!=NULL)fclose(fptr2);
      if(fptr3!=NULL)fclose(fptr3);
  }
  else{ // filetype = 0 a multimat file	      

  }
  if(fptr2!=NULL)fclose(fptr);    
  return(SUCCESS);
}
