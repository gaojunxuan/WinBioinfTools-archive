//############################################################################
//Formatdb.c
//It needs MPICH2 to be installed before. 
//This is Done By: Hisham Adel Hasan
//Research Assistant, Nile University,Cairo,Egypt.
//Bioinformatics group
//E-mail:hisham.mohamed@nileu.edu.eg
//      :hosham2004@yahoo.com
//Date:15 Dec 2008
//############################################################################
#define _LARGEFILE_SOURCE 
#define _LARGEFILE64_SOURCE 
#define _CRT_SECURE_NO_DEPRECATE
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>
#include<string.h>
void Readfile();
void SlaveWork();
//Blast Program
//********************************************************************************
unsigned long long nvalues;          //num of elements in the file/number of nodes
   unsigned long long i;
   unsigned long long j;
	unsigned long long l;
   int node;            //to determine the node
  unsigned long long inc;            //to determine start of each node
    FILE *infile;          //to open the data file
    FILE *outfile;        //save data on master node
    unsigned long long f;
    int Ready=0;
    int error=0;
    char *revbuff;
  unsigned long long rbuffsize=0;
//  rbuffsize=0;
    int chec=0;
    FILE *SlaveOutfile;
    int FS=0;
   int rank,size;        //rank of each node  size-Number of nodes
    char *b;
    char temp; 
  unsigned long long start;
  unsigned long long end;
 unsigned long long nvf=0;            //used to count elements in the file
    char *ExtraArray;
    int dbt;
 unsigned long long exarrsz;// size of extra array;
	char ch;
unsigned long long countt;
unsigned long long intdiv=0;
	int coIntDiv=0;
	char *SubExarray;
	//int intsize=1000;
	int intsize=2147483640;
	int IntIncr=0;
	unsigned long long Intexdata=0;
    MPI_Status status;
    MPI_Request request;
//*********************************************************************************
main (int argc,char* argv[])
{
 
    MPI_Init(&argc,&argv);
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    MPI_Comm_size(MPI_COMM_WORLD,&size);

	system("mkdir C:\\db");

 printf("\n\n");
fflush (stdout);
    if(rank==0)
    {
       Readfile();


    }//end if rank=0

 //BroadCast the querry size to all nodes
//**********************************************************************************
//The following code done on all nodes Except Master nodes
//**********************************************************************************
    if(rank!=0)
    {
	SlaveWork();

    }// if rank!= 0

//***********************************************************************************
//Start Format the data base and start search..........
//***********************************************************************************

//Format DataBase
MPI_Barrier(MPI_COMM_WORLD);               //to make sure that all nodes reach to this Point
    
 if(rank==0)
      {

if(dbt==2)
{
printf("\n Database is DNA sequence.....\n");
fflush (stdout);
error=system("D:\\Blastout\\bin\\formatdb -i C:\\BLASTDB\\data.nt -p F -o T");
	if (error!=0)
	{
	    printf("\nCannot Format Database on the Node %d ....\n",rank);
		fflush (stdout);
	}//end if
	else 
	    printf("\nData formated on node %d ...\n",rank);
fflush (stdout);
            free(b);
}
else if (dbt==1)
{
printf("\n Database is protein sequence.....\n");
fflush (stdout);
error=system("D:\\Blastout\\bin\\formatdb -i C:\\BLASTDB\\data.nt -p T -o T");
	if (error!=0)
	{
	    printf("\nCannot Format Database on the Node %d ....\n",rank);
		fflush (stdout);
	}//end if
	else 
	    printf("\nData formated on node %d ...\n",rank);
	fflush (stdout);

            free(b);

}

      }

     if(rank!=0)
      {
if(dbt==2)
{
printf("\n Database is DNA sequence.....\n");
fflush (stdout);
error=system("D:\\Blastout\\bin\\formatdb -i C:\\BLASTDB\\data.nt -p F -o T");
	if (error!=0)
	{
          printf("\nCannot Format Database on the Node %d ....\n",rank);
		  fflush (stdout);
	}
	else 
	    printf("\nData formated on node %d ...\n",rank);
	fflush (stdout);

}
else if (dbt==1)
{
printf("\n Database is protein sequence.....\n");
fflush (stdout);
error=system("D:\\Blastout\\bin\\formatdb -i C:\\BLASTDB\\data.nt -p T -o T");
	if (error!=0)
	{
          printf("\nCannot Format Database on the Node %d ....\n",rank);
		  fflush (stdout);
	}
	else 
	    printf("\nData formated on node %d ...\n",rank);
fflush (stdout);

}

      }
   

     MPI_Finalize();        //end all MPI Functions
}//End main
//*************************Function Readfile()**********************************************
void Readfile()
{
FILE *infile1;
char filePath[200];
int rep;
temp='x';
printf("\n*****************************************\n");
printf("\n***Welcome to Our Program***\n");
//printf("\nPlease Enter '1' if file size less than 2GB\nand '2' if file size more than 2GB:");
//fflush (stdout);
//scanf("%d",&FS);
FS=1;
printf("\nPlease Enter '1' if it is Protien sequence \nand '2' if it is DNA sequence:");
fflush (stdout);
scanf("%d",&dbt);
MPI_Bcast(&dbt,1,MPI_INT,0,MPI_COMM_WORLD); 
for(rep=0;rep<1;rep++)
{
printf("\nPlease Enter the Path of the File:");
fflush (stdout);
	scanf("%s",filePath);
	//printf("\n file path= %s",filePath);
	//fflush (stdout);
if(FS==2)
{
//infile=fopen64(filePath,"r");      
}
else if (FS==1)
{
infile=fopen(filePath,"r");

}
	if(infile==NULL)
	{
		printf("\nError: %s ", strerror(errno));
            //printf("\nError:Can not Open The File.....\n");
			fflush (stdout);
           rep=0;

	}//end if infile=null

}  //end for Repeat

//***************Open the file and determine number of elements in the file**************************************
while(!feof(infile))                   //determine number of char in the file
	{
	    fscanf(infile,"%c",&temp);    //read element by element and but them in temp
	  
	    nvf=nvf+1;                     //count element in the file
	
	}//end while
	printf("\nThere is ' %d ' Char in the file......\n",nvf);       //show how many element in the file
	fflush (stdout);
	fclose(infile);                                    //close the file
//**************************************************************************************************************
if(FS==2)
{
//infile=fopen64(filePath,"r");      
}
else if (FS==1)
{
infile=fopen(filePath,"r");
}
//*****************************Start reading the file******************************
if (size>1)
{
nvalues=nvf/size;   //divide number of char in the file on number of nodes
printf("nvalues=%d",nvalues);
printf("Size=%d",size);
fflush (stdout);
for(i=0;i<nvalues;++i)
{
temp=fgetc(infile);

}//end for
countt=nvalues;
while (temp!='>')
{
temp=fgetc(infile);
countt++;
}
infile1=fopen(filePath,"r");
outfile=fopen("C:\\BLASTDB\\data.nt","w");
for(i=0;i<countt-1;++i)
{
  fscanf(infile1,"%c",&temp);  //start reading elements from file and put them in b array
fprintf(outfile,"%c",temp);

}//end for
    fclose(outfile);
printf("\nData saved on master node....\n");
fflush (stdout);
}

else
{
infile1=fopen(filePath,"r");
outfile=fopen("C:\\BLASTDB\\data.nt","w");
for(i=0;i<nvf;++i)
{
  fscanf(infile1,"%c",&temp);  //start reading elements from file and put them in b array
fprintf(outfile,"%c",temp);

}//end for
    fclose(outfile);
fclose(infile1);
printf("\nData saved on master node....\n");
fflush (stdout);

}

//***************************************************************************************
//*************************************************************************
//Type:for loop
//            Move on all Nodes and send the data to each node
//
//
//*************************************************************************
f=0;
inc=0;
temp='x';
countt=0;




//*******************************NEW PART***********************************************
intdiv=nvalues/intsize;
printf("\Intsize div= %d ....\n",intsize);
fflush (stdout);
while(intdiv>0)
{
coIntDiv++;

intdiv=intdiv/intsize;
printf("\Inint div= %d ....\n",intdiv);
fflush (stdout);

}

printf("\coIntDiv= %d ....\n",coIntDiv);
fflush (stdout);


 MPI_Bcast(&coIntDiv,1,MPI_INT,0,MPI_COMM_WORLD); 

//*****************************************************************************************

for(node=1;node<size;node++)       //to pass on all nodes
{
 error=MPI_Recv(&Ready,1,MPI_INT,node,13,MPI_COMM_WORLD,&status);             //Check if the required node is ready or not
 if(error==MPI_SUCCESS)     //check the error 
 {
    if(Ready==123)                       //if node is ready it will send 123 if not ready will =0
    {
	printf("\nNode %d Ready to get data.......",node);
	fflush (stdout);


    if(node!=size-1)                // to make last node take all data till the end of the array
    {

    
		for(i=0;i<nvalues;i++)
				{
					temp=fgetc(infile);
						countt++;
				}
		while (temp!='>')
				{
					temp=fgetc(infile);
						countt++;
				}
  if(coIntDiv<=0)
	  {
					exarrsz=countt;           //determine the size of the array of data

				MPI_Send(&exarrsz,1,MPI_UNSIGNED_LONG_LONG,node,50,MPI_COMM_WORLD);  //start sending size of  the array to the node

				ExtraArray=(char *)malloc(exarrsz * sizeof(char));     //save data you want to send in extraarray 

		for(i=0;i<countt;i++)  //start to save wanted data in extraarray 
				 {
				fscanf(infile1,"%c",&ExtraArray[f]);
					f++;
				 }//end for 

				MPI_Send(ExtraArray,exarrsz,MPI_CHAR,node,86,MPI_COMM_WORLD); //now start to send the extraarray to a node

		f=0;
		temp='x';
		free( ExtraArray);
		countt=0;
	  
	  }//end if coIntDiv<=0
//***************************new PArt*********************************************
	  else if(coIntDiv>0)
	  {

exarrsz=countt;           //determine the size of the array of data

ExtraArray=(char *)malloc(exarrsz * sizeof(char));     //save data you want to send in extraarray 
SubExarray=(char *)malloc(intsize * sizeof(char));
 for(i=0;i<countt;i++)  //start to save wanted data in extraarray 
    {
fscanf(infile1,"%c",&ExtraArray[f]);
	f++;
    }//end for 

IntIncr=0;
Intexdata=0;
for(j=0;j<coIntDiv;j++)
{
	
for(l=IntIncr;l<intsize;l++)
{

SubExarray[l]=ExtraArray[Intexdata];
IntIncr++;
Intexdata++;
}

 MPI_Send(&intsize,1,MPI_INT,node,50,MPI_COMM_WORLD);  //start sending size of  the array to the node

 MPI_Send(SubExarray,intsize,MPI_CHAR,node,86,MPI_COMM_WORLD); //now start to send the extraarray to a node

IntIncr=0;
}
if(Intexdata<exarrsz)
{
	//printf("Index data = %d",Intexdata);
	//printf("\nExarray = %d\n",exarrsz);
	//printf("Exarr-Intexdata= %d",exarrsz-Intexdata);
//fflush (stdout);
	IntIncr=exarrsz-Intexdata+1;
SubExarray=(char *)malloc(IntIncr * sizeof(char));
for(l=0;l<IntIncr;l++)
{

SubExarray[l]=ExtraArray[Intexdata-1];
Intexdata++;

}
//printf("IntIncr befor sending = %d",IntIncr);
//fflush (stdout);
MPI_Send(&IntIncr,1,MPI_INT,node,50,MPI_COMM_WORLD);  //start sending size of  the array to the node

 MPI_Send(SubExarray,IntIncr,MPI_CHAR,node,86,MPI_COMM_WORLD); //now start to send the extraarray to a node


}
else
{
IntIncr=-5;
MPI_Send(&IntIncr,1,MPI_INT,node,50,MPI_COMM_WORLD);  //start sending size of  the array to the node

}

    f=0;
temp='x';

countt=0;
free( ExtraArray);

	  }//end else if coIntDiv>0
//************************************************************************************************************


    }// node !=size-1;
    else                     //in case of it is last node it will take last start and end =nvf
    {
		if(coIntDiv<=0)
		{
countt=0;
f=0;
while(!feof(infile))                   //determine number of char in the file
	{
	    //fscanf(infile,"%c",&temp); 
//countt++;
temp=fgetc(infile);
countt++;
}
	exarrsz=countt;

	MPI_Send(&exarrsz,1,MPI_UNSIGNED_LONG_LONG,node,50,MPI_COMM_WORLD);             //send the size of the array

	ExtraArray=(char *)malloc(exarrsz * sizeof(char));  
   
while(temp=='\n')
{                 //save data you want to send in ExtraArray
temp=fgetc(infile1);
}
while(!feof(infile1))
    {
        fscanf(infile1,"%c",&ExtraArray[f]);
	f++;
    }//end for 
  
  MPI_Send(ExtraArray,exarrsz,MPI_CHAR,node,86,MPI_COMM_WORLD);              //send the Array

free( ExtraArray);

		}//end if coIntDiv<=0


//*********************************new part********************************************
		
		else if (coIntDiv>0)
		{
//printf("\n nvf>intsize\n");
//fflush (stdout);
			countt=0;
f=0;
while(!feof(infile))                   //determine number of char in the file
	{
	    //fscanf(infile,"%c",&temp); 
//countt++;
//temp=fgetc(infile);

	//printf("\nIn while 2\n");
temp=fgetc(infile);
//printf("temp=%c",temp);
countt++;
}
	exarrsz=countt;

ExtraArray=(char *)malloc(exarrsz * sizeof(char));     //save data you want to send in extraarray 
SubExarray=(char *)malloc(intsize * sizeof(char));


while(!feof(infile1))
    {
        fscanf(infile1,"%c",&ExtraArray[f]);
	f++;
    }//end for 


IntIncr=0;
Intexdata=0;
for(j=0;j<coIntDiv;j++)
{
	
for(l=IntIncr;l<intsize;l++)
{

SubExarray[l]=ExtraArray[Intexdata];
IntIncr++;
Intexdata++;
}
//printf("IntIncr befor sending = %d",IntIncr);
//fflush (stdout);
 MPI_Send(&intsize,1,MPI_INT,node,50,MPI_COMM_WORLD);  //start sending size of  the array to the node

 MPI_Send(SubExarray,intsize,MPI_CHAR,node,86,MPI_COMM_WORLD); //now start to send the extraarray to a node

IntIncr=0;
}

if(Intexdata<exarrsz)
{
	
//fflush (stdout);
free(SubExarray);
IntIncr=exarrsz-Intexdata+1;
SubExarray=(char *)malloc((exarrsz-Intexdata) * sizeof(char));
for(l=0;l<IntIncr;l++)
{

SubExarray[l]=ExtraArray[Intexdata-1];
Intexdata++;

}

MPI_Send(&IntIncr,1,MPI_INT,node,50,MPI_COMM_WORLD);  //start sending size of  the array to the node

 MPI_Send(SubExarray,IntIncr,MPI_CHAR,node,86,MPI_COMM_WORLD); //now start to send the extraarray to a node


}//end (Intexdata<exarrsz)
else
{
IntIncr=-5;
MPI_Send(&IntIncr,1,MPI_INT,node,50,MPI_COMM_WORLD);  //start sending size of  the array to the node

}


	//MPI_Send(&exarrsz,1,MPI_LONG_LONG,node,50,MPI_COMM_WORLD);             //send the size of the array

	//ExtraArray=(char *)malloc(exarrsz * sizeof(char));  
   
//while(temp=='\n')
//{                 //save data you want to send in ExtraArray
//temp=fgetc(infile1);
//}
	//infile1=infile;

  
  //MPI_Send(ExtraArray,exarrsz,MPI_CHAR,node,86,MPI_COMM_WORLD);              //send the Array

free( ExtraArray);

		}//end else if nvf>intsize
//************************************************************************************************
    }//end else node size
    printf("\nNode %d got the data.....",node); 
	fflush (stdout);
    Ready=0;
	}//end if ready==123
 }
    else
    {

	printf("\nNode %d Not Ready to get Data......Cluster will neglect Node %d in this querry\n",node,node);
       fflush (stdout);

    }
}//end for node

fclose(infile);
fclose(infile1);
}//end read file function

//******************************end Read File************************************************************

//**********************************Function Slave Work**************************************************
void SlaveWork()
{
	SlaveOutfile=fopen("C:\\BLASTDB\\data.nt","w");
        Ready=123;
MPI_Bcast(&dbt,1,MPI_INT,0,MPI_COMM_WORLD); 
MPI_Bcast(&coIntDiv,1,MPI_INT,0,MPI_COMM_WORLD); 
 MPI_Send(&Ready,1,MPI_INT,0,13,MPI_COMM_WORLD);
//printf("\nCoIntDiv=%d\n",coIntDiv);
//printf("\ndbt=%d\n",dbt);
//fflush (stdout);
if(coIntDiv<=0)
{
	
//printf("HHHHHHH in node %d count div less 0 and = %d",rank,coIntDiv);
//fflush (stdout);
       // MPI_Send(&Ready,1,MPI_INT,0,13,MPI_COMM_WORLD);

error=MPI_Recv(&rbuffsize,1,MPI_UNSIGNED_LONG_LONG,0,50,MPI_COMM_WORLD,&status);     //recieve size of the file

  revbuff=(char *)malloc((rbuffsize) * sizeof(char));    //save in memory space to save comming data

if(error==MPI_SUCCESS)
{

MPI_Recv(revbuff,rbuffsize,MPI_CHAR,0,86,MPI_COMM_WORLD,&status);     //save the data in revbuff

//SlaveOutfile=fopen("C:\\BLASTDB\\data.nt","w");

  for(i=0;i<rbuffsize-1;i++)  //To save data from 0 to nvalues+inc
    {
	fprintf(SlaveOutfile,"%c",revbuff[i]);

    }//end for loop
fclose(SlaveOutfile);
free(revbuff);
}//end if 
}//end if coIntDiv<0

else 
{
	
for(j=0;j<coIntDiv;j++)
{

error=MPI_Recv(&rbuffsize,1,MPI_LONG_LONG,0,50,MPI_COMM_WORLD,&status);     //recieve size of the file

revbuff=(char *)malloc((rbuffsize) * sizeof(char));    //save in memory space to save comming data
if(error==MPI_SUCCESS)
{
MPI_Recv(revbuff,rbuffsize,MPI_CHAR,0,86,MPI_COMM_WORLD,&status);     //save the data in revbuff
	//printf("\n\nHHHH in else state for loop-Start write\n\n");
fflush (stdout);
for(i=0;i<rbuffsize-1;i++)  //To save data from 0 to nvalues+inc
    {
	fprintf(SlaveOutfile,"%c",revbuff[i]);
//printf("\n%c\n",revbuff[i]);
//fflush (stdout);
    }//end for loop
//printf("Node %d in count div",rank);
}

free(revbuff);
}
MPI_Recv(&rbuffsize,1,MPI_LONG_LONG,0,50,MPI_COMM_WORLD,&status);
//printf("buff size after receving with rank %d= %d",rank,rbuffsize);
fflush (stdout);
if(rbuffsize!=-5)
{
	//printf("Node %d in rbuff!=-5",rank);
revbuff=(char *)malloc(rbuffsize * sizeof(char));    //save in memory space to save comming data

MPI_Recv(revbuff,rbuffsize,MPI_CHAR,0,86,MPI_COMM_WORLD,&status);     //save the data in revbuff
//printf("\nbuff size before start writing rank %d= %d",rank,rbuffsize);
fflush (stdout);
for(i=0;i<rbuffsize-1;i++)  //To save data from 0 to nvalues+inc
    {
	fprintf(SlaveOutfile,"%c",revbuff[i]);
//printf("\n%c\n",revbuff[i]);
    }//end for loop
//printf("\nuff size after writing rank%d=%d\n",rank,rbuffsize);
//fflush (stdout);
}

fclose(SlaveOutfile);
free(revbuff);


}



}//end function SlaveWork
//*********************************************************************************************************//
