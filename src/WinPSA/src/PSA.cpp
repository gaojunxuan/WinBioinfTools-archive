//############################################################################
//Parallel sequence Alignment
//It needs MPICH2 to be installed before. 
//This is Done By: Hisham Adel Hasan
//Research Assistant, Nile University,Cairo,Egypt.
//Bioinformatics group
//E-mail:hisham.mohamed@nileu.edu.eg
//      :hosham2004@yahoo.com
//Date:15 Dec 2008
//############################################################################
#include <iostream>
#include <pthread.h>
//#undef SEEK_END
//#undef SEEK_CUR
//#undef SEEK_SET
#include "mpi.h"
#include "PSA.h"
#include <string>
#include <fstream>
using namespace std;
using std::string;
//int xwz;
int* SData;
int startD;
int endD;
//int* X;
//int* Y;
int Xsize;
int Ysize;
int ranktemp;
  double ct1=0;
    double ct2=0;
double ct=0;
    double calct1=0;
    double calct=0;
double calctto=0;
 double calct2=0;
 
    double proct1=0;
    double proct2=0;
long double Glocation;
MPI_Status status;
//bool check(int x,int y,int& j,int start);
PSA::PSA()
{
    found=false;
    Max=0;
    askp=new int[2];

}

//##########################Open files to broad cast them####################
void PSA::openfiles(int size)
{

     int count=-1;
    char c;
    int inc=0;
    char path1[100];
    char path2[100];
    cout<<"\nEnter Path of file contains Sequence 1:";
    cin>>path1;
    cout<<"\nEnter Path of file contains Sequence 2:";
    cin>>path2;
	
    ifstream infile11(path1);
	
    while (!infile11.eof())
    {
	infile11>>c;
	count++;
    }//end while
    infile11.close();
    PSA::file1s=count;
    cout<<"There is "<<PSA::file1s<<" in file 1"<<endl;
    PSA::file1=new char[count];
    ifstream infile12;
		infile12.open(path1);
	while(!infile12.eof())
	{
	    infile12>>file1[inc];
		//cout<<" f1 = "<<file1[inc]<<endl; 
	    inc++;
	}//end while
    infile12.close();
    count=-1;
    inc=0;
	ifstream infile;
    infile.open(path2);
    while(!infile.eof())
    {
	infile>>c;
	count++;
    }//end while    
    PSA::file2s=count;
    cout<<"There is "<<PSA::file2s<<" in file 2"<<endl;
    PSA::file2=new char[count];
    infile.close();
    ifstream infile21;
		infile21.open(path2);
    while(!infile21.eof())
    {

	infile21>>file2[inc];
	inc++;
    }//end while
    infile21.close();

    PSA::NumOfDiag=PSA::file1s+PSA::file2s-1;
    PSA::DiaPN=PSA::NumOfDiag/size;
    PSA::DiaLN=PSA::NumOfDiag-(PSA::DiaPN*(size-1));
    for(int i=1;i<size-1;i++)
    {

	MPI_Send(&DiaPN,1,MPI_INT,i,10,MPI_COMM_WORLD);   //send num of diagonal that will work on it
    }

        MPI_Send(&DiaLN,1,MPI_INT,size-1,10,MPI_COMM_WORLD);  //send to last node
/*
//###############To display data in the arrays only######################
    for(int s=0;s<file1s;s++)
    {
	cout<<"file1 "<<file1[s]<<endl;

    }
    for(int z=0;z<file2s;z++)
    {
	cout<<"file 2 "<<file2[z]<<endl;

    }
//########################################################################
*/
}//end openfiles

//########################################################################





//###################Bcast to size of file1 and file2#####################
void PSA::Broadcast(int rank)
{
    //xwz=rank;
    ranktemp=rank;
    MPI_Bcast(&(PSA::file1s),1,MPI_INT,0,MPI_COMM_WORLD);
    MPI_Bcast(&(PSA::file2s),1,MPI_INT,0,MPI_COMM_WORLD);
    if(rank!=0)
    {
	file1=new char[file1s];
	file2=new char[file2s];
    }
    MPI_Bcast(file1,file1s,MPI_CHAR,0,MPI_COMM_WORLD);
    MPI_Bcast(file2,file2s,MPI_CHAR,0,MPI_COMM_WORLD);

    /*
    if(rank==2)
    {
	cout<<"File one = "<<file1s<<endl;
	cout<<"File two = "<<file2s<<endl;
	for(int i=0;i<file1s;i++)
	{
	    cout<<"File one contains = "<<file1[i]<<endl;


	}
	for(int i=0;i<file2s;i++)
	{
	    cout<<"File two contains = "<<file2[i]<<endl;


	}




    }//end rank==2
    */

}//end bcast
//########################################################################





//#######################Recv data########################################
void PSA::Recv(int Case,int rank,int size)
{

    switch(Case)
    {
	case 0:
	{
	    //  if(rank!=size-1)
		MPI_Recv(&DiaPN,1,MPI_INT,0,10,MPI_COMM_WORLD,&(PSA::status));
		// else 
		//	MPI_Recv(&DiaLN,1,MPI_INT,0,10,MPI_COMM_WORLD,&(PSA::status));
	    break;
	}//end case 0

    }//end switch

}//end Recv
//########################################################################



//#########################################Count#########################


void PSA::countt(int key,int rank,int size)
{

    int count=0;
  int start=rank*DiaPN;
  STARTD=start;
  ENDD=STARTD+DiaPN-1;
  if(rank==size-1)
	    {
		c=(file1s+file2s-1)/size;
	
	       	start=(c*(size-1));
		STARTD=start;
		ENDD=STARTD+DiaPN-1;
	    }
  for(int i=0;i<PSA::DiaPN;i++)
    {
        if(i==DiaPN-1)
  {
      endD=count;
  }

	if(start>=file1s)
       	{
	    c=start-(file1s-1);
	    r=0;
	    while(r<file1s && c<file2s)
	    {
		/*
		if(key==1)
		{
		    X[count]=r;
		    Y[count]=c;
		}
		*/
		if(r==0 && c==0)
		{

		    Stap=count;
		}
		
	        count++;
                r++;
		c++;
	    }//end while

//if(rank==3)
//{
//cout<<"DiaPn="<<DiaPN<<endl;
//cout<<"1-I am rank "<<rank<<" and i will save "<<count<<endl;
	    //}

	}//end if
	else 
	{
	    c=0;
	r=file1s-start-1;
            while(r<file1s && c<file2s)
	    {
		/*
		if(key==1)
		{
		    X[count]=r;
		    Y[count]=c;
		}
		*/
		if(r==0 && c==0)
		{

		    Stap=count;
		}
		
		count++;
		r++;
		c++;
	    }//end while
	  
	}//end else

  start++;

  if(i==0)
  {

      startD=count;
  }



    }//end for loop
/*
    if(rank==3)
    {

	cout<<"Start D num = "<<startD<<endl;
	cout<<"End D num = "<<endD<<endl;
    }
*/

    if(key==0)
    {
	//  X=new int[count];
//	  Y=new int [count];
    SData=new int [count];
    Xsize=count;
      Ysize=count;


    if(rank==0)
    {
	Glocation2=count;
        MPI_Send(&Glocation2,1,MPI_LONG_DOUBLE,1,502,MPI_COMM_WORLD);
	Glocation=0;
    }
 else if(rank<size-1 && rank!=0)
    {
    MPI_Recv(&Glocation,1,MPI_LONG_DOUBLE,rank-1,502,MPI_COMM_WORLD,&status);
	Glocation2=Glocation+count;
   MPI_Send(&Glocation2,1,MPI_LONG_DOUBLE,rank+1,502,MPI_COMM_WORLD);
    }
    else if (rank==size-1)
    {

    MPI_Recv(&Glocation,1,MPI_LONG_DOUBLE,rank-1,502,MPI_COMM_WORLD,&status);
    }


    }

//####################for display####################
/*
    if(key==1 && rank==6)
    {

	for(int j=0;j<count;j++)
	{
	    cout<<"X= "<<X[j]<<" Y= "<<Y[j]<<endl;


	}

    }
*/
//####################################################
}//end count
//########################################################################


//######################To save data that will be calculated##############
void PSA::Allocate(int rank,int size)
{

    countt(0,rank,size);    //Save location for positions
    //     countt(1,rank,size);    //save positions

}//end allocate
//#######################################################################



/*
extern "C" void *Threadstart(void *);
void *Threadstart(void *_rptr2)
{
    PSA *rptrs=(PSA *)_rptr2;
    cout<<"Run thread"<<endl;
    void *thread=rptrs->exchange(&xwz);




}
*/
void *exchange (void *rptr);
//##############################Start calculation########################
void PSA::Calc(int rank,int size)
{
/*
if(rank==1)
{
  for(int z=0;z<file2s;z++)
    {
	cout<<"file 2 "<<file2[z]<<endl;

    }
  for(int z=0;z<file1s;z++)
    {
	cout<<"file 1 "<<file1[z]<<endl;

    }
}
*/

  	 arr=new int[3];
int ptemp=0;
 int temp1=0;
		 int temp2=0;
		 int temp3=0;
int temp2b=0;
int temp3b=0;
int tempvalue1=0;
int tempvalue2=0;
//    pthread_t thread1;
//    int iret1;
    int incr=0;
    int antidia=0;
    int incCou=0;
    int NXr=0;
    int NXc=0;
    int NXRD=0;
    int NXCD=0;
 int gab;
		 int match;
		 int missmatch;
    int Nextposi=0;
    int Nextposi1=0;
	int X1=0;
	int X2=0;
	int X3=0;
	int X4=0;
int X5=0;
int checkD=0;
	int  tempr=0;
	int  tempc=0;
	int D=0;
int statposi=0;
    bool ache=false;
bool EXChange=false;
int tempvalue=0;
    // bool movecal=false;
    antidia=(DiaPN+1)/2;
    int currposi=0;
   // iret1=pthread_create(&thread1,NULL,exchange,&rank);
    NumOfDiag=file1s+file2s-1;
	bool kPath=false;
    if(rank==0)
    {
    cout<<"Start Calculating......."<<endl;
    cout<<"Please wait....."<<endl;
    }

	if(STARTD>file1s-1)
	{
	    STARTR=0;
	    STARTC=STARTD-file1s+1;
	    relvD=STARTR+STARTC;
	    position=0;

	}
	else if(STARTD<=file1s)
	{
	    if(ENDD>=file1s)
	    {
		STARTR=0;
		STARTC=0;
		relvD=0;
//	found=check(0,0,position,0);
//	Nextposi1=position;
		position=Stap;
		Nextposi1=Stap;
	    }
	    else
	    {
		STARTC=0;
		STARTR=file1s-ENDD-1;
		relvD=STARTR+STARTC;
		position=endD;
		Nextposi1=endD;

	    }

	}





  //if(rank==0)
/*
	    {

		cout<<"STARTD= "<<STARTD<<" ENDD= "<<ENDD<<" Relv D = "<<relvD<<" Iam node "<<rank<<" Xsize= "<<Xsize<<endl;
	    }
*/






ptemp=position;

    for( DCalc=0 ;DCalc<NumOfDiag;DCalc++)  //send diagonal num to  nodes to start calculation
    {
X1=0;
X2=0;
X3=0;
X4=0;
X5=0;


/*
if(relvD==DCalc)
{
cout<<" HEREEEEEEE = "<<rank<<" Diag= "<<DCalc<<endl;
if( ptemp<startD || ptemp>endD )
{
//if(position<Xsize)
{
cout<<"H here rank= "<<rank<<" position= "<<ptemp<<" Xsize= "<<Xsize<<" Diagonal = "<<DCalc<<endl;
if(rank==0)
{
cout<<"1-My rank = "<<rank<<endl;
MPI_Recv(&temp3,1,MPI_INT,rank+1,105,MPI_COMM_WORLD,&status);
cout<<"My rank = "<<rank<<" position = "<<ptemp<<" I got temp3 = "<<temp3<<endl;
}
else if(rank==size-1)
{
cout<<"2-My rank = "<<rank<<endl;
 MPI_Recv(&temp2,1,MPI_INT,rank-1,105,MPI_COMM_WORLD,&status);
cout<<"My rank = "<<rank<<" position = "<<ptemp<<" I got temp2 = "<<temp2<<endl;
}
else
{
cout<<"3-My rank = "<<rank<<" r= "<<r<<" c= "<<c<<endl;
if()
 MPI_Recv(&temp2,1,MPI_INT,rank-1,105,MPI_COMM_WORLD,&status);
cout<<"ouch"<<endl;
MPI_Recv(&temp3,1,MPI_INT,rank+1,105,MPI_COMM_WORLD,&status);
cout<<"My rank = "<<rank<<" position = "<<ptemp<<" I got temp3 = "<<temp3<<" temp2= "<<temp2<<endl;

}


}
}

}
*/


MPI_Barrier(MPI_COMM_WORLD);
	MPI_Bcast(&DCalc,1,MPI_INT,0,MPI_COMM_WORLD);
	int countt=0;

	if(relvD==DCalc)
	{

	    kPath=true;
	    NXCD=STARTC;
	    NXRD=STARTR;
EXChange=true;


	  


	}




	found=false;


	if(kPath==true)
	{
	r=NXRD;
	c=NXCD;
//	found=check(r,c,position,0);
	found=true;
	position=Nextposi1;

	}

	checkD=0;
  tempr=r-c;
	 tempc=c-r;
	 D=0;


//cout<<"r= "<<r<<" c= "<<c<<" Rank= "<<rank<<endl;
	if(tempc==0 && tempr==0)
	{
	checkD=file1s-1;
	D=checkD-1;
	}


else	if(tempr<0 && kPath==true && r<file1s && c<file2s)
	{
	    tempr=0;
	    checkD=tempc+file1s-1;
	    D=checkD-1;
//	position=position-(c-tempc);
	//  cout<<"1-D = "<<D<<endl;
	}
else	if(tempc<0 && kPath==true && r<file1s && c<file2s)
	{
	    tempc=0;
	    checkD=file1s-1-tempr;
	    D=checkD-1;
//	position=position-(r-tempr);
	//  cout<<"2-D = "<<D<<endl;
	}




	if(kPath==true && DiaPN>1 && r<file1s && c<file2s)    //Modified in Windows
	{
	    if(D>=STARTD)
	    {

		if(r+1<file1s && checkD<=ENDD)
		{
		    NXRD=r+1;
		    NXCD=c;

		    Nextposi1=reNP(r,c,2);
/*
		    if(rank==2)
		    {
			cout<<" nexp = "<<Nextposi1<<endl;

		    }


*/
		    Nextposi1=position-Nextposi1+1;

/*
	if(rank==2)
	{
	    cout<<"-NXr = "<<NXRD<<" -NXc= "<<NXCD<<" Nposition= "<<Nextposi1<<endl;
       
	}

*/

   
		}//end r+1<file1s
		else
		{
		    if(checkD==ENDD)
		    {
		    kPath=false;
		    //  found=false;
		    ache=true;
		    //   DCalc=NumOfDiag;
		    }
		    NXRD=r;
		    NXCD=c+1;
		    Nextposi1=reNP(r,c,1);
		    Nextposi1=Nextposi1+position+1;
/*
		    if(rank==1)
			cout<<" Next posi1 = "<<Nextposi1<<endl;
*/
		 
/*
	if(rank==2)
	{
	    //  cout<<"D= "<<D<<endl;
	    cout<<"--NXr = "<<NXRD<<" --NXc= "<<NXCD<<" position= "<<Nextposi1<<endl;
       
	}

*/
		}

	    }//end D>=STARTD
	    else
	    {

		NXRD=r;
		NXCD=c+1;
		Nextposi1=reNP(r,c,1);
/*
		if(rank==2)
		{
		    cout<<"r= "<<r<<" c= "<<c<<" PIO = "<<Nextposi1<<endl;
		}
*/
		if(r==0 && c==0)
		{
		    Nextposi1=Nextposi1+position;

		}
		else
		{
		Nextposi1=Nextposi1+1+position;
		}
/*
	if(rank==2)
	{
	    cout<<"NXr = "<<NXRD<<" NXc= "<<NXCD<<" Position = "<<position<<" Nexposition= "<<Nextposi1<<endl;
       
	}
*/

	
	    }

	}//end KPath==true

	else if(kPath==true && DiaPN<=1)
	{
	    NXRD=r+1;
	    NXCD=c+1;
	    Nextposi1=position+1;
/*
	    if(rank==2)
	    {
	    cout<<"hNXr = "<<NXRD<<" NXc= "<<NXCD<<" Position = "<<position<<" Nexposition= "<<Nextposi1<<endl;
	    }

*/

	}


calct1=MPI_Wtime();
	if(found==true && r<file1s && c<file2s )
	{

	 	incr=0;
	while(incr<antidia)
	{ 
	     {
		 incr++;
		 currposi=position;
		 // gab=0;
		 gab=-2;
		  match=2;
		 missmatch=-1;
	
		
		 if (position<endD && position >=startD)
		 {

		     if(r-1<0)
		     {
			 temp3=0;
			 found2=false;
		     }
		     else
		     {
			 position3=reNP(r,c,1);
			 position3=position+position3;
			 temp3=SData[position3];
			 found2=false;
		     }
		    
		     if(c-1<0)
		     {
			 temp2=0;
			 found2=false;
		     }
		     else
		     {
			 position2=reNP(r,c,2);
			 position2=position-position2;
			 temp2=SData[position2];
			 found2=false;

		     }


		     if(r-1>=0 && c-1>= 0)
		     {
		      position4=position-1;    
			  temp1=SData[position4];
			  found2=false;
		     }

		      else
		          {	
		     	 temp1=0;
		       found2=false;
		     }

		     Value(r,c,match,missmatch,gab,position,temp1,temp2,temp3);



		 }//end if--it was "end else"

		 else  //to start exchange between nodes
		 {

					if(position>=endD)
							{
								X5=1;
							}
					else if(position<startD)
							{
								X5=2;

							}

		     if(r-1<0)
		     {
			 temp3b=0;
			
			 found2=false;
		     }

				 else
		     {
			
			 found2=false;

			 if(position>=endD)
			 {
			
			     posi=reNP(r,c,1);
			     posi=Glocation+position+posi;
			     statposi=posi;
			     found2=false;
if(rank==0 || rank==size-1)
{
temp3b=tempvalue;
}
else
{
temp3b=tempvalue2;
}
			 }//end else found=true



				 else
			 {
			     position3=reNP(r,c,1);
			     position3=position3+position;
		  temp3b=SData[position3];
			 }


		  
                  
		     }//end else r-1<0

		     if(c-1<0)
		     {

            temp2b=0;
			 found2=false;
            statposi=position;

		     }
			     else
		     {

		 if(position<startD)
			 {		   
			     posi=reNP(r,c,2);
			     posi=Glocation+position-posi;
                 statposi=posi;
			     found2=false;
//temp2b=tempvalue1;


if(rank==0 || rank==size-1)
{
temp2b=tempvalue;
}
else
{
temp2b=tempvalue1;
}

			 }

			 else
			 {
			     position2=reNP(r,c,2);
			     position2=position-position2;

			     temp2b=SData[position2];
			 }

		     }//end else if c-1<0
		     if(r-1>=0 && c-1>= 0)
		     {
		
			 	 position4=position-1;
			  temp1=SData[position4];
			  found2=false;
		
		     }
				 else
		     {
		     	 temp1=0;
		       found2=false;

		     }


		     Value(r,c,match,missmatch,gab,position,temp1,temp2b,temp3b);

if(X5==1)
{
					if(rank<size-1)
						{
					MPI_Send(&SData[position],1,MPI_INT,rank+1,100,MPI_COMM_WORLD);

					X1=-1;
							}
X5=0;
}
else if(X5==2)
{
if(rank>0)
{
MPI_Send(&SData[position],1,MPI_INT,rank-1,100,MPI_COMM_WORLD);
//if(rank==1)
//{
//cout<<"X5==2 Iam Node "<<rank<<" I sent  "<<SData[position]<<" r= "<<r<<" c= "<<c<<" to node "<< rank-1<<endl;
//}
X2=-2;
}


X5=0;
}



		 }//end else 		
	     }//end found = true
	    countt++;

	   if(r-1<0 || c+1>=file2s)
	   {
	  
	       incr=DiaPN;
	  ptemp=position;
	   }
	   else
	   {

	       Nextposi=reNP(r,c,1);
	       Nextposi=position+Nextposi+1;
	       Nextposi=reNP(r,c+1,1)+Nextposi;     
			ptemp=position;

	       if(Nextposi>Xsize)
	       {
		   incr=DiaPN;
	
	       }


	        position=Nextposi;
	       r--;
	       c++;

	   }
	 
	   
	}//end while incr
	found=false;
	}//found true


calct2=MPI_Wtime();
calct=calct2-calct1;
calctto=calct+calctto;
	
MPI_Barrier(MPI_COMM_WORLD);



//cout<<"I am here rank= "<<rank<<" position= "<<position<<" Xsize= "<<Xsize<<" StartD= "<<startD<<" endD= "<<endD<<" STARTD= "<<STARTD<<" ENDD= "<<ENDD<<endl;

ct1=MPI_Wtime();
if(rank==0)
{
//cout<<"1-My rank = "<<rank<<endl;
MPI_Send(&X1,1,MPI_INT,1,777,MPI_COMM_WORLD);
//X1=0;


MPI_Recv(&X3,1,MPI_INT,1,888,MPI_COMM_WORLD,&status);

//cout<<" -I got "<<X3<<" from node 1"<<endl;

if(X3==-2)
{

 MPI_Recv(&tempvalue,1,MPI_INT,1,100,MPI_COMM_WORLD,&status);
//cout<<" Iam Node "<<rank<<" I got "<<tempvalue<<" from node 1"<<" Dcalc = "<<DCalc<<endl;
//temp2=tempvalue;
temp3b=tempvalue;

}



X1=0;
X3=0;
//cout<<"My rank = "<<rank<<" position = "<<ptemp<<" I got temp3 = "<<temp3<<endl;
}
else if(rank==size-1)
{
//cout<<"2-My rank = "<<rank<<endl;
MPI_Send(&X2,1,MPI_INT,rank-1,777,MPI_COMM_WORLD);
//X1=0;
 MPI_Recv(&X3,1,MPI_INT,rank-1,888,MPI_COMM_WORLD,&status);
//cout<<"Iam node = "<<rank<<" --I got "<<X3<<" from node "<<rank-1<<endl;
if(X3==-1)
{

  MPI_Recv(&tempvalue,1,MPI_INT,rank-1,100,MPI_COMM_WORLD,&status);
//cout<<" Iam Node "<<rank<<" I got "<<tempvalue<<" from node "<< rank-1<<" Dcalc = "<<DCalc<<endl;
temp2b=tempvalue;
//temp3=tempvalue;
}


//cout<<"My rank = "<<rank<<" position = "<<ptemp<<" I got temp2 = "<<temp2<<endl;
X2=0;
X3=0;
}
else
{
//cout<<"3-My rank = "<<rank<<" r= "<<r<<" c= "<<c<<endl;

MPI_Send(&X2,1,MPI_INT,rank-1,888,MPI_COMM_WORLD);
MPI_Send(&X1,1,MPI_INT,rank+1,888,MPI_COMM_WORLD);
//X1=0;



if(rank-1>0 && rank+1<size-1)
{
 MPI_Recv(&X3,1,MPI_INT,rank-1,888,MPI_COMM_WORLD,&status);
//cout<<" 1-I got "<<X3<<" from node "<<rank-1<<endl;
 MPI_Recv(&X4,1,MPI_INT,rank+1,888,MPI_COMM_WORLD,&status);


}


else if(rank-1==0 && rank+1<size-1)
{
 MPI_Recv(&X3,1,MPI_INT,rank-1,777,MPI_COMM_WORLD,&status);
//cout<<" 1-I got "<<X3<<" from node "<<rank-1<<endl;
 MPI_Recv(&X4,1,MPI_INT,rank+1,888,MPI_COMM_WORLD,&status);
}

else if(rank-1>0 && rank+1==size-1)
{
 MPI_Recv(&X3,1,MPI_INT,rank-1,888,MPI_COMM_WORLD,&status);
//cout<<" 1-I got "<<X3<<" from node "<<rank-1<<endl;
 MPI_Recv(&X4,1,MPI_INT,rank+1,777,MPI_COMM_WORLD,&status);


}


//cout<<" 2-I got "<<X4<<" from node "<<rank+1<<endl;


if(X3==-1)
{

  MPI_Recv(&tempvalue1,1,MPI_INT,rank-1,100,MPI_COMM_WORLD,&status);
//cout<<" Iam Node "<<rank<<" I got "<<tempvalue1<<" from node "<< rank-1<<" Dcalc = "<<DCalc<<endl;
temp2b=tempvalue1;


}
if(X4==-2)
{

  MPI_Recv(&tempvalue2,1,MPI_INT,rank+1,100,MPI_COMM_WORLD,&status);
//cout<<" Iam Node "<<rank<<" I got "<<tempvalue2<<" from node "<< rank+1<<" Dcalc = "<<DCalc<<endl;
temp3b=tempvalue2;
}
//cout<<"ouch"<<endl;
//MPI_Recv(&X3,1,MPI_INT,rank+1,777,MPI_COMM_WORLD,&status);
//cout<<"My rank = "<<rank<<" position = "<<ptemp<<" I got temp3 = "<<temp3<<" temp2= "<<temp2<<endl;
X1=0;
X2=0;
X3=0;
X4=0;
}

ct2=MPI_Wtime();
ct=ct+ct2-ct1;
//cout<<"I am in diagonal "<<DCalc<<" My rank= "<<rank<<endl;
    }//end for NumofDiag


/*
    if(rank!=0)
    {
    posi=-100;
    // MPI_Send(askp,2,MPI_INT,rank-1,100,MPI_COMM_WORLD);
MPI_Send(&posi,1,MPI_INT,rank-1,100,MPI_COMM_WORLD);
    }
    if(rank==size-2)
    {
MPI_Send(&posi,1,MPI_INT,rank+1,100,MPI_COMM_WORLD);
    }
*/
if(rank!=0)
{
 
MPI_Send(&calctto,1,MPI_DOUBLE,0,999,MPI_COMM_WORLD);
MPI_Send(&ct,1,MPI_DOUBLE,0,989,MPI_COMM_WORLD);
//cout<<"calc time=  "<<calct<<endl;;

}
if(rank==0)
{
   double *calctime= new double [size];
 double *commtime= new double [size];
calctime[0]=calctto;
commtime[0]=ct;
for(int i=1; i<size;i++)
{
MPI_Recv(&calctime[i],1,MPI_DOUBLE,i,999,MPI_COMM_WORLD,&status);
MPI_Recv(&commtime[i],1,MPI_DOUBLE,i,989,MPI_COMM_WORLD,&status);
}

 //getmax(calctime,calct);

double maxcal;
maxcal=calctime[0];
//cout<<"calc time=  "<<calctime[0]<<endl;;
    for(int f=1;f<3;f++)
    {
//cout<<"calc time=  "<<calctime[f]<<endl;;
	if(calctime[f]>maxcal)
	{
	    maxcal=calctime[f];

	}
    }
calct=maxcal;
//cout<<" max = "<<max<<endl;

maxcal=commtime[0];
//getmax(commtime,ct);
   for(int f=1;f<3;f++)
    {
//cout<<"calc time=  "<<commtime[f];
	if(commtime[f]>maxcal)
	{
	    maxcal=commtime[f];

	}
    }

ct=maxcal;
//cout<<"commtime=  "<<ct;






}


}//end calc
//########################################################################

//######################check for the data###############################
/*
bool check(int x,int y,int& j,int start)
{
     j=start;
     // if(rank==0)
     //cout<<"Start = "<<start<<endl;
    while(j<Xsize)
    {
	if(x!=X[j])
	{
	    j++;
	}
	else
	{
	    if(y!=Y[j])
		j++;
	    else
		return true;

	}

    }//end while

    return false;

}
*/
//########################################################################



//########################Get Max########################################

void PSA::getmax(int arr[],int& max)
{
    max=arr[0];
    for(int f=1;f<3;f++)
    {
	if(arr[f]>max)
	{
	    max=arr[f];
	}
    }
//cout<<" max = "<<max<<endl;




}
//######################################################################
//##################TO save Data##########################################
void PSA::MaxSave(int array[],int pos)
{

	 getmax(array,Max);
		 SData[pos]=Max;



}
//######################################################################
//#############################Value####################################
void PSA::Value(int row,int colom,int match,int missmatch,int gab,int posi,int temp1,int temp2,int temp3)
{


	     if(file1[row]==file2[colom])
		     {
			 arr[0]=temp1+match;
		
		     }
		     else if (file1[row]!=file2[colom])
		     {
			 arr[0]=temp1+missmatch;
			 // cout<<"2-arr[0] = "<<arr[0]<<endl;
		     }
		     arr[1]=temp2+gab;
		     arr[2]=temp3+gab;
		     MaxSave(arr,posi);






}
//################################################################################

//################################TO Exchange data between #####################
void * exchange(void* rptr)
{
    int sear;
    // static int sear1=0;
    //  static int sear2=endD;
    int position5=0;
bool found3=false;
    int *ToRank;
    ToRank=(int *)rptr;
    int rank1= *ToRank;
    int chk=-1;
 int *wanp=new int[2];
 //PSA sa;
 wanp[0]=0;
 wanp[1]=0;
 //cout<<"Thread run on Node "<<rank1<<endl;
 while(chk!=-100)
 {
     //cout<<"In While"<<endl;
    

     //   MPI_Recv(wanp,2,MPI_INT,MPI_ANY_SOURCE,100,MPI_COMM_WORLD,&status);
     MPI_Recv(&sear,2,MPI_INT,MPI_ANY_SOURCE,100,MPI_COMM_WORLD,&status);
     // chk=wanp[0];
     chk=sear;
     //cout<<"Node "<<status.MPI_SOURCE<<" ask for r= "<<wanp[0]<<" c= "<<wanp[1]<<" I am Node "<<rank1<<endl;
     if(chk!=-100)
     {

	 if(rank1<status.MPI_SOURCE)
	 {
	     //   if(rank1==0)
	     //  cout<<"sear2 = "<<sear2<<endl;



	     // found3=check(wanp[0],wanp[1],position5,sear2);
	     // sear2=position5;

	     // position5=reNP(wanp[0],wanp[1],2);


	     /*
	     if(rank1==1)
	     {
		 cout<<"-r= "<<wanp[0]<<" -c= "<<wanp[1]<<" position5= "<<position5<<" Data= "<<SData[position5]<<endl;
	     }
	     */
	     

	     sear=sear-Glocation;



	 }
	 else if(rank1>status.MPI_SOURCE)
	 {
	     // if(rank1==0)
		  //  cout<<"sear1 = "<<sear1<<endl;

	     //   sear=sear-Glocation;

	     // found3=check(wanp[0],wanp[1],position5,sear1);
	     // sear1=position5;
	     //   position5=reNP(wanp[0],wanp[1],1);

	     sear=sear-Glocation;

/*

   if(rank1==1)
	     {
		 cout<<"--r= "<<wanp[0]<<" --c= "<<wanp[1]<<" Position5= "<<position5<<" Data= "<<SData[position5]<<endl;
	     }
*/




	 }

	 //cout<<"Found = "<<found3<<" it is in position "<<position5<<"r= "<<wanp[0]<<"c= "<<wanp[1]<<" I am Node "<<rank1<<endl;

	 // MPI_Send(&SData[position5],1,MPI_INT,status.MPI_SOURCE,105,MPI_COMM_WORLD);
MPI_Send(&SData[sear],1,MPI_INT,status.MPI_SOURCE,105,MPI_COMM_WORLD);


     }//end 


 }//end while


 //cout<<"out of while Rank ="<<rank1<<endl;

//pthread_exit(NULL);
return 0;
}
//##############################################################################
//##################Gather Data####################################
void PSA::gatherdata(int rank,int size)
{
if(rank==0)
{
cout<<"Start collecting data....."<<endl;


}
delete [file1s] file1;
delete [file2s] file2;
//file1=NULL;
//file2=NULL;


if(rank==0)
{
 
	
    savearr=file1s*file2s;
cout<<"Start saving the array of results....."<<endl;

//ct=ct+ct2-ct1;
 ofstream outfile("./time.txt");
outfile<<"Calc time= "<<calct<<endl;
outfile<<"Comm Time= "<<ct<<endl;
resl=new int [savearr];

}

//=new int [savearr];


 
  




//MPI_Barrier(MPI_COMM_WORLD);
 //file2=NULL;
 //file1=NULL;

    int rwx=1;
  
    int tempsize;
    int *cluster;
    cluster=new int [size];
    int *dis;
    dis=new int [size];
    if(rank!=0)
    {

	MPI_Send(&Xsize,1,MPI_INT,0,200,MPI_COMM_WORLD);
    }
    else
    {
cluster[0]=Xsize;	
	while(rwx<size)
	{
	    //cout<<"3-size = "<<size<<endl;
MPI_Recv(&tempsize,1,MPI_INT,MPI_ANY_SOURCE,200,MPI_COMM_WORLD,&status);
cluster[status.MPI_SOURCE]=tempsize;
rwx++;

	}
	dis[0]=0;
	for(int h=1;h<size;h++)
	{
	    dis[h]=cluster[h-1]+dis[h-1];


	}

    }







ct1=MPI_Wtime();
MPI_Barrier(MPI_COMM_WORLD);
    MPI_Gatherv(SData,Xsize,MPI_INT,resl,cluster,dis,MPI_INT,0,MPI_COMM_WORLD);
ct2=MPI_Wtime();


if(rank==0)
{
ct=ct+ct2-ct1;
 ofstream outfile("./time.txt");
outfile<<"Calc time= "<<calct<<endl;
outfile<<"Comm Time= "<<ct<<endl;


}


    if(rank==0)
    {
	int **Output;
	Output=new int *[file1s];
	int k=0;
	for(int i=0;i<file1s;i++)
	{
	    Output[i]=new int [file2s];

	}



	int start=0;
	int count=0;
 for(int i=0;i<NumOfDiag;i++)
    {
 
	if(i>=file1s)
       	{
	    c=start-(file1s-1);
	    r=0;
	    while(r<file1s && c<file2s)
	    {
	
		Output[r][c]=resl[count];
	        count++;
                r++;
		c++;
	    }//end while



	}//end if
	else 
	{
	    c=0;
	r=file1s-start-1;
            while(r<file1s && c<file2s)
	    {
	
		Output[r][c]=resl[count];
		count++;
		r++;
		c++;
	    }//end while
	  
	}//end else
  start++;

    }//end for loop





/*

 for(int i=0;i<file1s;i++)
 {
     for(int j=0;j<file2s;j++)
     {
	 cout<< Output[i][j]<<"\t";

     }
     cout<<endl;
 }
*/
 ofstream outfile("./OutPutPSA.txt");

 for(int i=0;i<file1s;i++)
 {
     for(int j=0;j<file2s;j++)
     {
	 outfile<<Output[i][j]<<"\t";
//cout<<"j= "<<j<<endl;
     }
     outfile<<endl;
//cout<<"---------------i= "<<i<<" --------------------------"<<endl;
 }

    }//end rank =0


 
}//end gatherdata
//##########################################################################
/*
//########################Getposition#####################################3
int PSA::getposition(int r,int c,int position)
{

    int NXr;
    int NXc;
    int tempr;
    int tempc;
    NXr=r-1;
    NXc=c+1;
    int D;
    //   cout<<"NXr= "<<NXr<<endl;
    // cout<<"NXc= "<<NXc<<endl;
       if(NXr<0 || NXc<0)
       {
	   return 0;
       }

    tempr=r-c;
    tempc=c-r;
    //  cout<<"tempr = "<<tempr<<endl;
    //  cout<<"tempc = "<<tempc<<endl;
    if(tempr<0)
    {
	tempr=0;
	D=tempc+file1s-1;
	position=position-(c-tempc);
	//  cout<<"1-D = "<<D<<endl;
    }
    if(tempc<0)
    {
	tempc=0;
	D=file1s-1-tempr;
	position=position-(r-tempr);
	//  cout<<"2-D = "<<D<<endl;
    }
 
    if(	tempc==0 && tempr==0)
    {
	D=file1s-1;
	position=position-r;
	//  cout<<"3-D = "<<D<<endl;
    }



    while(r!=NXr ||  c!=NXc)
	    {
		if(D>=file1s)
		{
		    //   cout<<"hereeee"<<endl;
		    c=D-(file1s-1);
		    r=0;
		    // cout<<"**r= "<<r<<" **c= "<<c<<endl;
		    while(r<file1s && c<file2s)
		    {
			//	cout<<"1-r= "<<r<<" 1-c= "<<c<<" 1-Position = "<<position<<endl;

			if(r==NXr &&  c==NXc)
			{

			    return position;
			}
			position++;
			r++;
			c++;
//	cout<<"1-Position = "<<position<<endl;
		    }//end while

//if(rank==3)
//{
//cout<<"DiaPn="<<DiaPN<<endl;
//cout<<"1-I am rank "<<rank<<" and i will save "<<count<<endl;
	    //}

	}//end if
		else 
		{
		    //  int fxz;
		    //  cin>>fxz;
		    c=0;

		    r=file1s-D-1;

		    while(r<file1s && c<file2s)
		    {
			//     	cout<<"2-r= "<<r<<" 2-c= "<<c<<" 2-Position = "<<position<<endl;
		
			if(r==NXr &&  c==NXc)
			{
			   return position;
			}
			position++;
			r++;
			c++;
		
		    }//end while
		    
		}//end else
		
		D++;

	    }




		if(r==file1s-1)
		{
		    r=0;
		    D++;
		    c=D-(file1s-1);
		}
		if(c==file2s-1)
		{
		    c=0;
		    D++;
		    r=file1s-D-1;
		}
	
	        position++;
                r++;
		c++;
		cout<<"r = "<<r<<" c= "<<c<<endl;
		cout<<"D ++ = "<<D<<endl;




 


}//end getposition
*/

//*****************************************************************
int PSA::reNP(int r,int c,int key)
{
  int tempr;
    int tempc;
    int D;

  tempr=r-c;
    tempc=c-r;
    //  cout<<"tempr = "<<tempr<<endl;
    //  cout<<"tempc = "<<tempc<<endl;
    if(tempr<0)
    {
	tempr=0;
	D=tempc+file1s-1;
//	position=position-(c-tempc);
	//  cout<<"1-D = "<<D<<endl;
    }
    if(tempc<0)
    {
	tempc=0;
	D=file1s-1-tempr;
//	position=position-(r-tempr);
	//  cout<<"2-D = "<<D<<endl;
    }
 
    if(	tempc==0 && tempr==0)
    {
	D=file1s-1;
//	position=position-r;
//	  cout<<"3-D = "<<D<<endl;
    }
if(file1s<=file2s)
{

    if(D<file1s)
    {
	if(key==1)
	{

	 
	    if(c==0)
	    {
		return (file1s-r);

	    }//end c==0
	    else if(D==file1s-1)
	    {
//	  cout<<"4-D = "<<D<<endl;
		return(file1s-1);

	    }
	    else
	    {
		return (file1s-r+c);
	    }//end else

	}//end key ==1

	else if(key==2)
	{

	    return (file1s-tempr);


	}

    }//end D<file1s
    else if(D>=file1s && D<file2s)
    {
	if(key==1)
	{
	   
	    
	    return (file1s-1);

	}//end key==1
	else if(key==2)
	{
	    return (file1s);
	}

    }//end D>=file1s && D<file2s
    else
    {
	if(key==1)
	{
	    return(file2s-c+(r-1));

	}//end key==1
	else if(key==2)
	{
	    return (file2s-tempc+1);

	}//end else if key==2




    }//end else

}//end if file1<file2
else
{


	    if(key==1)
	    {
		if(c==0)
		{
		    if(D<file2s-1)
		    {

			return (file1s-r);
		    }

		    else if (D>=file2s-1 && D<=file1s-1)
		    {
			return(file2s);

		    }
		   
		  

		}//end c==0

		else
		{
		    if(D>=file1s-1)
		    {

			return (file2s-c+r-1);
		    }

		    else  if (D< file2s-1)
		    {
			return (file1s-r+c);
		    }
		    else
		    {
			return(file2s);

		    }
	



		}

	    }//end if key ==1

	    else if(key==2)
	    {


	
		if(D<=file2s-1)
		{
		    return(file1s-tempr);

		}//end if c==0
		else if(D>file2s-1 && D<=file1s-1)
		{


            	    return(file2s+1);


		}
		else if (D>file1s-1)
		{
		    return(file2s-tempc+1);

		}



	    }

}

}

//****************************************************************
