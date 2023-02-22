#include <iostream>
#include <fstream>
using std::cout;

#include <string>
#include "assert.h"
#include "math.h"

#include<iterator>
using std::ostream_iterator;

#include "Transformer.h"
#include "myglobaldef.h"

static unsigned int stat_dim;

#ifndef max
#define max(a, b)  ((a)>=(b) ? (a) : (b))
#endif
#ifndef min
#define min(a, b)  ((a)<=(b) ? (a) : (b))
#endif


PlainTransformer::PlainTransformer(char* datFile, unsigned int dim):Transformer(datFile, dim){
    
}

PlainTransformer::~PlainTransformer(){
}


int PlainTransformer::transformIt(){
    parseDatFilekGenomes(inDatFile);
    // add shiftvalue
    for(unsigned long i=0;i<Chain.size();i++){
	for(long j=0;j<dimensions;j++){
	    Chain[i].b.coordinates[j]=Chain[i].b.coordinates[j]+shiftvalue;
	    Chain[i].e.coordinates[j]=Chain[i].e.coordinates[j]+shiftvalue;
	}
    }
    //cout << "# Overlap flag= " << overlap_flag<< endl;
    if(filterrep_flag==1){
	filter_repeats_all_dimensions();
	//cerr << "HALOOOOOOOOOOOOOO\n";
	report_repeats();
    }
    

    chooseFragments();
    report_permutations();
//    report_chain();
    plot_compact_chain();        
/*
     for(int k=0;k<dimensions;k++){				
	 ContigSorting* contigobj=(*contig_ptr_array)[k];
	 long contig2=contigobj->
	     getrecordnum(contigobj->recordseps,contigobj->number_of_contigs,
			  contigobj->total_conitg_length,
			  957550);
	 cerr << "HALOOOOOOOOOOOOOOOOO : "<< contig2 << endl;
	
    }
*/
    return (0);
}



//////////////////////////////// Transformer Class //////////////////////////////////////////

Transformer::Transformer(char* datFile, unsigned int dim):inDatFile(datFile), dimensions(dim){
    /*
    try{
	//workingDir = get_path(inDatFile);
	//workingDir =".";
    }
    catch (exception& e)
    {
	cout << "Standard exception: " << e.what() << endl;
	workingDir = "";
    }
    */
    stat_dim=dimensions-1;
    overlap_ratio=0;

    for(long i=0;i<dimensions;i++){
	min_len.push_back(2000000000);
    }    
    overlap_function_type=1;
    overlap_flag=0;
    shiftvalue=3;
    SynFilePrefix.append(inDatFile);
    SynFilePrefix.append(".syn.");

    RepFile.append(inDatFile);
    RepFile.append(".rep.");

    Total_final_score=0;
    filter_rep_ratio=1;
    filterrep_flag=0;

   
    draft_flag=0;
    contig_ptr_array=NULL;


}

Transformer::~Transformer(){
//	delete Intervals;
//	delete Chain;
//	delete permutations;
    combinations.clear();
//    for(long i=0;i<combinations.size();i++){
//	combinations[i]
    //   }
} 



////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
// Function: quick sort list
// 1- creats two arrays sorted_list and point_type
// the array sorted_list store indices for the sorted start/end points of 
// the fragments and the array point_type store the corresponding 
// point type (start/end).
// 2- the functions sorts the start/end points w.r.t. the first dimension
// which is stored in blocks[i].region[dimension]. If two points have the
// same first dimension value, then they are recursively sorted w.r.t. the
// second dimension and so on.
// the sorted indices are in the two arrays sorted_list with their corresponding
// point type
//////////////////////////////////////////////////////////////////////

int Transformer::quick_sort_list(long in_dimension){
    
    // Sort w.r.t. the first dimension which is stored at the end of
    // array region
    quicksortkdlistkeysOndim(in_dimension,0,Points.size()-1);
    
    /// sort according to other directions recursive
    //   recursive_subdim_sort(in_dimension,1,Points.size()-1);
    return(1);
    
}
////////////////////////////////////////////////////////////////////////////////
// Recursive sorting
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
// Function recursive_subdim_sort
// sorts the indices w.r.t. the next dimension if they are equal w.r.t. the current one
// ////////////////////////////////////////////////////////////////////////////

void Transformer::recursive_subdim_sort(long in_dimension, long startkeys, long endkeys){
    long i;
    long end=0;
    long start=0;
    unsigned char rec_dimension;
    
    if(in_dimension>dimensions) rec_dimension=1;
    else rec_dimension=in_dimension+1;
    //in_dimension=in_dimension%(dimensions+1)+1;
    char flag=0;
    for(i=startkeys;i<endkeys;i++){	
	if(get_dimension_value(in_dimension,Points[i],i)==get_dimension_value(in_dimension, Points[i+1],(i+1))){
	    if(flag==0){
		start=i;
		flag=1;
	    }
	    end=i+1;
	    if(end<endkeys)continue;
	}
	if((flag==1)&&(end>start)&&(rec_dimension<(dimensions+1))){
	    quicksortkdlistkeysOndim(rec_dimension,start,end);
	    recursive_subdim_sort(rec_dimension,start,end);
	}
	flag=0;
    }
}
////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
// Function quicksortkdlistkeysOndim
// starts sorting w.r.t. the given dimension (dim) between the two array
// boundaries lo and hi
//////////////////////////////////////////////////////////////////////

void Transformer::quicksortkdlistkeysOndim(long dim,long lo, long hi)
{
    //  lo ist der unterste Index, hi ist der oberste Index des
    //  zu sortierenden (Teil-)Feldes a
  
    long i=lo, j=hi;
    //kdMasterPoint temp;
    long temp;
    char tempchar;
    long x=this->get_dimension_value(dim,Points[(long)(((double)lo+(double)hi)/2)],(long)(((double)lo+(double)hi)/2));  
    
    //  Aufteilung
    do
    {    
        while (this->get_dimension_value(dim,Points[(long)i],i)<x) i++;
        while (this->get_dimension_value(dim,Points[(long)j],j)>x) j--;
        
	if (i<=j){
            temp=Points[(unsigned int)i]; Points[(unsigned int)i]=Points[(unsigned int)j]; Points[(unsigned int)j]=temp;
	    tempchar=Types[(unsigned int)i]; 
	    Types[(unsigned int)i]=Types[(unsigned int)j]; 
	    Types[(unsigned int)j]=tempchar;
	    i++; 
	    j--;
        }
    } while (i<=j);
    
    // Rekursion
    if (lo<j) quicksortkdlistkeysOndim(dim,lo, j);
    if (i<hi) quicksortkdlistkeysOndim(dim,i, hi);    
}

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
// Function get_dimension_value
//  returns the coordinate of a point referenced by both 
//  the index (in sorted_list array) and its type (in point_type array) according to the
//  given dimension value.
//////////////////////////////////////////////////////////////////////

long Transformer::get_dimension_value(long dim, long index, long type_index){
    
    // modified to be as blocks data structure
    // Chek if the index is 
    if(Types[type_index]==1){
	return(Chain[index].b.coordinates[dim-1]);
    }
    else if(Types[type_index]==2){
	//long pos=Chain[index].b.coordinates[dim-1]+long (overlap_ratio*(float)(Chain[index].e.coordinates[dim-1]-Chain[index].b.coordinates[dim-1]));
	long pos=get_position_v_point(index,dim-1);
	return(pos);
    }
    else{
	return(Chain[index].e.coordinates[dim-1]);
    }
    
}

/////////////////////////////////////////////////////////////////////
// Handling input file /////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////

int Transformer::parseDatFilekGenomes(char* infile){

ifstream in (infile, ios::in);

if (!in) {
cerr << "Error in opening input file!\n";
exit (1);
}



vector<long> csb(dimensions); vector<long> cse(dimensions);

//parse infile (datfile)    
//bool isnewx, isnewy, isnewz;
 vector<bool> forwardMatches(dimensions); //direction for each dimension/genome 
 string currentline;  
// int flag=0;    
 vector< string> words; 
 
 for(long i=0;i<dimensions;i++){
     forwardMatches[i] = 1;
     csb[i]=0;
     cse[i]=0;
 }
 Point b1(csb); 
 Point e1(cse);
 Interval iv(b1,e1,forwardMatches);
 iv.intervalID=0;
 iv.s=0;
 iv.w=0;
 Intervals.push_back(iv);     
 Chain.push_back(iv);

//  while(infilestream.eof()==0){   
 string zeile;
 string testcomment;
 while (getline (in, zeile)) {
     if (zeile.substr (0, 1).compare ("#") == 0) {
	 fsplit(zeile, " ", words);
	 testcomment=words[words.size()-1];
	 for(unsigned long i=0;i<dimensions;i++){
	     if((testcomment[i]=='p')||(testcomment[i]=='+')){
		 forwardMatches[i] = 1;
	     }
	     else{
		 forwardMatches[i] = 0;
	     }
	 }	
	 continue;
	 
     }
     if (zeile.substr (0, 1).compare ("\n") == 0) {
	 continue;
     }  
     if (zeile.substr (0, 1).compare ("\t") == 0) {
	 continue;
     } 
     if (zeile.substr (0, 1).compare ("") == 0) {
	 continue;
     } 
     if (zeile.substr (0, 1).compare (" ") == 0) {
	 continue;
     } 
//fsplit(zeile, "\t", words);
     char* line_str=(char*)zeile.c_str();
     char* token=strtok(line_str," \t");
     csb[0]=atol(token);
     for(long i=1;i<dimensions;i++){
	 token=strtok(NULL," \t");
	 if(forwardMatches[i]==1){
	     csb[i]=atol(token);
	 }
	 else{
	     cse[i]=atol(token);
	 }
         //csb[i]=atol(words[i].c_str());
     }    	
     getline (in, zeile);
     //fsplit(zeile, "\t", words);
     line_str=(char*)zeile.c_str();
     token=strtok(line_str," \t");
     cse[0]=atol(token);
     for(long i=1;i<dimensions;i++){
     //cse[i]=atol(words[i].c_str());
	 token=strtok(NULL," \t");
	 if(forwardMatches[i]==1){
	     cse[i]=atol(token);
	 }
	 else{
	     csb[i]=atol(token);
	 }
     }
     Point b(csb); 
     Point e(cse);
     int origin_flag=0;
     for(long i=0;i<dimensions;i++){
	 if(csb[i]!=0){
	     origin_flag=1;
	     break;
	 }
	 if(cse[i]!=0){
	     origin_flag=1;
	     break;
	 }	 
     }
 
         
     if(origin_flag==0)continue;
     Interval iv(b,e,forwardMatches);
     
     for(long i=0;i<dimensions;i++){
	 if(cse[i]-csb[i])
	     min_len[i]=min(min_len[i],(cse[i]-csb[i]));
     }
     iv.intervalID=0;
     iv.w=0;
     for(long i=0;i<dimensions;i++){
	 iv.w=iv.w+(cse[i]-csb[i]+1);
     }    
     iv.s=iv.w;    
     iv.choosen_flag=1;
     Intervals.push_back(iv);     
     Chain.push_back(iv);
     //   cout << "start points: "<<csb[0] << " " << csb[1] <<  " " << csb[2] << " " << csb[3] << endl;
     //cout << "end points: "<<cse[0] << " " << cse[1] << " " << cse[2] << " " << cse[3] << endl;

 }
 in.close();
 for(unsigned long i=0;i<Intervals.size();i++){
     Intervals[i].intervalID=i;
     Chain[i].intervalID=i; // id w.r.t. the input intervals
     Intervals[i].index_perc=0;
     Chain[i].index_perc=0;
 }

 cout<<"# Number of input regions = "<<Intervals.size()<<endl;
 cout <<"# Min. region length in sequences 1..k: ";
 for(long i=0;i<dimensions;i++){     
     cout <<  min_len[i]<< " ";
 }
 cout << endl;
 return 0;
}


void Transformer::fsplit (string line, const char* sep, vector<string> & words) {
		
	string::size_type a = 0, e;
	while ((a = line.find_first_not_of (sep, a)) != string::npos) {
		e = line.find_first_of (sep, a);
		if (e != string::npos) {
			words.push_back (line.substr (a, e-a));
			a = e + 1;
		}
		else {
			words.push_back (line.substr (a));
			break;
		}
	}
}

//////////////////////////////////////////////////////////////////////////////
//Filtering Repeats Functions
////////////////////////////////////////////////////////////////////////////
int Transformer::filter_repeats_all_dimensions(){
    // indexing starts from 1...to k
    
    // cout << "# Chaining w.r.t. dimension 1" << endl;
    //chainOneDirectionWithOverlap(1);
    //return (0);
	
 
    for(long j=0;j<dimensions;j++){
	//cout << "# Chaining w.r.t. dimension " << j+1 << endl;
	//cout << "# Overlap flag= " << overlap_flag<< endl;
	filter_repeat(j+1);	
	//cout << "New Chain Size: " << NewChain.size() << endl;
	Chain.clear();	
	// add origing point
	add_origin_to_chain();
	for(unsigned long i=0;i<NewChain.size();i++){
	    NewChain[i].w=NewChain[i].e.coordinates[j]-NewChain[i].b.coordinates[j]+1;
	    NewChain[i].s=NewChain[i].w;
	    //cout << "New Chain score: " << NewChain[i].s << endl;
	    NewChain[i].index_perc=0;
	    Chain.push_back(NewChain[i]);	    
	}
	NewChain.clear();	
	//cout << ", new chain size  "<<  NewChain.size() << endl;
    }
    //report_permutations();
    return 0;
}

long Transformer::filter_repeat(long in_dimension){

    //Apply one dimensional chaining over the Chain list
    long dimension=in_dimension-1;
    // construct point and tyoe array 
    //cout << "The Chain size: " << Chain.size() << endl;
    Points.clear();
    Types.clear();
    for(unsigned long i=0;i<Chain.size();i=i+1){
	//Points.push_back(Chain[i].b.coordinates[dimension]);
	Points.push_back(i);
	Types.push_back(1); // consider start points only	
    }
    
    // sort points w.r.t. 
    quick_sort_list(in_dimension);
    //cout << "Sorting  done. Chaining starts.....\n";

    //long max_interval_index_in_chain=0;
    //float max_score=0;
    for(unsigned long i=1;i<Chain.size();i++){
	Chain[i].choosen_flag=1;
    }
    //cout << "Point Sizeeeeeeeeee "<<Points.size() <<"\n";
    //for(long i=0;i<Points.size();i++){
    //	cout << "Sorted on dimension " << i<< " " << Chain[Points[i]].b.coordinates[dimension]<< endl;
    //   }

    // report sorted list
    for(unsigned long i=1;i<Points.size();i++){	
	if(Types[i]==1){ // it is a start point	    
	    //   check amount of overlap
	    float dif=Chain[Points[i-1]].e.coordinates[dimension]-Chain[Points[i]].b.coordinates[dimension]+1;
	    dif=min(dif,Chain[Points[i]].e.coordinates[dimension]-Chain[Points[i]].b.coordinates[dimension]+1);
	    	    	    
	    float ratio1=dif/double(Chain[Points[i]].e.coordinates[dimension]-Chain[Points[i]].b.coordinates[dimension]+1);
	    float ratio2=dif/double(Chain[Points[i-1]].e.coordinates[dimension]-Chain[Points[i-1]].b.coordinates[dimension]+1);
	    //cout << "dif"<<dif << " R1 "<<ratio1 << " R2 " << ratio2 << " pos "<< Chain[Points[i]].b.coordinates[dimension] <<"\n";
	    if(ratio1 >=(filter_rep_ratio)){
		//if(ratio1 >0.5){
		//cout << "Ratio1: " << ratio1 << endl;
		Chain[Points[i]].choosen_flag=0;
		Intervals[Chain[Points[i]].intervalID].choosen_flag=0;
	    }
	    if(ratio2 >=(filter_rep_ratio)){
		//if(ratio2 >0.5){
		//cout << "Ratio2: " << ratio2 << endl;
		Chain[Points[i-1]].choosen_flag=0;
		Intervals[Chain[Points[i-1]].intervalID].choosen_flag=0;
	    }
	}
	else if(Types[i]==0){
	}
    }
    
    // report best chain    
    //  cout << "# After repeats are filtered: " << endl;
//    while(index!=0){
    //while(1){	
    long index;
    for(unsigned long i=1;i<Chain.size();i++){
	index=i;
	if(Chain[index].choosen_flag==0)continue;
	//cout << "# " << Chain[index].w<< endl;
	//for(int i=0;i<dimensions;i++){
	//    cout << "["<<Chain[index].b.coordinates[i]<< ","<< Chain[index].e.coordinates[i] << "] ";
	//}
	//cout << Chain[index].intervalID <<endl;
	NewChain.push_back(Chain[index]);
	if(Chain[index].intervalID==0)break;
	index=Chain[index].index_perc;
    }
    for(unsigned long i=1;i<Chain.size();i++){
	Chain[i].choosen_flag=0;
    }
    //cout << endl;
    //cout << endl;
    return 0;  
}


/////////////////////////////////////////////////////////////////////
/// Chaining Functions
///////////////////////////////////////////////////////////////////
int Transformer::chooseFragments(){
    // indexing starts from 1...to k
    
    // cout << "# Chaining w.r.t. dimension 1" << endl;
    //chainOneDirectionWithOverlap(1);
    //return (0);
	
 
    for(long j=0;j<dimensions;j++){
	//cout << "# Chaining w.r.t. dimension " << j+1 << endl;
	//cout << "# Overlap flag= " << overlap_flag<< endl;
	//filter_repeat(j+1);
	
	if(overlap_flag==0){
	    chainOneDirection(j);	    
	}
	else{
	    chainOneDirectionWithOverlap(j+1);
	}
	//cout << "New Chain Size: " << NewChain.size() << endl;
	Chain.clear();	
	// add origing point
	add_origin_to_chain();
	for(unsigned long i=0;i<NewChain.size();i++){
	    //NewChain[i].w=NewChain[i].e.coordinates[j]-NewChain[i].b.coordinates[j]+1;
	    NewChain[i].s=NewChain[i].w;
	    //cout << "New Chain score: " << NewChain[i].s << endl;
	    NewChain[i].index_perc=0;
	    Chain.push_back(NewChain[i]);	    
	}
	NewChain.clear();	
	//cout << ", new chain size  "<<  NewChain.size() << endl;
    }
    // report_permutations();
    for(unsigned long j=0;j<Chain.size();j++){			
	for(unsigned long k=0;k<dimensions;k++){	
	    Total_final_score=Total_final_score+(Chain[j].e.coordinates[k]-Chain[j].b.coordinates[k]+1);
	}
    }
    
    cout << "# Total score: " << Total_final_score << "\n";
    return 0;
}



long Transformer::chainOneDirection(int dimension){

    //Apply one dimensional chaining over the Chain list
    
    // construct point and tyoe array 
    //cout << "The Chain size: " << Chain.size() << endl;
    Points.clear();
    Types.clear();
    for(unsigned long i=0;i<Chain.size();i=i+1){
	Chain[i].s=Chain[i].e.coordinates[dimension]-Chain[i].b.coordinates[dimension]+1;
	Chain[i].w=Chain[i].s;
	//Points.push_back(Chain[i].b.coordinates[dimension]);
	Points.push_back(i);
	Types.push_back(1);
	Points.push_back(i);
	//Points.push_back(Chain[i].e.coordinates[dimension]);
	Types.push_back(0);
	
    }
    
    // sort points w.r.t. 
    quick_sort_list(dimension+1);
    //cout << "Sorting  done. Chaining starts.....\n";

    long max_interval_index_in_chain=0;
    float max_score=0;

    // report sorted list
    for(unsigned long i=0;i<Points.size();i++){
	if(Types[i]==1){ // it is a start point
	    
	    //    cout << "Corrd: "<< Chain[Points[i]].b.coordinates[dimension-1]<< " "<<Types[i]<< " "
	    //	 << Points[i]<<" score: "<< Chain[Points[i]].s << endl;
	    Chain[Points[i]].index_perc=max_interval_index_in_chain;
	    Chain[Points[i]].s=Chain[Points[i]].s+max_score;
	}
	else if(Types[i]==0){

	    //  cout << "Coord: " << Chain[Points[i]].e.coordinates[dimension-1]<< " "<<Types[i]<< " "
	    //	 << Points[i]<<" score: "<< Chain[Points[i]].s <<endl;
	    
	    if(Chain[Points[i]].s>max_score){
		max_score=Chain[Points[i]].s;
		max_interval_index_in_chain=Points[i];		
	    }
	}
    }
    
    // report best chain    
    long index=max_interval_index_in_chain;
    //cout << "# Total Score: " << Chain[index].s<< endl;
    cout << "# Total Score in dimension: " << dimension<< " = "<< Chain[index].s<< endl;
    //cout << "# Total Score in dimension: " << dimension<< " = "<< max_score<< endl;
    while(index!=0){
	//while(1){	
	//cout << "# " << Chain[index].w<< endl;
	//cout << "# " << Chain[index].s;
	//for(int i=0;i<dimensions;i++){
	//      cout << " ["<<Chain[index].b.coordinates[i]<< ","<< Chain[index].e.coordinates[i] << "] ";
	//}
	//cout << endl;
//	cout << Chain[index].intervalID <<endl;
	NewChain.push_back(Chain[index]);
	if(Chain[index].intervalID==0)break;
	index=Chain[index].index_perc;
    }
    //  cout << endl;
    //cout << endl;
    return 0;  
}





void Transformer::add_origin_to_chain(){
    vector<long> csb(dimensions); vector<long> cse(dimensions);
    vector<bool> forwardMatches(dimensions); //direction for each dimension/genome 
 
    for(long i=0;i<dimensions;i++){
	forwardMatches[i] = 1;
	csb[i]=0;
	cse[i]=0;
    }
    Point b1(csb); 
    Point e1(cse);
    Interval iv(b1,e1,forwardMatches);
    iv.intervalID=0;
    iv.s=0;
    iv.w=0;  
    Chain.push_back(iv);
}

void Transformer::report_permutations(){
 
    cout << endl<< "# Permutaions in terms of chain id's in input file: " << endl<< endl;

    Points.clear();
    Types.clear();
//    Chain.clear();
    for(unsigned long j=0;j<Chain.size();j=j+1){
	
	Points.push_back(j);
	Types.push_back(1);
	//Points.push_back(j);
	
	//Types.push_back(0);
    }
    for(long i=0;i<dimensions;i++){
       // sort w.r.t. genome and report ids, we will report its index later

       // sort points w.r.t. 
       quick_sort_list(i+1);
       cout << "# Genome " << i+1 << ": ";
       for(unsigned long k=0;k<Points.size();k++){
	   if(Types[k]==1){ // it is a start point
	       
	       // rport chromosome boundary
	       if((draft_flag)&&(k>0)){
		   ContigSorting* contigobj=(*contig_ptr_array)[i];
		   long contig1=contigobj->
		       getrecordnum(contigobj->recordseps,contigobj->number_of_contigs,
				    contigobj->total_conitg_length,
				    Chain[Points[k-1]].b.coordinates[i]);
		   long contig2=contigobj->
		       getrecordnum(contigobj->recordseps,contigobj->number_of_contigs,
				    contigobj->total_conitg_length,
				    Chain[Points[k]].b.coordinates[i]);
		   if(contig1!=contig2){
		       cout << ">"<< contig1 <<","  << contig2<< "< ";
		   }
	       }
	       // report direction
	       if(Chain[Points[k]].isForward[i]==1){
		   cout << "+";
	       }
	       else{
		   cout << "-";
	       }
	       cout << Chain[Points[k]].intervalID+1<<" ";	       
	   }
	   else if (Types[k]==0){	       
	       //cout << "Coord: " << Chain[Points[i]].e.coordinates[dimension-1]<< " "<<Types[i]<< " "
	       //	 << Points[i]<<" score: "<< Chain[Points[i]].s <<endl;	       	       
	   }
       }
       cout << endl;
   }
    cout << endl;
    
    cout << "\n# Permutations w.r.t. identity permutation: "<< endl<< endl;
    // store k-tuple of ids;
    vector <Point> ktuple;
    vector <Point> ktuple_dir;
    vector<long> csb(dimensions);
    vector<long> dirp(dimensions);
    vector<long> inverse_index(dimensions);
    vector <Point> ktuple_inverse;
    vector<long> sorted_wrt_x;
    //Point id_s;

    int set_identity_flag=0;
    
    for(long i=0;i<dimensions;i++){
       // sort w.r.t. genome and report ids, we will report its index later
       // sort points w.r.t. 
       quick_sort_list(i+1);
       if(set_identity_flag==0){
	   set_identity_flag=1;
	   long counter=0;
	   for(unsigned long x=0;x < Points.size();x++){
	       if(Types[x]){
		   Chain[Points[x]].interval_identity=counter+1;
		   if(i==0){
		       sorted_wrt_x.push_back(Points[x]);
		   }
		   counter++;
		   csb[0]=Chain[Points[x]].interval_identity;
		   Point id_s(csb);	
		   ktuple.push_back(id_s);

		   dirp[0]=1;	   		   
		   Point dir_s(dirp);
		   ktuple_dir.push_back(id_s);	
		   
		   inverse_index[0]=counter+1;
		   ktuple_inverse.push_back(inverse_index);
	       }
	   }
       }
       cout << "# Genome " << i+1 << ": ";
       for(unsigned long k=0;k<Points.size();k++){
	   if(Types[k]==1){ // it is a start point
	       
	       // report chromosome boundary
	       if((draft_flag)&&(k>0)){
		   ContigSorting* contigobj=(*contig_ptr_array)[i];
		   long contig1=contigobj->
		       getrecordnum(contigobj->recordseps,contigobj->number_of_contigs,
				    contigobj->total_conitg_length,
				    Chain[Points[k-1]].b.coordinates[i]);
		   long contig2=contigobj->
		       getrecordnum(contigobj->recordseps,contigobj->number_of_contigs,
				    contigobj->total_conitg_length,
				    Chain[Points[k]].b.coordinates[i]);
		   if(contig1!=contig2){
		       cout << ">"<< contig1 <<","  << contig2<< "< ";
		   }
	       }
	       /// end report chromosome boundary

	       if(Chain[Points[k]].isForward[i]==1){
		   cout << "+";
		   ktuple_dir[k].coordinates[i]=1;
	       }
	       else{
		   cout << "-";
		   ktuple_dir[k].coordinates[i]=0;
	       }
	       cout << Chain[Points[k]].interval_identity<<" ";	       
	       if(1){
		   ktuple[k].coordinates[i]=Chain[Points[k]].interval_identity;
		   ktuple_inverse[Chain[Points[k]].interval_identity-1].coordinates[i]=k;
	       }
	   }
	   else if(Types[k]==0){
	       
	       //cout << "Coord: " << Chain[Points[i]].e.coordinates[dimension-1]<< " "<<Types[i]<< " "
	       //	 << Points[i]<<" score: "<< Chain[Points[i]].s <<endl;	       	       
	   }
       }
       cout << endl;
   }

    cout << endl;

    report_compact_chain(&ktuple_inverse);

    return;

    cout << "\n# Compact permutations w.r.t. identity permutation: "<< endl<< endl;
   
    
    
    long* selected_tuples=new long[ktuple.size()];
    for(unsigned long i=0;i<ktuple.size();i++){
	selected_tuples[i]=0;
    }
    selected_tuples[0]=1;
//    cout << "Inverse: "<<endl;
    for(unsigned long i=1;i<ktuple.size();i++){
	//int colinearity_flag=0;
	for(long k=0;k<dimensions;k++){
	    //cout << ktuple_inverse[i].coordinates[k] << " ";
	    //cout << ktuple[i].coordinates[k] << " ";
	    //cout << ktuple_dir[ktuple_inverse[i].coordinates[k]].coordinates[k]<< " ";
	    if(abs(ktuple_inverse[i].coordinates[k]-
	       ktuple_inverse[i-1].coordinates[k])!=1){		
		//colinearity_flag=1;
		selected_tuples[i]=1;		
		//break;
	    }
	    else{
		if(ktuple_dir[ktuple_inverse[i].coordinates[k]].coordinates[k]!=
		    ktuple_dir[ktuple_inverse[i-1].coordinates[k]].coordinates[k]){
		    //colinearity_flag=1;
		    selected_tuples[i]=1;		
		    //break;
		}
		else{ // to check for transposition
		    if((ktuple_dir[ktuple_inverse[i].coordinates[k]].coordinates[k]==1)
		       &&((ktuple_inverse[i].coordinates[k]-ktuple_inverse[i-1].coordinates[k])<0)){
			//colinearity_flag=1;
			selected_tuples[i]=1;		
			//break;
		    }
		    if((ktuple_dir[ktuple_inverse[i].coordinates[k]].coordinates[k]==0)
		       &&((ktuple_inverse[i].coordinates[k]-ktuple_inverse[i-1].coordinates[k])>0)){
			//colinearity_flag=1;
			selected_tuples[i]=1;		
			//break;
		    }
		}
		// tp check for chromosome boundary
	       if((draft_flag)&&(i>0)){
		   long contig1, contig2;
		   
		   
		   ContigSorting* contigobj=(*contig_ptr_array)[k];
		   contig1=contigobj->
		       getrecordnum(contigobj->recordseps,contigobj->number_of_contigs,
				    contigobj->total_conitg_length,
				    Chain[sorted_wrt_x[ktuple_inverse[i-1].coordinates[k]]].b.coordinates[k]);
		   contig2=contigobj->
		       getrecordnum(contigobj->recordseps,contigobj->number_of_contigs,
				    contigobj->total_conitg_length,
				    Chain[sorted_wrt_x[ktuple_inverse[i].coordinates[k]]].b.coordinates[k]);		   

		   if(contig1!=contig2){
		       //cout << ">"<< contig1 <<","  << contig2<< "< ";
		       selected_tuples[i]=1;
		   }
		   /*
		   cout << "HHHHHHHHHHI>"<< contig1 <<","  << contig2<< " id: "
		       << sorted_wrt_x[ktuple_inverse[i-1].coordinates[k]]<<" "<<
		       sorted_wrt_x[ktuple_inverse[i].coordinates[k]]<< 
		       " coord "
			<<Chain[sorted_wrt_x[ktuple_inverse[i-1].coordinates[k]]].b.coordinates[k]<<" "<<
		       Chain[sorted_wrt_x[ktuple_inverse[i].coordinates[k]]].b.coordinates[k]<< "\n";
		   cout << "KTUPLE>"
			<<ktuple_inverse[i-1].coordinates[k]<<" "<<
		       ktuple_inverse[i].coordinates[k]<< "\n";
		   */
	       }
	    }
	}
	//if(selected_tuples[i]==1)cout << " selected";
	//cout << endl;
    }
    //cout << endl;
    long* rename_array=new long[ktuple.size()];
    long counter=0;
    for(unsigned long i=0;i<ktuple.size();i++){
	if(selected_tuples[ktuple[i].coordinates[0]-1]==1){
	    rename_array[i]=counter;
	    counter++;
	}
	else{
	    rename_array[i]=0;
	}
    }
    
    for(long k=0;k<dimensions;k++){
	cout << "# Genome " << k+1 << ":  ";
	for(unsigned long i=0;i<ktuple.size();i++){
	    if(selected_tuples[ktuple[i].coordinates[k]-1]==1){
		// report chromosome boundary
	       if((draft_flag)&&(i>0)){
		   ContigSorting* contigobj=(*contig_ptr_array)[k];
		   long contig1=contigobj->
		       getrecordnum(contigobj->recordseps,contigobj->number_of_contigs,
				    contigobj->total_conitg_length,
				    Chain[sorted_wrt_x[ktuple_inverse[i-1].coordinates[k]]].b.coordinates[k]);
		   long contig2=contigobj->
		       getrecordnum(contigobj->recordseps,contigobj->number_of_contigs,
				    contigobj->total_conitg_length,
				    Chain[sorted_wrt_x[ktuple_inverse[i].coordinates[k]]].b.coordinates[k]);
		   if(contig1!=contig2){
		       cout << ">"<< contig1 <<","  << contig2<< "< ";
		       cout << ">"<< Chain[sorted_wrt_x[ktuple_inverse[i-1].coordinates[k]]].b.coordinates[k] 
			    <<","  << Chain[sorted_wrt_x[ktuple_inverse[i].coordinates[k]]].b.coordinates[k]<< "< ";
		   }
	       }
	       /// end report chromosome boundary

		if(ktuple_dir[i].coordinates[k]){
		    cout << "+" ;
		    cout << ktuple[i].coordinates[k]+1<< " ";
		    //cout << rename_array[ktuple[i].coordinates[k]-1]+1<< " ";
		}
		else{
		    cout << "-";
		    cout << ktuple[i].coordinates[k]<< "  ";
		    //cout << rename_array[ktuple[i].coordinates[k]-1]+1<< " ";
		}
	    }
	}
	cout << endl;
    }
    
    delete[] rename_array;
    delete[]selected_tuples;

    



  
}

/////////////Functions for handling Overlaps /////////////////

long Transformer::chainOneDirectionWithOverlap(long in_dimension){

    long dimension=in_dimension-1;
    // allowed overlap
    //float overlap_ratio=0.5;

    //Apply one dimensional chaining over the Chain list
    
    // construct point and tyoe array 
    //cout << "The Chain size: " << Chain.size() << endl;
    Points.clear();
    Types.clear();
    for(unsigned long i=0;i<Chain.size();i=i+1){
	Chain[i].s=Chain[i].e.coordinates[dimension]-Chain[i].b.coordinates[dimension]+1;
	Chain[i].w=Chain[i].s;
	// start point
	Points.push_back(i);
	Types.push_back(1);
	// intermediate v point
	Points.push_back(i);
	Types.push_back(2);
	// end point
	Points.push_back(i);
	//Points.push_back(Chain[i].e.coordinates[dimension]);
	Types.push_back(0);
	
    }
    
    // sort points w.r.t. 
    quick_sort_list(dimension+1);
    //cout << "Sorting  done. Chaining starts.....\n";
    
    // create 1D kdTree to store for end points of intervals updated with certain
    // priority
    //My kdTree code  
    kdTree kdtreeobj;
    long numofblocks=Chain.size();
    //Kdelem_ptr *perm;    //stores the elements of the tree in order
    //Kdelem_ptr nntarget, first, tail; 
    //Kdnode_ptr *bucket_ptr;  //stores the bucket location of each element
    //Kdnode_ptr root;
    //Kdelem_ptr blocks;   //permanent storage of the elements

    //blocks = (Kdelem_ptr) malloc(numofblocks* sizeof(Kdelem));
    blocks=new Kdelem[numofblocks];
    for(long i=0;i<numofblocks;i++){
	blocks[i].region=(Region*)malloc(sizeof(Region));
	if(blocks[i].region==NULL){
	    cerr<<"Not Enough Memory \n";
	    exit(-3);
	}	
	blocks[i].region[0].start=Chain[i].b.coordinates[dimension];	
	blocks[i].region[0].end=Chain[i].e.coordinates[dimension];
	blocks[i].weight=Chain[i].e.coordinates[dimension]-Chain[i].b.coordinates[dimension]+1;
	blocks[i].globalscore=(double)blocks[i].weight-blocks[i].region[0].end;
	//cout << "Blocks global score: "<< ((float)blocks[i].globalscore) << endl;
	if(i>0){
	    blocks[i-1].next = &blocks[i];
	}
    }
    

    bucket_ptr = (Kdnode_ptr *) malloc(numofblocks * sizeof(Kdnode_ptr));
    if(bucket_ptr==NULL)
    {
	printf("Not Enough Memory\n");
	return(-1);
    }    
    
    perm = (Kdelem_ptr *) malloc(numofblocks * sizeof(Kdelem_ptr));
    if(perm==NULL)
    {
	printf("Not Enough Memory\n");
	return(-1);
    } 
    for(long i=0;i<numofblocks;i++){
	perm[i] = &blocks[i];
    }
 
    kdtreeobj.dimensions=1;
    
    kdtreeobj.blocks=blocks;
    kdtreeobj.bucket_ptr=bucket_ptr;
    kdtreeobj.perm=perm;
    kdtreeobj.numofblocks=numofblocks;
    kdtreeobj.root = kdtreeobj.build(0, numofblocks-1, 0);
    if(kdtreeobj.root==NULL){
	printf("Not enough memory for building a kdtree \n");
	return(-1);
    }
    for(long i=0;i<numofblocks;i++)
    {
	perm[i]->num = i; //setting num according to their order
	perm[i]->start_child_list=NULL;
    }

    kdtreeobj.DeleteAll();

    Kdelem_ptr lower_bound=new Kdelem;
    Kdelem_ptr upper_bound=new Kdelem;
    lower_bound->region=(Region*)malloc(sizeof(Region)*(1));
    if(lower_bound->region==NULL){
      return(-1);
    }
    upper_bound->region=(Region*)malloc(sizeof(Region)*(1));
    if(upper_bound->region==NULL){
      return(-1);
    }
    // for testing
    /*
    kdtreeobj.gapped_undelete(&blocks[Points[0]]);
    kdtreeobj.gapped_undelete(&blocks[Points[3]]);
    lower_bound->region[0].start=507776;
    lower_bound->region[0].end=lower_bound->region[0].start;
    upper_bound->region[0].start=507865;
    upper_bound->region[0].end=upper_bound->region[0].start;
    long resindex=kdtreeobj.region_gapped_nearestNeighbor(upper_bound,lower_bound);
    cout << "Perm index "<<resindex << endl;
    cout << "Corrd: "<< Chain[perm[resindex]-blocks].b.coordinates[0]<< " "
	 <<" score: "<< Chain[perm[resindex]-blocks].s << endl <<" GlobalScoreBlock: "
	 << perm[resindex]->globalscore <<endl;
    //return 0;
    //end for testing
    */
    long max_interval_index_in_chain=0;
    float max_score=0;

    

    
    long permindex;
    float score1;
    for(unsigned long i=0;i<Points.size();i++){
	if(Types[i]==1){ // it is a start point
	    
	    //    cout << "Corrd: "<< Chain[Points[i]].b.coordinates[dimension-1]<< " "<<Types[i]<< " "
	    //	 << Points[i]<<" score: "<< Chain[Points[i]].s << endl;
	    Chain[Points[i]].index_perc=max_interval_index_in_chain;
	    Chain[Points[i]].s=Chain[Points[i]].w+max_score;
	}
	else if(Types[i]==2){ // it is the overlap limit point
	    lower_bound->region[0].start=Chain[Points[i]].b.coordinates[dimension];
	    lower_bound->region[0].end=lower_bound->region[0].start;
	    //long pos=Chain[Points[i]].b.coordinates[dimension]+
	    //	long (overlap_ratio*(float)(Chain[Points[i]].e.coordinates[dimension]-
	    //				    Chain[Points[i]].b.coordinates[dimension]));
	    long pos=get_position_v_point(Points[i],dimension);
	    upper_bound->region[0].start=pos;
	    upper_bound->region[0].end=upper_bound->region[0].start;
	    //permindex = kdtreeobj.gapped_nearestNeighbor(&blocks[Points[i]]);
	    permindex=kdtreeobj.region_gapped_nearestNeighbor(upper_bound,lower_bound);
	    //cout << "Permindex: " << permindex << endl;
	    blocks[Points[i]].globalscore=Chain[Points[i]].s-blocks[Points[i]].region[0].end;
	    if((permindex>=0)){
		long index_in_chain=(perm[permindex]-blocks);
		//cout << "Index in chain: " << index_in_chain << endl;
		//cout << "Start: " << Chain[index_in_chain].b.coordinates[dimension]<< " "
		//    <<Chain[Points[i]].b.coordinates[dimension] << endl;
		//cout << "Score: " << kdtreeobj.nndist << endl;
	
	
		float val=(float)(Chain[index_in_chain].e.coordinates[dimension]-
				  Chain[Points[i]].b.coordinates[dimension])/
		    (float)(Chain[index_in_chain].e.coordinates[dimension]-
			    Chain[index_in_chain].b.coordinates[dimension]);
		if((Chain[index_in_chain].b.coordinates[dimension]<Chain[Points[i]].b.coordinates[dimension])&&
			   (val<overlap_ratio)){
		    score1=(float)kdtreeobj.nndist+Chain[index_in_chain].e.coordinates[dimension];
		    // subtract overlap
		    //  cout << "Score1: " << score1 <<endl;
		    // cout << " Index in Chain " << index_in_chain << endl;
		    score1=score1+Chain[Points[i]].w
			-(Chain[index_in_chain].e.coordinates[dimension]-Chain[Points[i]].b.coordinates[dimension]+1);  
		    // cout << "Score: " << score1 << " stored score: "<< Chain[Points[i]].s << endl;
		    if(score1>Chain[Points[i]].s){
			Chain[Points[i]].s=score1;
			Chain[Points[i]].index_perc=index_in_chain;
			// cout << "PrecedingIndex: " << Chain[Points[i]].index_perc << endl;
		    }
		    blocks[Points[i]].globalscore=Chain[Points[i]].s-blocks[Points[i]].region[0].end;
		    continue;
		}
		else{
		    
		    long permindex2=permindex;
		    long index_in_chain2=index_in_chain;
		    vector<long> list;
		    while(permindex2>0){
			list.push_back(index_in_chain2);
			kdtreeobj.DeleteElem(&blocks[index_in_chain2]);  
			permindex2=kdtreeobj.region_gapped_nearestNeighbor(upper_bound,lower_bound);
			if(permindex2<=0)break;
			index_in_chain2=(perm[permindex2]-blocks);
			//index_in_chain=(perm[permindex2]-blocks);
			float val=(float)(Chain[index_in_chain2].e.coordinates[dimension]-
					  Chain[Points[i]].b.coordinates[dimension])/
			    (float)(Chain[index_in_chain2].e.coordinates[dimension]-
				    Chain[index_in_chain2].b.coordinates[dimension]);
			
			if((Chain[index_in_chain2].b.coordinates[dimension]<Chain[Points[i]].b.coordinates[dimension])&&
			   (val<overlap_ratio)){
			    score1=(float)kdtreeobj.nndist+Chain[index_in_chain2].e.coordinates[dimension];
			    // subtract overlap
			    //  cout << "Score1: " << score1 <<endl;
			    // cout << " Index in Chain " << index_in_chain << endl;
			    score1=score1+Chain[Points[i]].w
				-(Chain[index_in_chain2].e.coordinates[dimension]-Chain[Points[i]].b.coordinates[dimension]+1);  
			    // cout << "Score: " << score1 << " stored score: "<< Chain[Points[i]].s << endl;
			    if(score1>Chain[Points[i]].s){
				Chain[Points[i]].s=score1;
				Chain[Points[i]].index_perc=index_in_chain2;
				// cout << "PrecedingIndex: " << Chain[Points[i]].index_perc << endl;
			    }
			    blocks[Points[i]].globalscore=Chain[Points[i]].s-blocks[Points[i]].region[0].end;
			    break;
			}
		    }
		    for(unsigned long x=0;x<list.size();x++){
			kdtreeobj.gapped_undelete(&blocks[list[x]]);  
		    }
		    
		}
	    }
	}
	else{

	    //  cout << "Coord: " << Chain[Points[i]].e.coordinates[dimension-1]<< " "<<Types[i]<< " "
	    //	 << Points[i]<<" score: "<< Chain[Points[i]].s <<endl;
	    // check if it make sense to update it
	    kdtreeobj.gapped_undelete(&blocks[Points[i]]);
	    if(Chain[Points[i]].s>max_score){
		max_score=Chain[Points[i]].s;
		max_interval_index_in_chain=Points[i];		
	    }
	}
    }
    
    if(bucket_ptr!=NULL)free(bucket_ptr);
    free_blocks_array();
    free_blocks_array();
    if(perm!=NULL)free(perm); 
    if(lower_bound!=NULL){
	free(lower_bound->region);
	delete lower_bound;lower_bound=NULL;
    }
    if(upper_bound!=NULL){
	free(upper_bound-> region);
	delete upper_bound;upper_bound=NULL;
    }

    // report best chain
    long index=max_interval_index_in_chain;
    cout << "# Total Score in dimension: " << dimension<< " = "<< Chain[index].s<< endl;
    //cout << "# Total Score in dimension: " << dimension<< " = "<< max_score<< endl;
    NewChain.clear();
    while(index!=0){
	//while(1){	
//	cout << "# " << Chain[index].w<< endl;
//	cout << Chain[index].s;
//	for(int i=0;i<dimensions;i++){
//	    cout << " ["<<Chain[index].b.coordinates[i]<< ","<< Chain[index].e.coordinates[i] << "] ";
//	}
//	cout << endl;
//	cout << Chain[index].intervalID <<endl;
	
	// here we can resolve overlap and correct boundary	

	if(NewChain.size()>0){
	    if(Chain[index].e.coordinates[dimension]>=NewChain[NewChain.size()-1].b.coordinates[dimension]){
		//cout << "HALOOOOOOOOOOOOOOOOOOOOOOOO\n";
		//for(int i=0;i<dimensions;i++){
		//   cout << "["<<Chain[index].b.coordinates[i]<< ","<< Chain[index].e.coordinates[i] << "] ";
		//}
		Chain[index].e.coordinates[dimension]=NewChain[NewChain.size()-1].b.coordinates[dimension]-1;
		//cout << "\nAfteeeeeeeeerHALOOOOOOOOOOOOOOOOOOOOOOOO\n";
		//for(int i=0;i<dimensions;i++){
		//   cout << "["<<Chain[index].b.coordinates[i]<< ","<< Chain[index].e.coordinates[i] << "] ";
		//}
		//cout << endl;
	    }
	}
	NewChain.push_back(Chain[index]);
       
	if(Chain[index].intervalID==0)break;
	index=Chain[index].index_perc;
    }
    //  cout << endl;
    // cout << endl;
    return 0;
  

}



void Transformer::free_blocks_array(){
	
//	Kdelem_ptr elem;
//    child* child_ptr;
    unsigned long j;
    if(blocks==NULL)return;
    for(j=0;j<Chain.size();j++){
//	child_ptr=blocks[j].start_child_list;
//	while(child_ptr!=NULL){
//	    blocks[j].start_child_list=blocks[j].start_child_list->next;
//	    free(child_ptr);
//	    child_ptr=blocks[j].start_child_list;	    
//	}
	free(blocks[j].region); blocks[j].region=NULL;
	
    }
//    free(blocks);blocks=NULL;
    delete[] blocks;
    blocks=NULL;
}



long Transformer::get_position_v_point(long index,long dimension){
    long pos;
    
    if(overlap_function_type==1){
	pos=Chain[index].b.coordinates[dimension]+long (overlap_ratio*(float)(Chain[index].e.coordinates[dimension]-Chain[index].b.coordinates[dimension]));
    }
    else{
	pos=Chain[index].b.coordinates[dimension]+ long (min_len[dimension]/2)-1;
    }
    return(pos);
}

void Transformer::report_compact_chain(vector <Point> * ktuple_inverse){

    Points.clear();
    Types.clear();
    for(unsigned long j=0;j<Chain.size();j=j+1){	
	Points.push_back(j);
	Types.push_back(1);
	for(long k=0;k<dimensions;k++){	
	    Total_final_score=Total_final_score+(Chain[j].e.coordinates[k]-Chain[j].b.coordinates[k]+1);
	}
    }
    
    //cout << "Total score: " << Total_final_score << "\n";
    // sort w.r.t. 1st dimension
    quick_sort_list(1);

    //vector <Interval> compact_chain;

    compact_chain.clear();
    compact_chain.push_back(Chain[Points[0]]);
    
    for(unsigned long i=1;i<Chain.size();i++){
	// construct file extension
	char colinearity_flag=1;
	for(long j=0;j<dimensions;j++){
	    //cout << "["<<Chain[index].b.coordinates[i]<< ","<< Chain[index].e.coordinates[i] << "] ";
	    // check polarity
	    if(Chain[Points[i]].isForward[j]!=Chain[Points[i-1]].isForward[j]){
		colinearity_flag=0;
		break;
	    }	
	    // check colinearity in other genomes
	    if(abs((*ktuple_inverse)[i].coordinates[j]-(*ktuple_inverse)[i-1].coordinates[j])!=1){
		colinearity_flag=0;
		break;
	    }	    	    

	}
	 // check colinearity in the same genome (transposition)
	for(long j=0;j<dimensions;j++){	    
	    if(Chain[Points[i]].isForward[j]==1){
		if(Chain[Points[i]].e.coordinates[j]<Chain[Points[i-1]].e.coordinates[j]){
		    colinearity_flag=0;
		    break;
		}
		if(Chain[Points[i]].e.coordinates[j]<Chain[Points[i-1]].e.coordinates[j]){
		    colinearity_flag=0;
		    break;
		}
	    }
	    else{
		if(Chain[Points[i]].b.coordinates[j]>Chain[Points[i-1]].b.coordinates[j]){
		    colinearity_flag=0;
		    break;
		}
		if(Chain[Points[i]].b.coordinates[j]>Chain[Points[i-1]].b.coordinates[j]){
		    colinearity_flag=0;
		    break;
		}
 		 
	    }
	    
	}
	// check for chromsome boundary
	     
	if((draft_flag)&&(i>0)){
	    for(long j=0;j<dimensions;j++){	
		ContigSorting* contigobj=(*contig_ptr_array)[j];
		long contig1=contigobj->
		    getrecordnum(contigobj->recordseps,contigobj->number_of_contigs,
				 contigobj->total_conitg_length,
			     Chain[Points[i-1]].b.coordinates[j]-shiftvalue);
		long contig2=contigobj->
		    getrecordnum(contigobj->recordseps,contigobj->number_of_contigs,
				 contigobj->total_conitg_length,
				 Chain[Points[i]].b.coordinates[j]-shiftvalue);
		//cout << "HALOO "<< contig1 <<","  << contig2 << " "<< Chain[Points[i-1]].b.coordinates[j]-shiftvalue <<
		//   " " << Chain[Points[i]].b.coordinates[j]-shiftvalue << "\n";
		if(contig1!=contig2){
		    //cout << ">"<< contig1 <<","  << contig2<< "< ";
		    colinearity_flag=0;
		    break;
		}
	    }
	}
	      
	// end check chromosome boundary
	if(colinearity_flag==1){ // this polarity colinearlity
	   	    
	    // update boundary in compact_chain
	    
	    for(long j=0;j<dimensions;j++){
		//cout << "["<<Chain[index].b.coordinates[i]<< ","<< Chain[index].e.coordinates[i] << "] ";
		if(Chain[Points[i]].isForward[j]==1){
		    compact_chain[compact_chain.size()-1].e.coordinates[j]=Chain[Points[i]].e.coordinates[j];
		}	
		else{
		    compact_chain[compact_chain.size()-1].b.coordinates[j]=Chain[Points[i]].b.coordinates[j];
		}
	    }
	    
	}
	else{
	    compact_chain.push_back(Chain[Points[i]]);
	}
    }

    cout << "\n# Reporting Optimal Non-compact Chain: \n"; 
    //cout << "# Score "<< Chain[Chain.size()-1].s << "\n"; 
    for(unsigned long i=0;i<Chain.size();i++){
	cout << "# "; 
	for(long j=0;j<dimensions;j++){
	    if(Chain[Points[i]].isForward[j]==1)
		cout << "+ ";
	    else
		cout << "- ";
	}
	cout << endl;
	for(long j=0;j<dimensions;j++){
	    if(Chain[Points[i]].interval_identity!=1){
		cout << "["<< Chain[Points[i]].b.coordinates[j]-shiftvalue << ","<<Chain[Points[i]].e.coordinates[j]-shiftvalue << "] ";
	    }
	    else{
		cout << "["<< Chain[Points[i]].b.coordinates[j] << ","<<Chain[Points[i]].e.coordinates[j] << "] ";
	    }
	}
	cout << endl;
    }

    cout << "\n# Reporting Compact Chain: \n"; 
    //cout << "# Score "<< compact_chain[compact_chain.size()-1].s << "\n"; 
    for(unsigned long i=0;i<compact_chain.size();i++){
	cout << "# "; 
	for(unsigned long j=0;j<dimensions;j++){
	    if(compact_chain[i].isForward[j]==1)
		cout << "+ ";
	    else
		cout << "- ";
	}
	cout << endl;

	for(long j=0;j<dimensions;j++){
	    if(compact_chain[i].b.coordinates[j]>shiftvalue){
		cout << "["<< compact_chain[i].b.coordinates[j]-shiftvalue;
	    }
	    else{
		cout << "["<< compact_chain[i].b.coordinates[j];
	    }
	    if(compact_chain[i].e.coordinates[j]>shiftvalue)
	    {
		cout << ","<<compact_chain[i].e.coordinates[j]-shiftvalue << "] ";
	    }       
	    else{
		cout << ","<<compact_chain[i].e.coordinates[j] << "] ";
	    }
	}
	cout << endl;
    }
    //cout << "###DOOOOOOOOOOOOOOOOON\n";
    report_permutations_from_compact_chain();

}

/// this function creates gnuplot files

void Transformer::report_chain(){
    
    // loop over choosen fragments in chain
    // mark them and report them in the correct file
    // extension
    //string command="rm -f "+workingDir+"syntenies.*";
    string command="rm -f "+SynFilePrefix+"*";
    //cerr << "HALOOOOOOOOOOOOOOO: " << command << endl;
    system(command.c_str());
    string extension="";
    //string synteny_file_name_prefix=workingDir+"syntenies.";
    string synteny_file_name_prefix=SynFilePrefix;

    generate_combinations();
    create_files_for_all_suffixes(synteny_file_name_prefix);
    
    FILE* fptr=NULL;
    for(unsigned long i=0;i<Chain.size();i++){
	// construct file extension
	for(unsigned long j=0;j<dimensions;j++){
	    //cout << "["<<Chain[index].b.coordinates[i]<< ","<< Chain[index].e.coordinates[i] << "] ";
	    if(Chain[i].isForward[j]==1){
		extension.append("p");
	    }
	    else{
		extension.append("m");
	    }	    
	}
	string synteny_file_name=synteny_file_name_prefix+extension;
	fptr=FOPEN((char*)synteny_file_name.c_str(),"r");
	if(fptr==NULL){
	    fptr=FOPEN((char*)synteny_file_name.c_str(),"w");
	    if(fptr==NULL){
		cerr << "Error: Unable to open output file"<< synteny_file_name << endl;
		exit(-1);
	    }
	    fprintf(fptr,"#-- MATCHES %s\n",extension.c_str());
	    for(long j=0;j<dimensions;j++){
		fprintf(fptr,"0\t");
	    }
	    fprintf(fptr,"\n");
	    for(long j=0;j<dimensions;j++){
		fprintf(fptr,"0\t");
	    }
	    fprintf(fptr,"\n\n");
	    fclose(fptr);
	}
	else{
	    fclose(fptr);
	}
	fptr=FOPEN((char*)synteny_file_name.c_str(),"a");
	if(fptr==NULL){
	    cerr << "Error: Unable to open output file"<< synteny_file_name << endl;
	    exit(-1);
	}
	// check origin
	int origin_flag=0;
	for(long j=0;j<dimensions;j++){
	    if(Chain[i].b.coordinates[j]!=0){
		origin_flag=1;
		break;
	    }
	    if(Chain[i].e.coordinates[j]!=0){
		origin_flag=1;
		break;
	    }
	    
	}
	if(origin_flag){
	    report_combinations(Chain[i], synteny_file_name,extension,0);
	    for(long j=0;j<dimensions;j++){
		if(Chain[i].isForward[j]==1){
		    if(Chain[i].b.coordinates[j]<shiftvalue)
			fprintf(fptr,"%lu ",Chain[i].b.coordinates[j]);
		    else fprintf(fptr,"%lu ",Chain[i].b.coordinates[j]-shiftvalue);
		}
		else{
		    if(Chain[i].e.coordinates[j]<shiftvalue)
			fprintf(fptr,"%lu ",Chain[i].e.coordinates[j]);
		    else fprintf(fptr,"%lu ",Chain[i].e.coordinates[j]-shiftvalue);
		}
	    }
	    fprintf(fptr,"\n");
	    for(long j=0;j<dimensions;j++){
		if(Chain[i].isForward[j]==1){
		    if(Chain[i].e.coordinates[j]<shiftvalue)
			fprintf(fptr,"%lu ",Chain[i].e.coordinates[j]);
		    else
			fprintf(fptr,"%lu ",Chain[i].e.coordinates[j]-shiftvalue);
		    
		}
		else{
		    if(Chain[i].b.coordinates[j]<shiftvalue)
			fprintf(fptr,"%lu ",Chain[i].b.coordinates[j]);
		    else
			fprintf(fptr,"%lu ",Chain[i].b.coordinates[j]-shiftvalue);
		}
	    }
	    fprintf(fptr,"\n\n");
	}
	fclose(fptr);
	extension.clear();
	
    }
    for(unsigned long i=0;i<Intervals.size();i++){
	Intervals[i].choosen_flag=0;
    }
    for(unsigned long i=0;i<Chain.size();i++){
	Intervals[Chain[i].intervalID].choosen_flag=1;
    }
    //string repeat_file_name=repeat_file_name_prefix+extension;
    
    // for plotting
    command.clear();
    for(unsigned long j=0;j<combinations.size();j++){
	for(unsigned long i=0;i<extension_array.size();i++){
	    command.clear();
	    string tmp_file=synteny_file_name_prefix+extension_array[i]+"."+combinations[j];
	    FILE* tmp_fptr=FOPEN(tmp_file.c_str(),"a");
	    fprintf(tmp_fptr,"\n");
	    fclose(tmp_fptr);
	    command="cat "+synteny_file_name_prefix+extension_array[i]+"."+combinations[j]+" >> "
		+synteny_file_name_prefix+combinations[j]+".dat";
	    //cerr << command<<endl;
	    system(command.c_str());
	    string x=synteny_file_name_prefix+extension_array[i]+"."+combinations[j];
	    unlink(x.c_str());
	}
	
    }

    
//    command="cat "+workingDir+"syntenies.* > ";
//    command=command+workingDir+"syntenies.dat";

    command="cat "+SynFilePrefix+"p* > ";
    command=command+SynFilePrefix+"dat";
    
    system(command.c_str());
//    command="rm -f "+workingDir+"syntenies.p*";
    command="rm -f "+SynFilePrefix+"p*";

    system(command.c_str());
    generate_gp_files(synteny_file_name_prefix);
}


void Transformer::generate_combinations(){
    string comb="";
    long digits=(long)log10((double)dimensions*2)+2;  
    digits=digits*2+2;
    char* str=new char[digits];
    //cerr << "HALOOOOOOOOOOOOOO " << digits << "\n";
    for(long i=0;i<(dimensions-1);i++){	
	for(long j=i+1;j<dimensions;j++){
	    
	    sprintf(str,"%lux%lu",i+1,j+1);
	    comb.append(str);
	    combinations.push_back(comb);
	    proj1.push_back(i);
	    proj2.push_back(j);
	    comb.clear();
	}
    }
    delete[] str;
}

void Transformer::report_combinations(Interval input_interval, string in_file_name,string extension, int header_flag){

    FILE* fptr;
    string working_file_name;
    for(unsigned long i=0;i<combinations.size();i++){
	working_file_name.clear();
	working_file_name=in_file_name+"."+combinations[i];
	fptr=FOPEN((char*)working_file_name.c_str(),"r");
	if(fptr==NULL){
	    fptr=FOPEN((char*)working_file_name.c_str(),"w");
	    if(fptr==NULL){
		cerr << "Error: Unable to open output file"<< in_file_name << endl;
		exit(-1);
	    }
	    fprintf(fptr,"#-- MATCHES %s\n",extension.c_str());
	    for(long j=0;j<2;j++){
		fprintf(fptr,"0\t");
	    }
	    fprintf(fptr,"\n");
	    for(long j=0;j<2;j++){
		fprintf(fptr,"0\t");
	    }
	    fprintf(fptr,"\n\n");
	    fclose(fptr);
	}
	else{
	    fclose(fptr);
	}
    }
    if(header_flag)return;

    //long digits=(long)log10((double)dimensions*2);
    long digits=(long)log10((double)dimensions*2)+2;  
    digits=digits*2+2;
    char* str=new char[digits];
    for(unsigned long i=0;i<combinations.size();i++){

	working_file_name.clear();
	working_file_name=in_file_name+"."+combinations[i];
	fptr=FOPEN((char*)working_file_name.c_str(),"a");
	if(fptr==NULL){
	    cerr << "Error: Unable to open output file"<< in_file_name << endl;
	    exit(-1);
	}
	string comb=combinations[i];
	//long proj2=atol(comb.c_str()[1]);
	if((input_interval.isForward[proj1[i]]==1)&&(input_interval.isForward[proj2[i]]==1)){
	    fprintf(fptr,"%lu %lu\n",max(0,input_interval.b.coordinates[proj1[i]]-shiftvalue),max(0,input_interval.b.coordinates[proj2[i]]-shiftvalue));
	    fprintf(fptr,"%lu %lu\n\n",max(0,input_interval.e.coordinates[proj1[i]]-shiftvalue),max(0,input_interval.e.coordinates[proj2[i]]-shiftvalue));
	}
	else if((input_interval.isForward[proj1[i]]==1)&&(input_interval.isForward[proj2[i]]==0)){
	    fprintf(fptr,"%lu %lu\n",max(0,input_interval.b.coordinates[proj1[i]]-shiftvalue),max(0,input_interval.e.coordinates[proj2[i]]-shiftvalue));
	    fprintf(fptr,"%lu %lu\n\n",max(0,input_interval.e.coordinates[proj1[i]]-shiftvalue),max(0,input_interval.b.coordinates[proj2[i]]-shiftvalue));
	}
	else if((input_interval.isForward[proj1[i]]==0)&&(input_interval.isForward[proj2[i]]==1)){
	    fprintf(fptr,"%lu %lu\n",max(0,input_interval.b.coordinates[proj1[i]]-shiftvalue),max(0,input_interval.e.coordinates[proj2[i]]-shiftvalue));
	    fprintf(fptr,"%lu %lu\n\n",max(0,input_interval.e.coordinates[proj1[i]]-shiftvalue),max(0,input_interval.b.coordinates[proj2[i]]-shiftvalue));
	}
	else if((input_interval.isForward[proj1[i]]==0)&&(input_interval.isForward[proj2[i]]==0)){
	    fprintf(fptr,"%lu %lu\n",max(0,input_interval.b.coordinates[proj1[i]]-shiftvalue),max(0,input_interval.b.coordinates[proj2[i]]-shiftvalue));
	    fprintf(fptr,"%lu %lu\n\n",max(0,input_interval.e.coordinates[proj1[i]]-shiftvalue),max(0,input_interval.e.coordinates[proj2[i]]-shiftvalue));
	}
	fclose(fptr);
    }
    delete[] str;
}

void Transformer::create_files_for_all_suffixes(string synteny_file_name_prefix){
    
    long no_files=(long)exp2((double)(dimensions-1));
    long buffer;
    //int shift=1;
    //cerr << "No. of files: " << no_files << endl;
   
    for(long i=0;i<no_files;i++){
	buffer=i; 
	//cerr << "buffer: "<< buffer << endl;
	string tmp="p";
	//for(long j=0;j<dimensions-1;j++){
	for(long j=dimensions-2;j>=0;j--){
	    //if(buffer & (1 << j)) {
	    if(buffer & (1 << j)) {
		//cerr << "1";
		tmp.append("p");
	    } else {
		//cerr << "0";
		tmp.append("m");
	    }	    
	}
	//cerr << tmp <<endl;
	extension_array.push_back(tmp);
    }
    
    for(long i=0;i<no_files;i++){
	string extension=extension_array[i];
	string synteny_file_name=synteny_file_name_prefix+extension;
	
	FILE* fptr=FOPEN((char*)synteny_file_name.c_str(),"r");
	if(fptr==NULL){
	    fptr=FOPEN((char*)synteny_file_name.c_str(),"w");
	    if(fptr==NULL){
		cerr << "Error: Unable to open output file"<< synteny_file_name << endl;
		exit(-1);
	    }
	    fprintf(fptr,"#-- MATCHES %s\n",extension.c_str());
	    for(long j=0;j<dimensions;j++){
		fprintf(fptr,"0\t");
	    }
	    fprintf(fptr,"\n");
	    for(long j=0;j<dimensions;j++){
		fprintf(fptr,"0\t");
	    }
	    fprintf(fptr,"\n\n");
	    fclose(fptr);
	}
	else{
	    fclose(fptr);
	}
	Interval tmp_interval;
	report_combinations(tmp_interval, synteny_file_name,extension,1);
    }
}

///////////////////////////Reporting Repeates

void Transformer::report_repeats(){
    
    // loop over choosen fragments in chain
    // mark them and report them in the correct file
    // extension
    //string command="rm -f "+workingDir+"repeats.*";
    string command="rm -f "+RepFile+"*";
    system(command.c_str());
    //cerr << command << endl;
    string extension="";
    
    //string repeat_file_name_prefix=workingDir+"repeats.";
    string repeat_file_name_prefix=RepFile;
    string synteny_file_name_prefix=repeat_file_name_prefix;
    FILE* fptr=NULL;
    for(unsigned long i=0;i<Intervals.size();i++){
	if(Intervals[i].choosen_flag==1)continue;
	// construct file extension
	for(long j=0;j<dimensions;j++){
	   
	    if(Intervals[i].isForward[j]==1){
		extension.append("p");
	    }
	    else{
		extension.append("m");
	    }	    
	}
	string synteny_file_name=synteny_file_name_prefix+extension;
	fptr=FOPEN((char*)synteny_file_name.c_str(),"r");
	if(fptr==NULL){
	    fptr=FOPEN((char*)synteny_file_name.c_str(),"w");
	    if(fptr==NULL){
		cerr << "Error: Unable to open output file"<< synteny_file_name << endl;
		exit(-1);
	    }
	    fprintf(fptr,"#-- MATCHES %s\n",extension.c_str());
	    for(long j=0;j<dimensions;j++){
		fprintf(fptr,"0\t");
	    }
	    fprintf(fptr,"\n");
	    for(long j=0;j<dimensions;j++){
		fprintf(fptr,"0\t");
	    }
	    fprintf(fptr,"\n\n");
	    fclose(fptr);
	}
	else{
	    fclose(fptr);
	}
	fptr=FOPEN((char*)synteny_file_name.c_str(),"a");
	if(fptr==NULL){
	    cerr << "Error: Unable to open output file"<< synteny_file_name << endl;
	    exit(-1);
	}
	// check origin
	int origin_flag=0;
	for(long j=0;j<dimensions;j++){
	    if(Intervals[i].b.coordinates[j]!=0){
		origin_flag=1;
		break;
	    }
	    if(Intervals[i].e.coordinates[j]!=0){
		origin_flag=1;
		break;
	    }
	    
	}
	if(origin_flag){
	    for(long j=0;j<dimensions;j++){
		if(Intervals[i].isForward[j]==1){
		    fprintf(fptr,"%lu ",max(0,Intervals[i].b.coordinates[j]));
		}
		else{
		    fprintf(fptr,"%lu ",max(0,Intervals[i].e.coordinates[j]));
		}
	    }
	    fprintf(fptr,"\n");
	    for(long j=0;j<dimensions;j++){
		if(Intervals[i].isForward[j]==1){
		    fprintf(fptr,"%lu ",max(0,Intervals[i].e.coordinates[j]));
		}
		else{
		    fprintf(fptr,"%lu ",max(0,Intervals[i].b.coordinates[j]));
		}
	    }
	    fprintf(fptr,"\n\n");
	}
	fclose(fptr);
	extension.clear();
	
    }
    //string repeat_file_name=repeat_file_name_prefix+extension;
//    command="cat "+workingDir+"repeats.* > "+workingDir+"repeats.dat";
    command="cat "+RepFile+"* > "+RepFile+"dat";
    
    system(command.c_str());
//    command="rm -f "+workingDir+"repeats.p*";
    command="rm -f "+RepFile+"p*";
    
    system(command.c_str());
    
}

void Transformer::generate_gp_files(string synteny_file_name_prefix){


    string command;
    for(unsigned long j=0;j<combinations.size();j++){

	//
	string tmp_file=synteny_file_name_prefix+combinations[j]+".gp";
	string ps_file=synteny_file_name_prefix+combinations[j]+".ps";
	string dat_file=synteny_file_name_prefix+combinations[j]+".dat";
	FILE* fptr=FOPEN(tmp_file.c_str(),"a");
	if(fptr==NULL){
	    cerr << "Unable to open gnuplot output file\n";
	    return;
	}
	//cerr << "HIIIIIIIIIIII\n";
	fprintf(fptr,"set terminal postscript landscape color\n");
	fprintf(fptr,"set out '%s'\n",ps_file.c_str());
	fprintf(fptr,"set nokey\n");
	fprintf(fptr,"set nogrid\n");
	fprintf(fptr,"set xlabel \"\"\n");
	fprintf(fptr,"set ylabel \"\"\n");
	//fprintf(fptr,"set xrange [0:*]\n");
	//fprintf(fptr,"set yrange [0:*]\n");
	int color=1;

	if(draft_flag){
	    //for(long i=0;i<dimensions;i++){
	    long poss=0;
	    ContigSorting* contigobj=(*contig_ptr_array)[proj1[j]];
	    fprintf(fptr,"set xrange [0:%lu]\n",
		    contigobj->recordseps[contigobj->number_of_contigs-1]);
	    fprintf(fptr,"set x2tics \(");
	    for(long kk=0;kk<contigobj->number_of_contigs;kk++){
		poss=contigobj->recordseps[kk];
		fprintf(fptr," \" \" %lu",poss);
		if(kk<contigobj->number_of_contigs-1){
		    fprintf(fptr,",");
		}
	    }
	    fprintf(fptr,")\n");
	    contigobj=(*contig_ptr_array)[proj2[j]];
	    fprintf(fptr,"set yrange [0:%lu]\n",
		    contigobj->recordseps[contigobj->number_of_contigs-1]);
	    fprintf(fptr,"set y2tics \(");
	    for(long kk=0;kk<contigobj->number_of_contigs;kk++){
		poss=contigobj->recordseps[kk];
		fprintf(fptr," \" \" %lu",poss);
		if(kk<contigobj->number_of_contigs-1){
		    fprintf(fptr,",");
		}
	    }
	    fprintf(fptr,")\n");
	    //long contig1=contigobj->
	    //   getrecordnum(contigobj->recordseps,contigobj->number_of_contigs,
	    //		 contigobj->total_conitg_length,
	    //		 Chain[Points[k-1]].b.coordinates[i]);
	    //}
	    fprintf(fptr,"set grid  noxtics noytics  x2tics y2tics\n");
	    
	}

	fprintf(fptr,"plot ");
	for(unsigned long i=0;i<extension_array.size();i++){
	    if((extension_array[i].c_str()[proj1[j]]=='p')&&(extension_array[i].c_str()[proj2[j]]=='p')){
		color=1;
	    }
	    else if((extension_array[i].c_str()[proj1[j]]=='p')&&(extension_array[i].c_str()[proj2[j]]=='m')){
		color=2;
	    }
	    else if((extension_array[i].c_str()[proj1[j]]=='m')&&(extension_array[i].c_str()[proj2[j]]=='p')){
		color=2;
	    }
	    else if((extension_array[i].c_str()[proj1[j]]=='m')&&(extension_array[i].c_str()[proj2[j]]=='m')){
		color=1;
	    }
	    
	    fprintf(fptr,"\"%s\" index %lu title \"%s\" w lp lt %d lw 2 pt 1 ps 0.75",dat_file.c_str(),i,extension_array[i].c_str(),color);
	    if(i<extension_array.size()-1){
		fprintf(fptr,", ");
	    }
	    else{
		fprintf(fptr,"\n");
	    }
	}
	fclose(fptr);
	//command="gnuplot "+tmp_file;
	//system(command.c_str());
    }

}



////////////////////////// For compact chain


void Transformer::plot_compact_chain(){
    

        
    for(unsigned long i=0;i<Intervals.size();i++){
	Intervals[i].choosen_flag=0;
    }
    for(unsigned long i=0;i<Chain.size();i++){
	Intervals[Chain[i].intervalID].choosen_flag=1;
    }

    string command="rm -f "+SynFilePrefix+"*";
    //cerr << "HALOOOOOOOOOOOOOOO: " << command << endl;
    system(command.c_str());
    string extension="";
    //string synteny_file_name_prefix=workingDir+"syntenies.";
    string synteny_file_name_prefix=SynFilePrefix;

    generate_combinations();
    create_files_for_all_suffixes(synteny_file_name_prefix);
    
    FILE* fptr=NULL;
    for(unsigned long i=0;i<compact_chain.size();i++){
	// construct file extension
	for(long j=0;j<dimensions;j++){
	    //cout << "["<<Chain[index].b.coordinates[i]<< ","<< Chain[index].e.coordinates[i] << "] ";
	    if(compact_chain[i].isForward[j]==1){
		extension.append("p");
	    }
	    else{
		extension.append("m");
	    }	    
	}
	string synteny_file_name=synteny_file_name_prefix+extension;
	fptr=FOPEN((char*)synteny_file_name.c_str(),"r");
	if(fptr==NULL){
	    fptr=FOPEN((char*)synteny_file_name.c_str(),"w");
	    if(fptr==NULL){
		cerr << "Error: Unable to open output file"<< synteny_file_name << endl;
		exit(-1);
	    }
	    fprintf(fptr,"#-- MATCHES %s\n",extension.c_str());
	    for(long j=0;j<dimensions;j++){
		fprintf(fptr,"0\t");
	    }
	    fprintf(fptr,"\n");
	    for(long j=0;j<dimensions;j++){
		fprintf(fptr,"0\t");
	    }
	    fprintf(fptr,"\n\n");
	    fclose(fptr);
	}
	else{
	    fclose(fptr);
	}
	fptr=FOPEN((char*)synteny_file_name.c_str(),"a");
	if(fptr==NULL){
	    cerr << "Error: Unable to open output file"<< synteny_file_name << endl;
	    exit(-1);
	}
	// check origin
	int origin_flag=0;
	for(long j=0;j<dimensions;j++){
	    if(compact_chain[i].b.coordinates[j]!=0){
		origin_flag=1;
		break;
	    }
	    if(compact_chain[i].e.coordinates[j]!=0){
		origin_flag=1;
		break;
	    }
	    
	}
	if(origin_flag){
	    report_combinations(compact_chain[i], synteny_file_name,extension,0);
	    for(long j=0;j<dimensions;j++){
		if(compact_chain[i].isForward[j]==1){
		    if(compact_chain[i].b.coordinates[j]<shiftvalue)
			fprintf(fptr,"%lu ",compact_chain[i].b.coordinates[j]);
		    else fprintf(fptr,"%lu ",compact_chain[i].b.coordinates[j]-shiftvalue);
		}
		else{
		    if(compact_chain[i].e.coordinates[j]<shiftvalue)
			fprintf(fptr,"%lu ",compact_chain[i].e.coordinates[j]);
		    else fprintf(fptr,"%lu ",compact_chain[i].e.coordinates[j]-shiftvalue);
		}
	    }
	    fprintf(fptr,"\n");
	    for(long j=0;j<dimensions;j++){
		if(compact_chain[i].isForward[j]==1){
		    if(compact_chain[i].e.coordinates[j]<shiftvalue)
			fprintf(fptr,"%lu ",compact_chain[i].e.coordinates[j]);
		    else
			fprintf(fptr,"%lu ",compact_chain[i].e.coordinates[j]-shiftvalue);
		    
		}
		else{
		    if(compact_chain[i].b.coordinates[j]<shiftvalue)
			fprintf(fptr,"%lu ",compact_chain[i].b.coordinates[j]);
		    else
			fprintf(fptr,"%lu ",compact_chain[i].b.coordinates[j]-shiftvalue);
		}
	    }
	    fprintf(fptr,"\n\n");
	}
	fclose(fptr);
	extension.clear();
	
    }
/*
    for(long i=0;i<Intervals.size();i++){
	Intervals[i].choosen_flag=0;
    }
    for(long i=0;i<compact_chain.size();i++){
	Intervals[compact_chain[i].intervalID].choosen_flag=1;
    }
*/
    //string repeat_file_name=repeat_file_name_prefix+extension;
    
    // for plotting
    command.clear();
    for(unsigned long j=0;j<combinations.size();j++){
	for(unsigned long i=0;i<extension_array.size();i++){
	    command.clear();
	    string tmp_file=synteny_file_name_prefix+extension_array[i]+"."+combinations[j];
	    FILE* tmp_fptr=FOPEN(tmp_file.c_str(),"a");
	    fprintf(tmp_fptr,"\n");
	    fclose(tmp_fptr);
	    command="cat "+synteny_file_name_prefix+extension_array[i]+"."+combinations[j]+" >> "
		+synteny_file_name_prefix+combinations[j]+".dat";
	    //cerr << command<<endl;
	    system(command.c_str());
	    string x=synteny_file_name_prefix+extension_array[i]+"."+combinations[j];
	    unlink(x.c_str());
	}
	
    }

    
//    command="cat "+workingDir+"syntenies.* > ";
//    command=command+workingDir+"syntenies.dat";

    command="cat "+SynFilePrefix+"p* > ";
    command=command+SynFilePrefix+"dat";
    
    system(command.c_str());
//    command="rm -f "+workingDir+"syntenies.p*";
    command="rm -f "+SynFilePrefix+"p*";

    system(command.c_str());
    generate_gp_files(synteny_file_name_prefix);


}



void Transformer::report_permutations_from_compact_chain(){    

    cout << "\n# Compact Permutations w.r.t. identity permutation: "<< endl<< endl;
 


   // store k-tuple of ids;
    vector <Point> ktuple;
    vector <Point> ktuple_dir;
    vector<long> csb(dimensions);
    vector<long> dirp(dimensions);
    vector<long> inverse_index(dimensions);
    vector <Point> ktuple_inverse;
    vector<long> sorted_wrt_x;
    //Point id_s;

    int set_identity_flag=0;
    Chain.clear();
    //cout << "Compact chain size "<< compact_chain.size() << endl;
    for(unsigned long i=0;i<compact_chain.size();i++){
	Chain.push_back(compact_chain[i]);
    }

//    cout << "Chain size "<< Chain.size() << endl;

    Points.clear();
    Types.clear();
    for(unsigned long j=0;j<Chain.size();j=j+1){	
	Points.push_back(j);
	Types.push_back(1);	
    }

    for(long i=0;i<dimensions;i++){
       // sort w.r.t. genome and report ids, we will report its index later
       // sort points w.r.t. 
       quick_sort_list(i+1);
       if(set_identity_flag==0){
	   set_identity_flag=1;
	   long counter=0;
	   for(unsigned long x=0;x < Points.size();x++){
	       if(Types[x]){
		   Chain[Points[x]].interval_identity=counter+1;
		   if(i==0){
		       //  sorted_wrt_x.push_back(Points[x]);
		   }
		   counter++;
		   
		   csb[0]=Chain[Points[x]].interval_identity;
		   Point id_s(csb);	
		   ktuple.push_back(id_s);

		   dirp[0]=1;	   		   
		   Point dir_s(dirp);
		   ktuple_dir.push_back(id_s);	
		   
		   inverse_index[0]=counter+1;
		   ktuple_inverse.push_back(inverse_index);
		   
	       }
	   }
       }
       cout << "> Genome " << i+1 << "\n";
       for(unsigned long k=0;k<Points.size();k++){
	   if(Types[k]==1){ // it is a start point
	       
	       // report chromosome boundary
	       if((draft_flag)&&(k>0)){
		   ContigSorting* contigobj=(*contig_ptr_array)[i];
		   long contig1=contigobj->
		       getrecordnum(contigobj->recordseps,contigobj->number_of_contigs,
				    contigobj->total_conitg_length,
				    Chain[Points[k-1]].b.coordinates[i]);
		   long contig2=contigobj->
		       getrecordnum(contigobj->recordseps,contigobj->number_of_contigs,
				    contigobj->total_conitg_length,
				    Chain[Points[k]].b.coordinates[i]);
		   if(contig1!=contig2){
		       //cout << ">"<< contig1 <<","  << contig2<< "< ";
		       cout << " | ";
		   }
	       }
	       /// end report chromosome boundary

	       if(Chain[Points[k]].isForward[i]==1){
		   cout << "+";
		   ktuple_dir[k].coordinates[i]=1;
	       }
	       else{
		   cout << "-";
		   ktuple_dir[k].coordinates[i]=0;
	       }
	       cout << Chain[Points[k]].interval_identity<<" ";	       
	       if(1){
		   ktuple[k].coordinates[i]=Chain[Points[k]].interval_identity;
		   ktuple_inverse[Chain[Points[k]].interval_identity-1].coordinates[i]=k;
	       }
	   }
	   else if(Types[k]==0){
	       
	       //cout << "Coord: " << Chain[Points[i]].e.coordinates[dimension-1]<< " "<<Types[i]<< " "
	       //	 << Points[i]<<" score: "<< Chain[Points[i]].s <<endl;	       	       
	   }
       }
       cout << endl;
   }

    cout << endl;
}
