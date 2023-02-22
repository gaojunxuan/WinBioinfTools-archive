#include "Transformer.h"
#include <iostream>
#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include <fstream>

#include "ContigSorting.h"

#include <sys/types.h>
#include <sys/stat.h>

//#include <sys/dir.h>
#include <sys/param.h>


using namespace std;




//start testing boost
//#include <boost/test/unit_test.hpp>
//using namespace boost::unit_test;

/*
test_suite*
init_unit_test_suite( int argc, char* argv[] )
{
    test_suite* test= BOOST_TEST_SUITE( "const_string test" );

    return test;
}
*/

//end test



int main(int argc, char* argv[]){


    long number_of_genomes; 
    char* infilename;

	if (argc < 3){
		cerr<<"usage: \n"<<endl;
		cerr<<"> chainer2permutation.x infilename.dat number_of_genomes [options]\n" <<endl;
		cerr<<"Options: \n"<<endl;
		cerr<<"-overlap1 x: allow percentage overlap, this requires overlap ratio 0<x<1\n";
		cerr<<"-overlap2:   allow overlap as half the minimum block length\n";
		cerr<<"-filterrep x:  filter repeats overlapping with percentage larger than 0<x<1\n";
		cerr<< "              Default:1D chain without overlapping are computed and no repeats will be filtered\n";
		cerr<<"-draft seq1.seqinfo ...seqk.seqinfo: process  multichromosomal or draft genomes.\n";
		cerr<<"           Info file about chromomsome (contig) lengths is needed for each genome\n";
		
		//cerr<< "The highr the overlap ratio, the higher the allowed overlap in chain\n";
		cerr<< "\nExample1:\n";
		cerr<< "> chainer2permutation.x infile.dat 3 -overlap1 0.5 -filterrep 0.8\n";
		cerr<< "\nExample for two multi-chrom. genomes:\n";
		cerr<< "> chainer2permutation.x infile.dat 2 -overlap1 0.5 -filterrep 0.8 -draft seq1.info seq2.info\n";
		return -1;
	}
	infilename = argv[1];
	string full=infilename;
	
	ifstream input (infilename,ios::in);
	if (!input) {
	    cerr << "Error in opening: " << infilename << " \n";
		exit (1);
	}
	input.close();

	string working_dir;
	//std::string::size_type idx = full.find_last_of("/");
	cout << "# Input file: "<< infilename << endl;
	unsigned long idx=full.find_last_of("/");
	if (idx == std::string::npos){
	    //cout << "# Input file: "<< full.substr(idx+1) << endl;
	    working_dir="./";
	    //cout << "# Working dir: " << working_dir<< endl; ;	    	    
	}
	else{
	    //cout << full.substr(idx+1) << endl;
	    working_dir= full.substr(0,idx+1);
	    //cout << "# Working dir: " <<  working_dir<<endl;
	}

	int overlap_function=1;
	float overlap_ratio=0;
	number_of_genomes = atoi(argv[2]);
	int overlap_flag=0;
	int filterrep_flag=0;
	float filter_rep_ratio=1;
	int draft_flag=0;
	vector <string> seqinfo_array;
	if(argc >3){
	    for(int i=3;i<argc;i++){
		if(strcmp(argv[i],"-overlap1")==0){
		    overlap_flag=1;
		    if(argc >i+1){
			overlap_ratio=atof(argv[i+1]);
			if((overlap_ratio>=1)||(overlap_ratio<0)){
			    cerr << "Error: overlap ratio must be greater than 0 and less than 1\n";
			    exit(-1);				
			}
			i++;
			overlap_function=1;	
		    }
		    else{
			cerr << "Error: No argument for -overlao option\n";
			return(-1);
		    }
		}
		else if(strcmp(argv[i],"-overlap2")==0){
		    overlap_flag=1;
		    overlap_function=0;		    
		}
		else if(strcmp(argv[i],"-draft")==0){
		    draft_flag=1;
		    if(i+number_of_genomes>=argc){
			cerr << "Error: " << number_of_genomes <<" of seqinfo files are required\n";
			return(-1);
		    }
		    for(int k=0;k<number_of_genomes;k++){
			seqinfo_array.push_back(argv[i+1]);
			i++;
		    }
		}
		else if(strcmp(argv[i],"-filterrep")==0){
		    //overlap_flag=1;
		    filterrep_flag=1;
		    if(argc >i+1){			
			filter_rep_ratio=atof(argv[i+1]);
			if((filter_rep_ratio>=1)||(filter_rep_ratio<0)){
			    cerr << "Error: overlap ratio must be greater than 0 and less than 1\n";
			    exit(-1);			    
			}
			i++;
		    }
		    else{
			cerr << "Error: No argument for -overlao option\n";
			exit(-1);
		    }
		}	    	
		else{
		    cerr << "Error: unknown option: "<< argv[i] << endl;
		    return(-1);
		}
	    }
	}

	vector <ContigSorting*> contig_ptr_array;
	if(draft_flag==1){	    	    
	    for(int k=0;k<number_of_genomes;k++){
		ifstream input ((char*)seqinfo_array[k].c_str(),ios::in);
		if (!input) {
		    cerr << "Error in opening: " << seqinfo_array[k] << " \n";
		    exit (1);
		}
		input.close();
	    }
	    input.close();
	    for(int k=0;k<number_of_genomes;k++){		
		ContigSorting* Contigptr= new ContigSorting;	    
		Contigptr->read_contigs((char*)seqinfo_array[k].c_str(),0);
		contig_ptr_array.push_back(Contigptr);
		/*
		ContigSorting* contigobj=contig_ptr_array[contig_ptr_array.size()-1];
		long contig2=contigobj->
		    getrecordnum(contigobj->recordseps,contigobj->number_of_contigs,
				 contigobj->total_conitg_length,
				 857550);
		cerr << "HALOOOOOOOOOOOOOOOOO : "<< contig2 << endl;
		*/
	    }
	}


	//TransformerWithoutOverlap myTransformer(infilename, number_of_genomes);
	//myTransformer.transformIt();

	PlainTransformer myPlainTransformer(infilename,number_of_genomes);
	myPlainTransformer.workingDir=working_dir;
	if(overlap_flag==1){
	    myPlainTransformer.set_overlap_flag(overlap_flag);
	    myPlainTransformer.set_filterrep_flag(filterrep_flag,filter_rep_ratio);
	    myPlainTransformer.set_overlap_function_type(overlap_function);
	    myPlainTransformer.set_overlap_ratio(overlap_ratio);
	    if(draft_flag==1){
		myPlainTransformer.set_draft_flag(draft_flag,&contig_ptr_array);	    
	    }
	}
	
	myPlainTransformer.transformIt();
	
	

	if(draft_flag==1){	    
	    for(int k=0;k<number_of_genomes;k++){	
		delete contig_ptr_array[k];	    
	    }
	    contig_ptr_array.clear();
	}

/*
	char pathname[MAXPATHLEN];
	struct stat fstatus;
	stat(infilename,&fstatus);
*/

    

  return 0;
}
