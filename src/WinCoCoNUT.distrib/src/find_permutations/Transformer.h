#include "Interval.h"
#include "kdTree.h"

#include <map>
#include "ContigSorting.h"

//#include <boost/regex.hpp> //gerade herauskommmentiert
//using namespace boost::regex;

class Transformer{
	private:
	
	public:
		char* inDatFile;
		std::string workingDir;
		unsigned int dimensions;
		string SynFilePrefix;
		string RepFile;
		vector<class Interval> Intervals;
		vector<class Interval> Chain;
		vector<class Interval> NewChain;
		vector <long> Points; // stores indices for intervals
		vector <long> Types; // stores type of point
		vector <Interval> compact_chain;

	public:
		Transformer(char* datFile, unsigned int dim);
		virtual ~Transformer();
		virtual int transformIt() = 0;
		

		/// my functions
		float overlap_ratio;
		Kdelem_ptr *perm;    //stores the elements of the tree in order
		Kdelem_ptr nntarget, first, tail; 
		Kdnode_ptr *bucket_ptr;  //stores the bucket location of each element
		Kdelem_ptr blocks;   //permanent storage of the elements
		
		int parseDatFilekGenomes(char* infile);
		void fsplit (string line, const char* sep, vector<string> & words);
		int quick_sort_list(long in_dimension);
		void recursive_subdim_sort(long in_dimension, long startkeys, long endkeys);
		void quicksortkdlistkeysOndim(long dim,long lo, long hi);
		long get_dimension_value(long dim, long index, long type_index);

		long chainOneDirection(int dimension);
		int chooseFragments();
		void add_origin_to_chain();
		void report_permutations();
		long chainOneDirectionWithOverlap(long);
		void free_blocks_array();
		long get_position_v_point(long index,long dimension);
		vector<long> min_len;
		int overlap_function_type;
		void set_overlap_function_type(int in){overlap_function_type=in;};
		void set_overlap_ratio(float in){overlap_ratio=in;};
		int overlap_flag;
		void set_overlap_flag(int in){overlap_flag=in;};
		void report_chain();
		void plot_compact_chain();
		void report_repeats();
		void generate_combinations();
		void report_combinations(Interval input_interval,string,string,int);
		void create_files_for_all_suffixes(string synteny_file_name_prefix);
		vector<string> combinations;
		vector <string> extension_array;
		vector<long> proj1;
		vector<long> proj2;
		void generate_gp_files(string);
		long filter_repeat(long);
		int filter_repeats_all_dimensions();
		void report_compact_chain(vector <Point>*);
		int shiftvalue;
		long Total_final_score;
		int filterrep_flag;
		float filter_rep_ratio;
		void set_filterrep_flag(int in_flag, float in_ratio){filterrep_flag=in_flag;
		filter_rep_ratio=in_ratio;};
		int draft_flag;
		vector <ContigSorting*>* contig_ptr_array;
		void set_draft_flag(int in_draft_flag,vector <ContigSorting*>* in_contig_ptr_array ){draft_flag=in_draft_flag;
		contig_ptr_array=in_contig_ptr_array;};

		void report_permutations_from_compact_chain();

};




class PlainTransformer : public Transformer{
	protected:
		int clearData();

	public:
		PlainTransformer(char* datFile, unsigned int dim);
		virtual ~PlainTransformer();
		int transformIt();
};

class TransformerWithoutOverlap : public Transformer{
	public:

	protected:
		int clearData();
		
		
		int getOverlapsAndChange(float percentage);
		int findOverlapsAndChange(unsigned idnc, unsigned idc, float percentage);
	public:
		TransformerWithoutOverlap(char* datFile, unsigned int dim);
		virtual ~TransformerWithoutOverlap();
		int transformIt();

};



class TransformerWithPraeprocessing: public Transformer{
	protected:
		int clearOverlaps(float percentage);
		
	public:
		TransformerWithPraeprocessing(char* datFile, unsigned int dim);
		virtual ~TransformerWithPraeprocessing();
		int transformIt();
	
};
