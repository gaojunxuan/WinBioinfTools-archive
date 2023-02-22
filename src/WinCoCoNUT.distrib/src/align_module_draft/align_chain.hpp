#ifndef cdna_h_INCLUDED
#define cdna_h_INCLUDED

#include <string>
#include "fragment.hpp"
#include "memmap.hpp"

#include "contigs.hpp"

using namespace std;


class Multialg
{
//public:
//    Multialg();
//    ~Multialg();
public:
    // ArrayPosition startpos;
//  ArrayGaptype gaps;
    long * startpos;
    char*  gaps;
  long alignedlength;
};



class align_chain {
public:
    /* Vektor der einzelnen Fragmente */
    vector<fragment *> fragmente;
    long number_of_genomes;

    /* ID von vmatch */
    long id;
    
    /* Bezeichnung der cDNA */
    string name;

    /* Laenge der cDNA Sequenz */
    int c_laenge;
    
    //int fasta_line_size;
    
    /* gDNA (Sub-)Sequenz */
    string gdna_subseq;
    
    string* genomic_sequences;
    interval* chain_boundaris;

    char*   clustalw_name;

    /* cDNA Sequenz */
    string cdna_seq;

    // genome file handles
    int* genome_files_handles;
    long* start_offset_array;

    long* genome_size_array;
    string direction;

    /* Teilsequenzen der (nicht gefundenen) Raender */
    string cdna_vorn, gdna_vorn, cdna_hinten, gdna_hinten;
    
    /* absoluter gDNA Start- und Endpunkt der SubSequenz */
    int g_start;
    int g_ende;

    // threshold under which gap is aligned by alignment programn (clustalw)
    long threshold_alignable_gap;

    // thrshold to decide if to read the region of the chain or not, defulat is 5 times page size
    long boundary_threshold;


    // output formatting
    int output_mode; // lineal or fasta
    int fragment_output_mode;  // sequence in addition or not 
    int seq_no_digit; // for formating digit output
    long* contig_offset_array; // for reporting relative position
    long* contig_id_array; // for reporting contig_id
    int palindrome_flag; // for reporting coordinae not w.r.t. reverse complement
    int relative_flag; // to descriminate between reporting draft or not, important for palindrome report

    // Statistics
    long unaligned_gap_counter;
    long* total_fragment_lengths_in_genomes;
    long* total_aligned_gap_lengths_in_genomes;
    long* total_unaligned_gap_lengths_in_genomes;
    long* fasta_line_size_array;
    long total_identity;

    contigs* all_contigs_array;

    /* Parst die eingelesene Zeile und splittet sie in Id, c_tsart, c_ende,
       g_start und g_ende auf
       line   [in]   string   aufzuteilende Zeichenkette
       sep [in]   char*	Menge von Trennzeichen
       words  [out]  vector<string> Vektor mit den getrennten Substrings */
    void fsplit (string line, const char* sep, vector<string> & words);
    

    
public:
    // Konstrukteur: wird gerufen, wenn Headerzeile gelesen wurde 
    align_chain();
    align_chain(long ,int* ,long*);
    
    // Destrukteur 
    ~align_chain ();
public:
    /* Fragment hinzufuegen: wird aufgerufen, wenn ein Fragmentzeile eingelesen
       wurde */
    void addFragment (string frag);
    
    

    int generatefilenames (char *fasta, char *gde, char *dnd);

    /* set chain id */    
    void setID(long in_id){id=in_id;};
    void setClustalw(char* in_clustalw_name){clustalw_name=in_clustalw_name;};
        
    void set_boundary_threshold(long in_thresh){boundary_threshold=in_thresh;};
    void set_fasta_line_size_array(long* in_fasta_line_size_array){fasta_line_size_array=in_fasta_line_size_array;};
    void set_contig_offset_array(int in_relative_flag,long* in_contig_offset_array,long* in_contig_id_array,contigs* in_all_contigs_array,int in_palindrome_flag);
//      ermittelt die gesamte Sequenz der cdna oder die entsprechende gDNA-
//      Subsequenz. Die Sequenz wird dabei aus dem entsprechenden cDNA File
//      oder gDNA Chromosom File mit Memory Mapping gelesen */
//      string GetSeq (long start, long ende, long such_ID);
    
    void GetSeqFromFastaFile (string* subseq,long start, long ende, int file_handle,long,long);
    string GetSeqFromFastaFileReverseComplementOld (long start, long ende, long genome_size, int file_handle, long in_start_offset,long);
    void GetSeqFromFastaFileReverseComplement(string* in_gdna,long start, long ende, long genome_size,int file_handle, long in_start_offset,long);
    // Wenn alle Fragmente der cDNA eingelesen wurden, wird die gesamte cDNA-
    //   Sequenz aus dem File gelesen und in cdna_seq gespeichert und der
    //   entsprechende gDNA Bereich wird gelesen und in gdna_subseq gespeichert
    
    void setSeqs ();
    
//	int miniKorrektur2 (string s, string t, int k);
    
    /* berechnet die Luecken zwischen den Fragmenten und bearbeitet diese
       entsprechender der zutreffenden Faelle */
    
    void process_chain(int,int,long);
      
    int callclustalw (const char *fasta);

    void echosymbolstring2file (FILE *outf, char *w, long wlen);
 
    // handling clustwalw output
    int readclustalwalignment (fragment*,fragment*,Multialg *alignment, char *gde, long numofsequences,long accumulated_align_len);
    void readmultiplegdefileagain(unsigned char *newseq, unsigned char *input,long inputlen);
    
    unsigned char * getsequencebynumber (long number, long lenalg, unsigned char *all, long numofbytes);
    int output_formatted_alignment_single_line(fragment* in_frag1, fragment* in_frag2,unsigned char* allseq, long lenalg, long numofbytes, long accumulated_align_len);
    int output_formatted_alignment_fasta(fragment* in_frag1,fragment* in_frag2,unsigned char* allseq, long lenalg, long numofbytes, long accumulated_align_len);
    void report_fragment(fragment* in_frag,int file_handle, long accumulated_len);
    void report_unaligned_gap(fragment* in_frag1,fragment* in_frag2,int file_handle, long accumulated_len);
    
    void check_identity(long position, long in_lenalg,unsigned char* in_allseq,long in_numofbytes);
    void report_statistics();
    int check_chain_file_consistency();

    // functionalities for the reverse complement
    void set_reverse_complement_parameters(long* in_size_array,string in_dir){genome_size_array=in_size_array;direction=in_dir;};
    void GetSeqFromFastaFileGeneralized(string* subseq,long in_start, long in_end, long genome_id);
   
    int check_chain_boundary_against_file_size (long start, long ende, int file_handle, long in_start_offset, long fasta_line_size); 
};

#endif
