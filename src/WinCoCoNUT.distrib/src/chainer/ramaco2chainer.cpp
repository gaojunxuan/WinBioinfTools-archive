#include <stdio.h>
//#include <string.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <vector>

using namespace std;

#define max(a, b)  ((a)>=(b) ? (a) : (b))

#define FAIL 1
#define SUCCESS 0

void fsplit (string line, const char* sep, vector<string> & words);

int main(int argc,char *argv[]){

  if(argc!=3)
    {
	cout << "Invalid options\n";
	cout << "Usage: ramaco2chainer filename orientation\n";
	cout << "orientation: strand orientations except for first seq in +|- \n";    
	cout << "Example for 3 sequences with the -ve strand of the 3d:\n";
	cout << "> ramaco2chainer Seqfor3genomes ++ \n";
	cout << "Output:\n ";
	cout << "-One file: 'filename.pm' in case of multimat input file \n";    
    return(FAIL);
  } 
  
  
  ifstream inputf(argv[1], ios::in);
  if (!inputf) {
      cerr << "Error in opening input file!\n";
      exit (1);
  }
  

  // create output file name
  string in_suffix="";
  in_suffix.append((char*)argv[2]);
  long number_of_genomes=in_suffix.size()+1;
  cout << "Input file: " << argv[1] << "\n";
  cout << "Number of genomes: " << number_of_genomes << endl;
    
  vector <long> query_size_array;
  string suffix="p";
  for(unsigned long i=0;i<in_suffix.size();i++){
      if(in_suffix.c_str()[i]=='+'){
	  suffix.append("p");
      }
      else{
	  suffix.append("m");
      }
  }
  string outfilename=argv[1];
  outfilename.erase(outfilename.size()-number_of_genomes+1);
  outfilename.append(".");
  outfilename.append(suffix);

  ofstream outputf(outfilename.c_str(), ios::out);
  if (!outputf) {
      cerr << "Error in opening file" << outfilename << "!\n";
      exit (1);
  }
  outputf << ">CHA " << number_of_genomes << "\n";

  //cout << "In Suffix: +" << in_suffix << "Stored Suf: "<< suffix << endl;
  string line;
  vector< string> words; 
  long length;
  vector<long> positions; 

  while(getline (inputf, line)){
      string direction="p";
      if(line.substr(0,1).compare("#")==0){	  
	  words.clear(); 
	  fsplit(line, " ", words);
	  if(words.size()>3){
	      //cout << line << endl;
	      if(words[1].substr(0,9).compare("queryfile")==0){
		  //cout <<  "Query Seq. Size: " << words[words.size()-1] << endl;
		  query_size_array.push_back(atol(words[words.size()-1].c_str()));
		  cout <<  "Query Seq. Size: " << query_size_array[query_size_array.size()-1] << endl;
	      }
	  }
      }
      else{ // handling matches
	  words.clear(); 
	  fsplit(line, " ", words);words[0].c_str();
	  length=atol(words[0].c_str());
	  outputf << "#" << length << endl;
	  for(unsigned long i=1;i<words.size();i++){
	      if(words[i].c_str()[0]!='-'){ // not reverse complement
		  //positions.push_back(atol(words[i].c_str()));
		  long tmpnum=atol(words[i].c_str());
		  outputf << "[" << tmpnum << "," << tmpnum+length-1 << "] ";
		  //direction.append("p");
	      }
	      else{		  
		  long tmp_pos1=atol(words[i].substr(1,words[i].size()-1).c_str());
		  long tmp_pos2=tmp_pos1+length-1;
		  tmp_pos1=query_size_array[i-2]-tmp_pos1-1;
		  tmp_pos2=query_size_array[i-2]-tmp_pos2-1;
		  outputf << "[" << tmp_pos2 << "," << tmp_pos1 << "] ";
		  //positions.push_back(tmp_pos);
		  //direction.append("m");
	      }
	  }
	  outputf << "\n";
	  // printing to file
	  for(unsigned long i=0;i<positions.size();i++){
	      //  outputf;
	  }
      }
  }

  inputf.close();
  outputf.close();

  return(0);
}


void fsplit (string line, const char* sep, vector<string> & words) {
		
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
