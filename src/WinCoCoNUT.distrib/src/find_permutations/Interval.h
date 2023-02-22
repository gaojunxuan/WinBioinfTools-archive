#include <vector>
#include <iostream>

#include "Point.h"

using namespace std;

#ifndef _INTERVAL_
#define _INTERVAL_

struct interval_one_dim{
  long int p1; 
  long int p2;
};

class Interval{
 public:
 public:
  Point b; //begin
  Point e; //end
  float w; //weight
  float s; //score
  //bool isForward;
  vector<bool> isForward;
  // Interval* precedingInterval;
  char choosen_flag;
  long index_perc;
 public:
  long intervalID;
  long interval_identity; // for reporting final permutations w.r.t. identity

 public:
  Interval();
  ~Interval();
  Interval(Point begp, Point endp);
  //Interval(Point begp, Point endp, bool iF);
  Interval(Point begp, Point endp, vector<bool> iF);
  int getDirection(int dim);
	bool isDimForward(int dim);
  float weight();
  float weight(unsigned int dimension);
  Point* begin();
  Point* end();
  float score();
  Interval* prec();
  int setPrecedingInterval(Interval*);
  int setScore(float score);
  int setCoordinates(Point begp, Point endp);
  int setWeight(float weight);
  int print();
  int print(unsigned int dimension);
  int getID();
  void setID(long id);
  //copy operator (=) implementieren!
  //bool operator< (const Interval i2);
  interval_one_dim getOneDimension(int dim);
  Interval* findOverlaps(Interval* i2);
  bool overlap(Interval*i2);

};

bool operator< (Interval leftI, Interval rightI);
bool operator< (interval_one_dim leftI, interval_one_dim rightI);
//bool overlap(Interval* i1, Interval* i2);

#endif
