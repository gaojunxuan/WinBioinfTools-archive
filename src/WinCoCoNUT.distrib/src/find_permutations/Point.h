 
#include <vector>

using namespace std;

#ifndef _POINT_
#define _POINT_

class Interval;

class Point{
 private:
 protected:
 public:
  
  vector<long> coordinates;
  Interval* myInterval;
  int myIntervalIndex;
  bool isStartPoint;
  Point();
  Point(vector<long> cs, Interval* myI, bool iSP);
  Point(vector<long> cs);
  //Point& operator= (const Point&);
  //bool operator<();
  bool smallerP (const Point *p, const int d) const;
	
};

#endif
