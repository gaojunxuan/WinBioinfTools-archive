#include "Point.h"
using namespace std;

Point::Point()
{

}

Point::Point(vector<long int> cs, Interval* myI, bool iSP){
  coordinates = cs;
  myInterval = myI;
  isStartPoint = iSP;
} 

Point::Point(vector<long int> cs){
  coordinates = cs;
}

/*
Point& Point::operator= (const Point& p){
  if (this != &p){
    coordinates.resize(p.coordinates.size());
    for (int i = 0; i < p.coordinates.size(); i++){
      coordinates[i] = p.coordinates[i];
    }
  }
  return *this;
}
*/

bool Point::smallerP (const Point *p, const int d) const{
  return coordinates[d] < p->coordinates[d];
}
