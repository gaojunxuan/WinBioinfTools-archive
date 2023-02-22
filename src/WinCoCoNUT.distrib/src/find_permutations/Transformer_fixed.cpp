#include <iostream>
using std::cout;

#include<iterator>
using std::ostream_iterator;

#include "Transformer.h"

Transformer::Transformer(char* datFile, unsigned int dim):inDatFile(datFile), dimensions(dim){

	Intervals = new vector<class Interval>;
	permutations = new vector < vector<int> >;
	xIntervalsIDs = new set<int>;
	yIntervalsIDs = new set<int>;
	zIntervalsIDs = new set<int>;
	sortedXIntervalsIDs = new vector<int>;
	sortedYIntervalsIDs = new vector<int>;
	sortedZIntervalsIDs = new vector<int>;
	choosenIntervalsIDs = new set<int, less<int> >;
	permutationY = new vector<int>;
	if (dimensions > 2)
		permutationZ = new vector<int>;
	notInPermutationY = new vector<int>;
	if (dimensions > 2 )
		notInPermutationZ = new vector<int>;
	InputID2ID = new map<int, int>;
	allNotChoosenIntervalsIDs = new set<int, less<int> >;
}

Transformer::~Transformer(){
	delete Intervals;
	delete permutations;

	delete xIntervalsIDs;
	delete yIntervalsIDs;
	delete zIntervalsIDs;

	delete sortedXIntervalsIDs;
	delete sortedYIntervalsIDs;
	delete sortedZIntervalsIDs;
	delete choosenIntervalsIDs;
	delete permutationY;
	if (dimensions > 2 ) delete permutationZ;
	delete notInPermutationY;
	if (dimensions > 2 ) delete notInPermutationZ;
	delete InputID2ID;
	delete allNotChoosenIntervalsIDs;
} 

/*int Transformer::parseDatFile(char* infile){
//could be about two or more genomes / dimensions
	ifstream infilestream(infile);
	if (!infilestream){
		cout<<"cannot open infile "<<infile<<endl;
		exit(0);
	}
	char* fileline = new char[128];
	vector<long> firstpoints(dimensions);
	vector<long> secondpoints(dimensions);
	vector<bool> isnew(dimensions);
	vector<bool> forwardMatches(dimensions);
	string currentline;
	
	while(!infilestream.eof()){
		infilestream >> fileline;
		if (strcmp(testcomment, "") == 0 || strcmp(testcomment, "\n") == 0)  continue;

	}
  //parse infile (datfile)
}*/

int Transformer::parseDatFile2Genomes(char* infile){
ifstream infilestream(infile);
  if (!infilestream){
    cerr<<"cannot open infile "<<infile<<endl;
    exit(0);
  }
  //parse infile (datfile)
  char* testcomment = new char[128];
  long x1,y1,x2,y2= 0;
  long ox1 = 0; long oy1 = 0; long ox2 = 0; long oy2 = 0;
  bool isnewx, isnewy;
  vector<bool> forwardMatches(2); //direction for each dimension/genome 
  string currentline;  

  forwardMatches[0] = 0;  
  while(!infilestream.eof())
  {
    if (infilestream.peek() == '\n'){
      infilestream.ignore(999, '\n');
      continue;
    }
    if (infilestream.peek() == '0'){
      infilestream.ignore(999, '\n');
      continue;
    }
    if (infilestream.peek() == '#'){
      if (forwardMatches[0] == 0) {
    	forwardMatches[0] = 1; forwardMatches[1] = 1; }
      else
        forwardMatches[1] = 0;  
      infilestream.ignore(999, '\n');  
      continue;
    }  
    infilestream >> x1;
    infilestream >> y1;
    infilestream >> x2;
    infilestream >> y2;

    if (infilestream.eof()) break;
    
//    cerr << x1 << "  " << y1 << "  " << x2 << "  " << y2 << endl;

    if (x1 == ox1 && x2 == ox2) isnewx = 0;
    else isnewx = 1;
    if(y1 == oy1 && y2 == oy2) isnewy = 0;
    else isnewy = 1;
    
    if (isnewx || isnewy){
      vector<long int> csb(2); vector<long int> cse(2);
      csb[0] = x1; csb[1] = y1; cse[0] = x2; cse[1] = y2; 

      /*cout<<"Startpoint = ("<<x1<<", "<<y1<<", "<<z1<<")\n";
      cout<<"Endpoint = ("<<x2<<", "<<y2<<", "<<z2<<")\n";
      cout<<"forwardMatches = "<<forwardMatches[0]<<", "<<forwardMatches[1]<<", "<<forwardMatches[2]<<endl;*/

      Point b(csb); Point e(cse);
      Interval iv(b,e,forwardMatches);
      Intervals->push_back(iv);
      ox1 = x1; oy1 = y1; ox2 = x2; oy2 = y2;
    }
  }
  infilestream.close();  
  cout<<"NUMBER of parsed intervals = "<<Intervals->size()<<endl;
  return 0;
}

int Transformer::parseDatFile3Genomes(char* infile){
  ifstream infilestream(infile);
  if (!infilestream){
    cout<<"cannot open infile "<<infile<<endl;
    exit(0);
  }
  //parse infile (datfile)
  char* testcomment = new char[128];
  long x1,y1,x2,y2,z1,z2 = 0;
  long ox1 = 0; long oy1 = 0; long ox2 = 0; long oy2 = 0; long oz1 = 0; long oz2 = 0;
  bool isnewx, isnewy, isnewz;
  vector<bool> forwardMatches(3); //direction for each dimension/genome 
  string currentline;  

  while(!infilestream.eof()){
    infilestream >>testcomment;
    if (strcmp(testcomment, "") == 0 || strcmp(testcomment, "\n") == 0)  continue;
    if (strcmp(testcomment,"#--")==0){
      infilestream>>testcomment; 
      infilestream>>testcomment;
      if (strcmp(testcomment,"ppp") == 0){
	forwardMatches[0] = 1; forwardMatches[1] = 1; forwardMatches[2] = 1;
	continue;
      }
      if (strcmp(testcomment,"ppm") == 0){
	forwardMatches[0] = 1; forwardMatches[1] = 1; forwardMatches[2] = 0;
	continue;
      }
      if (strcmp(testcomment,"pmp") == 0){
	forwardMatches[0] = 1; forwardMatches[1] = 0; forwardMatches[2] = 1;
	continue;
      }
      if (strcmp(testcomment,"pmm") == 0){
	forwardMatches[0] = 1; forwardMatches[1] = 0; forwardMatches[2] = 0;
	continue;
      }
    }
 
    x1 = atoi(testcomment);
    infilestream >> y1;
    infilestream >> z1;
    infilestream >> x2;
    infilestream >> y2;
    infilestream >> z2;

    if (x1 == ox1 && x2 == ox2) isnewx = 0;
    else isnewx = 1;
    if(y1 == oy1 && y2 == oy2) isnewy = 0;
    else isnewy = 1;
    if (z1 == oz1 && z2 == oz2) isnewz = 0;
    else isnewz = 1;
    
    if (isnewx || isnewy || isnewz){
      vector<long int> csb(3); vector<long int> cse(3);
      csb[0] = x1; csb[1] = y1; csb[2] = z1; cse[0] = x2; cse[1] = y2; cse[2] = z2;

      /*cout<<"Startpoint = ("<<x1<<", "<<y1<<", "<<z1<<")\n";
      cout<<"Endpoint = ("<<x2<<", "<<y2<<", "<<z2<<")\n";
      cout<<"forwardMatches = "<<forwardMatches[0]<<", "<<forwardMatches[1]<<", "<<forwardMatches[2]<<endl;*/

      Point b(csb); Point e(cse);
      Interval iv(b,e,forwardMatches);
      Intervals->push_back(iv);
      ox1 = x1; oy1 = y1; oz1 = z1; ox2 = x2; oy2 = y2; oz2 = z2;
    }
  }
  infilestream.close();
  //cout<<"NUMBER of parsed intervals = "<<Intervals->size()<<endl;
  return 0;
}

//gibt true zurueck, falls der Startpunkt im ersten Interval vor dem Startpunkt des zweiten Intervals in x-Richtung liegt
bool Transformer::startbefore(Interval i1, Interval i2){
	Point* startPoint1 = i1.begin();
	Point* startPoint2 = i2.begin();

	//just to make sure:
	Point* endPoint1 = i1.end();
	Point* endPoint2 = i2.end();
    
#if 0
    if ((*startPoint1).coordinates[0] >= (*endPoint1).coordinates[0])
        cerr << (*startPoint1).coordinates[0] << " >= " <<
            (*endPoint1).coordinates[0] << endl;
    if ((*startPoint2).coordinates[0] >= (*endPoint2).coordinates[0])   
        cerr << (*startPoint2).coordinates[0] << " >= " <<
            (*endPoint2).coordinates[0] << endl;
#endif
            
	assert((*startPoint1).coordinates[0] < (*endPoint1).coordinates[0]);
	assert((*startPoint2).coordinates[0] < (*endPoint2).coordinates[0]);
	return (*startPoint1).coordinates[0] < (*startPoint2).coordinates[0];
}

int Transformer::sortIntervalsXDirection(){
	sort(Intervals->begin(), Intervals->end(), &Transformer::startbefore);
  //ID in der Reihenfolge vergeben, in der Intervalle nach dem Sortieren bzgl. der X-Richtung vorliegen:
	for (unsigned int i = 0; i < Intervals->size(); i++){
		(*Intervals)[i].setID(i+1);
	}
	return 0;
}

int Transformer::getPermutation(){    //bei 3 Genomen mssen 2 Permutationen gefunden werden. Das erste Genom entspricht der Identit�spermutation, gefunden werden muss eine Permutation fr das zweite Genom bezglich des ersten Genoms und eine Permutation fr das dritte Genom bezglich des ersten Genoms

	int directionY;
	int directionZ;
	int intervalIDY;
	int intervalIDZ;

	permutationY->clear();
	notInPermutationY->clear();
	for (unsigned int i = 0; i < sortedYIntervalsIDs->size(); i++){
		intervalIDY = (*sortedYIntervalsIDs)[i];
		directionY = (*Intervals)[intervalIDY-1].getDirection(1);
		//checke, ob Interval auch in choosenIntervals liegt:
		if (choosenIntervalsIDs->find(intervalIDY) != choosenIntervalsIDs->end()){
			permutationY->push_back(directionY * (intervalIDY));
		}
		else{
			notInPermutationY->push_back(directionY * (intervalIDY));
		}
	}
	if(dimensions > 2){
		permutationZ->clear();
		notInPermutationZ->clear();
		for (unsigned int i = 0; i < sortedZIntervalsIDs->size(); i++){
			intervalIDZ = (*sortedZIntervalsIDs)[i];
			directionZ = (*Intervals)[intervalIDZ-1].getDirection(2);
	//checke, ob Interval auch in choosenIntervals liegt:
			if (choosenIntervalsIDs->find(intervalIDZ) != choosenIntervalsIDs->end()){
			permutationZ->push_back(directionZ * (intervalIDZ));
			}
			else{
			notInPermutationZ->push_back(directionZ * (intervalIDZ));
			}
		}
	}
  /*cout<<"Y-Permutation (size = "<<permutationY->size() <<"):\n";
  copy(permutationY->begin(),permutationY->end(), ostream_iterator<int>(cout," "));
  cout<<endl;
  
  cout<<"Z-Permutation (size = "<<permutationZ->size() <<"):\n";
  copy(permutationZ->begin(),permutationZ->end(), ostream_iterator<int>(cout," "));
  cout<<endl; */

  //dadurch dass Nummern rausfliegen bei der Auswahl, entstehen "Luecken". IDs so mappen, dass eine Permutation der Zahlen 1...number_of_choosen_fragments entsteht

	mapNumbers();

	/*cout<<"Y-Permutation (size = "<<permutationY->size() <<"):\n";
	copy(permutationY->begin(),permutationY->end(), ostream_iterator<int>(cout," "));
	cout<<endl;
  
	if (dimensions > 2){
		cout<<"Z-Permutation (size = "<<permutationZ->size() <<"):\n";
		copy(permutationZ->begin(),permutationZ->end(), ostream_iterator<int>(cout," "));
		cout<<endl; 
	}*/
  return 0;
}

int Transformer::mapNumbers(){
  //zuordnung berechnen
	int inputID;
	int realID = 1;
	set<int>::const_iterator myI;
	InputID2ID->clear();
	for(myI = choosenIntervalsIDs->begin(); myI != choosenIntervalsIDs->end() ; ++myI){
		inputID = *myI; 
		(*InputID2ID)[inputID]=realID;
		realID++;
	}
  //zuordnung ausfhren:
	int sign;
	
	for (unsigned int i = 0; i < permutationY->size(); i++){
		inputID = (*permutationY)[i];
		if (inputID < 0) sign = -1; else sign = 1;
		realID = (*InputID2ID)[abs(inputID)];
		(*permutationY)[i] = sign * realID;
	}
	if (dimensions > 2){
		for (unsigned int i = 0; i < permutationZ->size(); i++){
			inputID = (*permutationZ)[i];
			if (inputID < 0) sign = -1; else sign = 1;
			realID = (*InputID2ID)[abs(inputID)];
			(*permutationZ)[i] = sign * realID;
		}
	}
  return 0;
}



int Transformer::getAllNotChoosenIntervalIDs(){
	//allNotChoosenIntervalsIDs = new set<int, less<int> >;
	allNotChoosenIntervalsIDs->clear();
	for (unsigned int i = 0; i < Intervals->size(); i++){
		if (choosenIntervalsIDs->find(i+1) == choosenIntervalsIDs->end()){
			allNotChoosenIntervalsIDs->insert(i+1);
		}
	}
	return 0;
}

int Transformer::printIntervalsForPlot( set<int>* intervalIDs, char* outDatFile, unsigned dim1, unsigned dim2){
	assert(0 <= dim1 < 3);
	assert(0 <= dim2 < 3);
	long int start1, start2, end1, end2;

	ofstream outfilestream(outDatFile, ios::app);
	if (!outfilestream){
		cout<<"cannot open outfile "<<outDatFile<<endl;
		exit(0);
	}
  //intervall grenzen rausschreiben
	set<int>::const_iterator myI;
	int intervalID;
	Interval interval;
	Point* startP;
	Point* endP;
	for (myI = intervalIDs->begin(); myI != intervalIDs->end(); myI++){
		intervalID = *myI; //vorsicht: id ist eins gr�er als entsprchender Index in Intervals
		assert (intervalID-1 < Intervals->size());
		interval = (*Intervals)[intervalID-1];
		startP = interval.begin();
		endP = interval.end();
		start1 = (*startP).coordinates[dim1];
		start2 = (*startP).coordinates[dim2];
		end1 = (*endP).coordinates[dim1];
		end2 = (*endP).coordinates[dim2];
		outfilestream<<start1<<"\t"<<start2<<"\n"<<end1<<"\t"<<end2<<"\n\n";
	}
	outfilestream<<endl<<endl;
	outfilestream.close();
	return 0;
}

int Transformer::printForPlot(char* outDatFile, unsigned dim1, unsigned dim2){
	ofstream outfilestream(outDatFile, ios::out);
	if (!outfilestream){
		cout<<"cannot open outfile "<<outDatFile<<endl;
		exit(0);
	}
	outfilestream<<"#-- Choosen Intervals\n";
	outfilestream.close();
	printIntervalsForPlot(choosenIntervalsIDs, outDatFile, dim1, dim2);
	outfilestream.open(outDatFile, ios::app);
	if (!outfilestream){
		cout<<"cannot open outfile "<<outDatFile<<endl;
		exit(0);
	}
	outfilestream<<"#-- Not Choosen Intervals\n";
	outfilestream.close();
	printIntervalsForPlot(allNotChoosenIntervalsIDs, outDatFile, dim1, dim2);
	return 0;
}

int Transformer::printForPlot(char* filePraefix){
	char completeFileName01[80];
	strcpy(completeFileName01, filePraefix);
	strcat(completeFileName01, "_01.dat");
	printForPlot(completeFileName01, 0,1);	
	if(dimensions > 2){
		char completeFileName02[80];
		strcpy(completeFileName02, filePraefix);
		strcat(completeFileName02, "_02.dat");
		printForPlot(completeFileName02, 0,2);
		char completeFileName12[80];
		strcpy(completeFileName12, filePraefix);
		strcat(completeFileName12, "_12.dat");
		printForPlot(completeFileName12,1,2);
	}
	return 0;
}


int Transformer::printPermutations(){
    ofstream output;    // only for G2 respective G1
    
    // !!DEBUG
    int* debug = new int[permutationY->size() + 1];
    for (int i = 0; i <= permutationY->size(); i++)
        debug[i] = 0;
    for (int i = 0; i < permutationY->size(); i++)
    {
        if ((*permutationY)[i] > 0)
            debug[(*permutationY)[i]]++;
        else
            debug[-(*permutationY)[i]]++;
    }
    for (int i = 1; i <= permutationY->size(); i++)
    {
        if (debug[i] != 1)
            cerr << "error in permutation\n";
    }            
    
    output.open("postprocessed.dat");
    output << "# permutation size\n";
    output << permutationY->size() << endl;
	cout<<"Y-Permutation (second genome with respect to the first one):\n";
    output << "# Y-Permutation (second genome with respect to the first one):\n";
 	copy(permutationY->begin(),permutationY->end(), ostream_iterator<int>(cout," "));
 	copy(permutationY->begin(),permutationY->end(), ostream_iterator<int>(output," "));
 	cout<<endl;
    output << endl;
    output << "# First genome (= id)\n";
    for (int i = 0; i < permutationY->size(); i++)
        output << i+1 << " ";
    output << endl;
    output.close();
    
    
  	if (dimensions > 2){
		cout<<"Z-Permutation (third genome with respect to the first on):\n";
 		copy(permutationZ->begin(),permutationZ->end(), ostream_iterator<int>(cout," "));
 		cout<<endl; 
	}
    
    delete[] debug;
    
	return 0;
}


TransformerWithoutOverlap::TransformerWithoutOverlap(char* datFile, unsigned int dim):Transformer(datFile, dim){
}

TransformerWithoutOverlap::~TransformerWithoutOverlap(){

}

int TransformerWithoutOverlap::chainOneDirection(int dimension){
	//Ergebnis des 1-dimensionalen Chainings wird in (x,y,z)IntervalsIDs und sorted(X,Y,Z)IntervalsIDs gespeichert

	Chain *myChain = new Chain;
	myChain->setDim(dimension);
	Interval *lastInterval = myChain->oneDimensionalChaining(Intervals);

	if(dimension == 0){//x
		xIntervalsIDs->clear(); sortedXIntervalsIDs->clear();
		myChain->getIntervalIDs(xIntervalsIDs,sortedXIntervalsIDs);
    /*cout<<"NUMBER in xIntervalsIDs = "<<xIntervalsIDs->size()<<endl;
    cout<<"NUMBER in sortedXIntervalsIDs = "<<sortedXIntervalsIDs->size()<<endl;*/
	}

	if(dimension == 1){//y
		yIntervalsIDs->clear(); sortedYIntervalsIDs->clear();
		myChain->getIntervalIDs(yIntervalsIDs,sortedYIntervalsIDs);
    /*cout<<"NUMBER in yIntervalsIDs = "<<yIntervalsIDs->size()<<endl;
    cout<<"NUMBER in sortedYIntervalsIDs = "<<sortedYIntervalsIDs->size()<<endl;*/
	}

	if(dimension == 2){//z
		zIntervalsIDs->clear(); sortedZIntervalsIDs->clear();
		myChain->getIntervalIDs(zIntervalsIDs,sortedZIntervalsIDs);
    /*cout<<"NUMBER in zIntervalsIDs = "<<zIntervalsIDs->size()<<endl;
    cout<<"NUMBER in sortedZIntervalsIDs = "<<sortedZIntervalsIDs->size()<<endl;*/
	}
	
	delete myChain;
	return 0;
}

int TransformerWithoutOverlap::chooseFragments(){
	//IDs of choosen Fragments are stored in choosenIntervalsIDs

	//1-dimensionales Chaining in x-Richtung:   
	chainOneDirection(0);
	//1-dimensionales Chaining in y-Richtung: 
	chainOneDirection(1);
	if (dimensions > 2){
		//1-dimensionales Chaining in z-Richtung:
		chainOneDirection(2);
	}

	//Schnittmenge finden
	choosenIntervalsIDs->clear();
	if (dimensions == 2){
		insert_iterator<set<int, less<int> > > res_ins_xy(*choosenIntervalsIDs,choosenIntervalsIDs->begin());
		set_intersection(xIntervalsIDs->begin(), xIntervalsIDs->end(), yIntervalsIDs->begin(), yIntervalsIDs->end(), res_ins_xy);
	}
	if (dimensions ==3){
		set<int, less<int> > *intersectionXYIDs = new set<int, less<int> >;
		insert_iterator<set<int, less<int> > > res_ins_xy(*intersectionXYIDs,intersectionXYIDs->begin());
		set_intersection(xIntervalsIDs->begin(), xIntervalsIDs->end(), yIntervalsIDs->begin(), yIntervalsIDs->end(), res_ins_xy);
		insert_iterator<set<int, less<int> > > res_ins_xyz(*choosenIntervalsIDs,choosenIntervalsIDs->begin());
		set_intersection(intersectionXYIDs->begin(),intersectionXYIDs->end(),zIntervalsIDs->begin(), zIntervalsIDs->end(), res_ins_xyz);
		delete intersectionXYIDs;

  //Indizes finden, die NICHT in der Schnittmenge sind fuer drei Mengen: notChoosenIntervalsIDs  (dafuer alle drei vereinigen und dann die choosenIntervalsIndizes abziehen)
/*		insert_iterator<set<int, less<int> > > res_ins_xy2(*unionXYIDs,unionXYIDs->begin());
		set_union(xIntervalsIDs->begin(), xIntervalsIDs->end(), yIntervalsIDs->begin(), yIntervalsIDs->end(),res_ins_xy2);
		insert_iterator<set<int, less<int> > > res_ins_xyz2(*unionXYZIDs,unionXYZIDs->begin());
		set_union(unionXYIDs->begin(), unionXYIDs->end(),zIntervalsIDs->begin(), zIntervalsIDs->end(), res_ins_xyz2);
		cout<<"NUMBER in unionXYZIDs = "<<unionXYZIDs->size()<<endl;
		insert_iterator<set<int, less<int> > > res_ins_nc(*notChoosenIntervalsIDs,notChoosenIntervalsIDs->begin());
		set_difference(unionXYZIDs->begin(),unionXYZIDs->end(),choosenIntervalsIDs->begin(),choosenIntervalsIDs->end(), res_ins_nc);*/
	}

	getAllNotChoosenIntervalIDs();
	return 0;

}

int TransformerWithoutOverlap::findOverlapsAndChange(unsigned idnc, unsigned idc, float percentage){
	vector<long int> csb(3); vector<long int> cse(3);
//return Interval; contains for each dimension (x,y and z) the overlap or zero if there isn't any overlap
    float percOverlap;
	
	for (unsigned int i = 0; i < dimensions; i++){
		interval_one_dim i1_1dim = (*Intervals)[idnc].getOneDimension(i);
		interval_one_dim i2_1dim = (*Intervals)[idc].getOneDimension(i);
		long int start1, end1, start2, end2;
		//start und endpunkte 
		if(i1_1dim.p1 <= i1_1dim.p2){
			start1 = i1_1dim.p1; end1 = i1_1dim.p2;
		}
		else{
			start1 = i1_1dim.p2; end1 = i1_1dim.p1;
		}
		if(i2_1dim.p1 <= i2_1dim.p2){
			start2 = i2_1dim.p1; end2 = i2_1dim.p2;
		}
		else{
			start2 = i2_1dim.p2; end2 = i2_1dim.p1;
		}
	
		if ((end2 < start1) || (start2 > end1)){//no overlap
			csb[i] = 0; cse[i] = 0;
		}
		else{//overlap
			if (start2 <= start1 && end1 <= end2){
				//i1 komplett umschlossen von i2
				csb[i] = start1; cse[i] = end1;
			}
			if (start1 <= start2 && end2 <= end1){
				//komplett umschlossen andersherum, rest zerfaellt in zwei teile!
				csb[i] = start2; cse[i] = end2;
			}
			if (start2 >= start1 && start2 <= end1 && end2 > end1){
				//Fall 3a
				csb[i] = start2; cse[i] = end1;
				//check size of overlap compared with size of not choosen interval
                percOverlap = (end1-start2+1) / (end1-start1+1);
				if (percOverlap <= percentage){//if there is only a small overlap compared to the size of the interval: resize interval so that there is no overlap any more
                    if (i1_1dim.p1 <= i1_1dim.p2){   // overlapping dimension not reversed -> move endpoints               
    					Point* startP1 = (*Intervals)[idnc].begin();
					    Point* endP1 = (*Intervals)[idnc].end();
    					for(unsigned dim = 0; dim < dimensions; dim++){
                            endP1->coordinates[dim] -= (long int)(percOverlap * (endP1->coordinates[dim] - startP1->coordinates[dim]));
                        }
                    }
                    else{   // overlapping dimension reversed -> move startpoints               
    					Point* startP1 = (*Intervals)[idnc].begin();
					    Point* endP1 = (*Intervals)[idnc].end();
    					for(unsigned dim = 0; dim < dimensions; dim++){
                            startP1->coordinates[dim] += (long int)(percOverlap * (endP1->coordinates[dim] - startP1->coordinates[dim]));
                        }
                    }
				}
			}
			if(start1 >= start2 && start1 <= end2 && end2 < end1){
				//Fall 3b
				csb[i] = start1; cse[i] = end2;
				//check size of overlap compared with size of not choosen interval
                percOverlap = (end2-start1+1) / (end1-start1+1);
				if (percOverlap <= percentage){//if there is only a small overlap compared to the size of the interval: resize interval so that there is no overlap any more
                    if (i1_1dim.p1 <= i1_1dim.p2){   // overlapping dimension not reversed -> move startpoints               
    					Point* startP1 = (*Intervals)[idnc].begin();
					    Point* endP1 = (*Intervals)[idnc].end();
    					for(unsigned dim = 0; dim < dimensions; dim++){
                            startP1->coordinates[dim] += (long int)(percOverlap * (endP1->coordinates[dim] - startP1->coordinates[dim]));
                        }
                    }
                    else{   // overlapping dimension reversed -> move endpoints               
    					Point* startP1 = (*Intervals)[idnc].begin();
					    Point* endP1 = (*Intervals)[idnc].end();
    					for(unsigned dim = 0; dim < dimensions; dim++){
                            endP1->coordinates[dim] -= (long int)(percOverlap * (endP1->coordinates[dim] - startP1->coordinates[dim]));
                        }
                    }
				}
			}
			else {}
		}
	}
	return 0;
}

int TransformerWithoutOverlap::getOverlapsAndChange(float percentage){
	set<int, less<int> > *overlaps_c = new set<int, less<int> >;
	set<int, less<int> > *overlaps_nc = new set<int, less<int> >;

	Interval nci; //not choosen interval
	set<int>::const_iterator myCI;
	int intervalCID;//choosen interval id 
	set<int>::const_iterator ncii;//not choosen interval iterator
	unsigned int intervalID;
	Interval* overlaps;
	for(ncii = allNotChoosenIntervalsIDs->begin(); ncii != allNotChoosenIntervalsIDs->end() ; ncii++){
		intervalID = *ncii; //vorsicht: id ist eins gr�er als entsprchender Index in Intervals
		assert (intervalID-1 < Intervals->size());
		//cout<<"\ntreating not Choosen Interval with id "<<intervalID<<endl;
		nci = (*Intervals)[intervalID-1];
		for(myCI = choosenIntervalsIDs->begin(); myCI != choosenIntervalsIDs->end(); myCI++){
			intervalCID = *myCI; //vorsicht: id ist eins gr�er als entsprchender Index in Intervals
			Interval ci = (*Intervals)[intervalCID-1];
			//if(overlap(&nci, &ci)){
			if(nci.overlap(&ci)){
				//cout<<"found overlap\n";
				overlaps_c->insert(intervalCID);
				findOverlapsAndChange(intervalID-1, intervalCID-1, percentage);
			}
		}
	}
	delete overlaps_c;
	delete overlaps_nc;

	return 0;
}

int TransformerWithoutOverlap::clearData(){
	xIntervalsIDs->clear();
	yIntervalsIDs->clear();
	zIntervalsIDs->clear();
	sortedXIntervalsIDs->clear();
	sortedYIntervalsIDs->clear();
	sortedZIntervalsIDs->clear();
	//choosenIntervalsIDs->clear();
	permutationY->clear();
	permutationZ->clear();
	notInPermutationY->clear();
	notInPermutationZ->clear();
	InputID2ID->clear();
	//allNotChoosenIntervalsIDs->clear();
	return 0;
}

int TransformerWithoutOverlap::transformIt(){
	if (dimensions ==2){
		parseDatFile2Genomes(inDatFile);
	}
	if (dimensions == 3){
		parseDatFile3Genomes(inDatFile);
	}
	//data from inDatFile is now in Intervals
	sortIntervalsXDirection();
	chooseFragments();
	getPermutation();
//	printPermutations();
	printForPlot("choosenAndNotChoosen");
	//Postprocessing:
	//checking if intervals are not choosen, because they overlap only little (<= 10 percent of its size) with others. In this case resize the interval so that it doesn't overlap any more and can be included
	cout<<"\nAfter Postprocessing (before not choosen Intervals, which overlap only little with another one are resized, so they fit in):\n";
	getOverlapsAndChange(0.1);

	sortIntervalsXDirection();
	chooseFragments();
	getPermutation();
	printPermutations();
	printForPlot("choosenAndNotChoosenPostprocessed");

	return 0;
}
