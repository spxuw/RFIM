#pragma warning (disable:4786) 
#pragma warning (disable:4267) 
#pragma warning (disable:4018) 
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <map>
#include <ctime>
#include <iomanip>
#include <stdio.h>
#include <errno.h>
#include <algorithm>
using namespace std;

#include "gBase.h"
#include "adj_matrix.h"

void removeFiles(string f);

class gscomp
{
public:
  bool operator()(  AdjMatrix * const & g1,  AdjMatrix *  const & g2){
    if (g1->size() > g2->size() ) return true;
    else if( g1->size() < g2->size() ) return false;
    else if( g1->esize() > g2->esize() ) return true;
    else return false;
  }
};

class dencomp
{
public:
  bool operator()(  AdjMatrix * const & g1,  AdjMatrix *  const & g2){
     
    float t1 = g1->esize() *1.0 / g1->size();
    float t2 = g2->esize() *1.0 / g2->size();
    if( t1 > t2 ) return true;
    else if( t1 < t2 ) return false;
    else if( g1->size() > g2->size() ) return true;
    return false;
  }
};

int main (int ARGC, char * ARGV[]){

	if( ARGC < 6){
	    cout << "feature_atom feature_bond target_atom " << 
                     "target_bond choice (0: size 1: density)\n";
	    return 0;
	}	
	
	string f1, f2, f3, f4;
        float freq = 1.0;
	f1 = ARGV[1];
	f2 = ARGV[2];
	f3 = ARGV[3];
	f4 = ARGV[4];
	int choice = atoi(ARGV[5]);

	cout << "constructing first\n";
	gBase * featureDB = new gBase(f1, f2, 0);

	ofstream outfile(f3.c_str());
	outfile.close();
	outfile.open(f4.c_str());
	outfile.close();

	vector<AdjMatrix *> results;
        int i;
        for(i=0; i< featureDB->size(); i++){
	  results.push_back(featureDB->graph(i));
	  //cout << "index " << i << " size " << featureDB->graph(i)->size() << " esize "
          //     << featureDB->graph(i)->esize() << endl;
	}
	cout << "end\n" ;
        cout.flush();
	if( choice == 0) 
	  sort(results.begin(), results.end(), gscomp());
	else
	  sort(results.begin(), results.end(), dencomp());

        //output the selected features
	cout << "resutls size" << results.size() << endl;
	for(i=0; i< results.size(); i++){
	  results[i]->print(i, f3, f4);
	}
	return 1;

}

void removeFiles(string f){
#ifdef LINUX
     string com = "rm ";
#else
     string com = "del ";
#endif

     string com1;
     FILE * fh;
     
     fh = fopen (f.c_str(), "r");
     if ( fh!= NULL || (fh == NULL && errno != ENOENT ) ){
        if( fh) fclose(fh);
        com1 = com + f;
        system(com1.c_str());
     }
}



