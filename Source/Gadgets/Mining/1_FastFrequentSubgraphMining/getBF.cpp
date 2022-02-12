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

#include "common.h"
#include "gBase.h"
#include "adj_matrix.h"

//void removeFiles(string f);

class gpcomp
{
public:
  bool operator()(  AdjMatrix * const & g1,  AdjMatrix *  const & g2){
    if (g1->size() > g2->size() ) return true;
    else if( g1->size() < g2->size() ) return false;
    //else if( g1->esize() >= g2->esize() ) return true;
    else return false;
  }
};

int main (int ARGC, char * ARGV[]){

	if( ARGC < 11){
		cout << "feature_atom feature_bond feature_info target_atom target_bond " <<
			"select_atom select_bond  select_info select_stat (optional) uplimit \n";
		return 0;
	}	
	
	string f1, f2, f3, f4, f5 , f6, f7, f8, f9, f10;
    float freq = 1.0;
	f1 = ARGV[1];
	f2 = ARGV[2];
	f3 = ARGV[3];
	f4 = ARGV[4];
	f5 = ARGV[5];
    f6 = ARGV[6];
    f7 = ARGV[7];
	f8 = ARGV[8];
	f9 = ARGV[9];
	freq = atof(ARGV[10]);


	int totalbins = 10;
	vector<int> hisg(totalbins+1, 0);
	
        
	cout << "constructing first\n" << f1 << f2 << endl;
	gBase * featureDB = new gBase(f1, f2, 0);
	featureDB->setAdditionalInfo(f3);

	cout << "constructing two now:\n";
        cout.flush();
	gBase * targetDB = new gBase(f4, f5, 0);

    //removeFiles(f5); removeFiles(f6);
	//remove f5. f6;
	ofstream outfile(f6.c_str());
	outfile.close();
	outfile.open(f7.c_str());
	outfile.close();
	outfile.open(f8.c_str());
	outfile.close();

	//get the feature vector 
    int uplimit = (int) (targetDB->size() * freq);
	cout << "check background frequency feature size: " << 
		featureDB->size() << " database size: " << targetDB->size() << 
        " uplimit: " << uplimit << endl;
    cout.flush();

	vector<int> selected;
	GRAPHS1 matches;

	int i, j;
    for(i=0; i< featureDB->size(); i++){
            
        if( (i+1) % 10 == 0) cout << "." ;
        if( (i+1) % 100 == 0 ) cout << (i+1) << endl;

        int tc =0;
		AdjMatrix * thisg = featureDB->graph(i);
		matches.clear();
		for( j=0; j< targetDB->size() && tc <= uplimit; j++){
          int k = thisg->isSubgraphOf(targetDB->graph(j), 0);
		  if( k ){
			  matches.push_back(j);
			  tc++; 
		  }
		}

		//cout << "set matches " << endl;
		thisg->setBackgroundF(tc);
		thisg->setTIDs(matches);
        if( tc <= uplimit ){
			selected.push_back(i);
			//cout << "select i " << i << endl;
			int index =0;
			if (tc == 0) index = 0;
			else index = (int) (tc / (targetDB->size()*freq*1.0/totalbins)); 
			if( index > totalbins) index = totalbins;
			hisg[index]++;
	    }
	    //cout << "i " << i << "tc is " << tc << " uplimit " << uplimit << endl;
	}

	//output the statistics
	int sum=0;
    outfile.open(f9.c_str());
    for( i=0; i<=totalbins; i++){
	  sum += hisg[i];
	  outfile << sum << " " ;
	}
	outfile << endl;
	outfile.close();

    //sort the graph by size, edge number, support, background support
    vector<AdjMatrix *> results;
	cout << "total selected features " << selected.size() << endl;

    for(i=0; i< selected.size(); i++){
	  results.push_back(featureDB->graph(selected[i]));
	}
	//random shufflying
	random_shuffle(results.begin(), results.end());
	sort(results.begin(), results.end(), gpcomp());

    //output the selected features
	cout << "resutls size" << results.size() << endl;
	map<int, int> gMap;
	targetDB->getGIndex(gMap);
	for(i=0; i< results.size(); i++){
	     if( results[i]->size() >= 4)
		results[i]->print(i, f6, f7, f8, gMap);
	}
	return 1;

}
