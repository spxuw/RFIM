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
#ifdef _WIN32
	#include <strstream>
#else
	#include <sstream>
#endif
using namespace std;

#include "common.h"
#include "gBase.h"
#include "adj_matrix.h"

void removeFiles(string f);

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

	if( ARGC < 9){
		cout << "feature_atom feature_bond target_atom " << 
			 "target_bond feature_act svm_out binary sub/auto(0/1)\n";
		return 0;
	}	
	
	string f1, f2, f3, f4, f5 , f6, f7;
        float freq = 1.0;
	f1 = ARGV[1];
	f2 = ARGV[2];
	f3 = ARGV[3];
	f4 = ARGV[4];
	f5 = ARGV[5];
        f6 = ARGV[6];
        int binary = atoi(ARGV[7]);
	int sub    = atoi(ARGV[8]);

	cout << "constructing first: binary " << binary << endl;
	gBase * featureDB = new gBase(f1, f2, 0);
	
	cout << "constructing two now:\n";
	gBase * targetDB = new gBase(f3, f4, 0);
	  
	//read in the activity
	char buff[81];
	vector<int> activity;
	char a1[80], a2[80], a3[80];
        vector<int> result;

	ifstream ifile(f5.c_str());
	
	while ( ifile.getline(buff, 80) ){
	  sscanf(buff, "%s %s %s", a1, a2, a3);
	  activity.push_back(atoi(a3));
          //cout << a1 << " " << a2 << " " << a3 << endl;
	}
	ifile.close();

	if(activity.size() != targetDB->size() ){
	  cout << "act size " << activity.size() << " " << featureDB->size() << endl;
	}

	ofstream outfile(f6.c_str());
	outfile.close();
	outfile.open(f6.c_str());

	cout << "check frequency and output SVM format" << 
	  featureDB->size() << " " << targetDB->size() << endl;
        cout.flush();

	int i, j;
        for(i=0; i< targetDB->size(); i++){
            
            if( (i+1) % 10 == 0) cout << "." ;
            if( (i+1) % 100 == 0 ) cout << (i+1) << endl;

            int tc =0;
	    AdjMatrix * thisg = targetDB->graph(i);
	    outfile << activity[i] << " ";
	    
	    int k= 0;
            result.clear();
	    for( j=0; j< featureDB->size(); j++){ 
	       if( !binary ){
		 //if( (k=thisg->findAllIsomophisms( featureDB->graph(j))) ){
                 if ( featureDB->graph(j)->isSubgraphOf(thisg, sub) ){
                   result.push_back(j+1);
		   tc++;
		 }
	       }
	       else{
		 if ( featureDB->graph(j)->isSubgraphOf(thisg, sub) ){
                   result.push_back(j+1);
 		   tc++;
		 }
	       }  
            }
            for(j=0; j< result.size(); j++)
	       outfile << result[j] << ":" << 1 << " ";
	    outfile << endl;
            cout << "total " << tc << endl;
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



