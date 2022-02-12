/*######################################################################################################################
  Permute the gene-wise pvalue file 
  ######################################################################################################################*/

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cmath>
#include <vector>
#include <list>
#include <map>

#include "find.h"
#include "timer.h"
#include "Splitter.h"


using namespace std;

void Char2String(char& c, string& s) {
  stringstream ss;
  ss << c;
  ss >> s;
}


string Int2String(int i) {
  string s;
  stringstream ss;
  ss << i;
  ss >> s;
  return s;
}


// random generator function:
int myrandom (int i) { return std::rand()%i;}


int main (int argc, char** argv)
{
  if(argc!=4)
    {
      cout <<"\nNumber of arguments is incorrect.";
      cout <<"\nCommand format:  ./Permute pvaluefile seedi seedf\n";
      exit(0);
    }

  string pfile = argv[1];
  int seedi = atoi(argv[2]);
  int seedf = atoi(argv[3]);

  timer time;
  time.start();


  /////////////////////////////////////////////////////////////////////////
  cout << "Read the gene pvalue file: " << pfile  << endl << endl;
  ifstream fin1(pfile.c_str(), ios_base::in);
  if(!fin1) {cout << "Cannot open " << pfile << " for read."; exit(0);}
  
  vector<string> nodenamelist; 
  string name; 
  double pvalue; 
  map<string,int> NodeMAP;
  vector<double> P;

  int count=0;
  int repeat=0;

  string S;
  while(getline(fin1,S)) {
    Splitter split(S, " ");
    int k = (int)split.size();
    if(k==2) {
      name = split[0];
      pvalue = atof(split[1].c_str());
   
      if(name!="") { // ignore those genes with missing names
	if(!find(nodenamelist, name)) 
	  {
	    nodenamelist.push_back(name);
	    NodeMAP[name]=count;
	    count++;
	    P.push_back(pvalue);
	  }
	else {
	  repeat++;
	  cout << name << ' ';
	}
      }
      else 
	cout << "Problematic line: " << S << endl; 
    }
    else 
      cout << "Problematic line: " << S << endl;
  }

  fin1.close();
  int N1= nodenamelist.size();
  cout << "\nIn total there are " << N1 << " nodes in the p-value files.\n";
  if(repeat>0) {
    cout << " And "<< repeat << " of them have repeated/multiple p-values.";
  }
  cout << endl; 


  /////////////////////////////////////////////////////////////////////////
  cout << "Permute the pvalues, and save them to new files : " << endl << endl;

  for(int seed=seedi; seed<=seedf; seed++) {
    cout << seed << endl; //test

    // set up the seed of the random generator std::rand()
    //std::srand(std::time(0)); //use current time as seed for random generator
    std::srand(seed); 

    random_shuffle(P.begin(), P.end(), myrandom);
  
    //cout << "Save the permutated pvalues to a new file:  " << ofile << endl;
    string ofile = pfile + ".seed." + Int2String(seed);

    ofstream fout(ofile.c_str(), ios_base::out);
    if(!fout) {cout << "Cannot open " << ofile  << " for write."; exit(0);}
    for(int i=0;i<N1;i++) {
      fout << nodenamelist[i] << ' ' << P[i] << endl;
    }
    fout.close();
  }


  cout << "\n It takes " << time << " s.\n\n";
  time.restart();

  return 1;
}

