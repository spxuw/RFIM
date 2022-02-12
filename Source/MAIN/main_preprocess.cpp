/*######################################################################################################################
  Pre-process the node-weight file and edge-weight file into a standard single file in the following format 

  #N E pmin p1min pmax wmin wmax
  N E pmin p1min pmax wmin wmax
  #nodeindex nodename nodeweight
  .
  .
  .
  #source target weight
  .
  .
  .

  here   
  N      # of nodes
  E      # of edges
  pmin   minimum p-value 
  p1min  minimum non-zero p-value 
  pmax   maximum p-value
  wmin   minimum edge weight
  wmax   maximum edge weight

  nodeweight = p-value
  edgeweight = confidence level
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

int main (int argc, char** argv)
{
  if(argc!=4)
    {
      cout <<"\nNumber of arguments is incorrect.";
      cout <<"\nCommand format:  ./PreProcess pvaluefile edgelistfile outputfile\n";
      exit(0);
    }

  string pfile = argv[1];
  string efile = argv[2];
  string ofile = argv[3];

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
    
      //while (fin1 >> name >> pvalue) {
      //if(!find(nodenamelist, name)) 

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
  cout << "Read the weighted PPI edgelist file: " << efile << endl << endl;
  ifstream fin2(efile.c_str(), ios_base::in);
  if(!fin2) {cout << "Cannot open " << efile  << " for read."; exit(0);}

  string source, target;
  double weight;  // Note that weight = confidence level
  double wmin = 1e300;
  double wmax = -1e300;
  int E=0;
  int Nselfloop=0;

  while (fin2 >> source >> target >> weight) {
    if(source==target) { 
      cout << source << "--" << target << "; "; 
      Nselfloop++; 
    }
    else {
      E++;
      if(weight<wmin) wmin = weight;
      if(weight>wmax) wmax = weight;

      if(!find(nodenamelist, source)) { 
	nodenamelist.push_back(source);
	NodeMAP[source]=count;
	count++;
	P.push_back(1.0);
      }
      if(!find(nodenamelist, target)) { 
	nodenamelist.push_back(target);
	NodeMAP[target]=count;
	count++;
	P.push_back(1.0);
      }
    }
  }
  fin2.close();

  if(Nselfloop>0) {
    cout << "\nIn total " << Nselfloop << " self-loops are ignored.\n\n";
  }
 
  int N = nodenamelist.size();
  cout << "There are " << N << " nodes: " << N1 << " of them have pre-assigned pvalues. Others are assigned pvalue=1.\n";
  double pmin =1e300;  // the smallest pvalue
  double p1min=1e300;  // the smallest non-zero pvalue
  double pmax =-1e300;
  for(int i=0;i<N;i++) {
    if(P[i]<pmin) 
      pmin = P[i];
    if(P[i]<p1min && P[i] >0) 
      p1min = P[i];
    if(P[i]>pmax)
      pmax = P[i];
  }
  cout << "pmin= " << pmin << endl;
  cout << "p1min= " << p1min << endl;
  cout << "pmax= " << pmax << endl << endl;
 
  cout << "There are " << E << " edges.\n";
  cout << "wmin= " << wmin << endl;
  cout << "wmax= " << wmax << endl << endl;

 

  ////////////////////////////////////////////////////////////////////////////////////////////////////////
  cout << "Save the standard network file to " << ofile << endl;
  ofstream fout(ofile.c_str(), ios_base::out);
  if(!fout) {cout << "Cannot open " << ofile  << " for write."; exit(0);}
  fout << "#N E pmin p1min pmax wmin wmax\n";
  fout << N << ' ' << E << ' ' << pmin << ' ' << p1min << ' ' << pmax << ' ' << wmin << ' ' << wmax << endl;
  fout << "#nodeindex nodename nodeweight\n";
  for(int i=0;i<N;i++) {
    fout << i << ' ' << nodenamelist[i] << ' ' << P[i] << endl;
  }

  ifstream fin3(efile.c_str(), ios_base::in);
  if(!fin3) {cout << "Cannot open " << efile  << " for read."; exit(0);}
  fout << "#source target edgeweight\n";
  while (fin3 >> source >> target >> weight) {
    if(source!=target) { 
      int i = NodeMAP[source];
      int j = NodeMAP[target];
      double wij = weight;
      fout << i << ' ' << j << ' ' << wij << endl;
    }
  } 
  fout.close();

  cout << "\n It takes " << time << " s.\n\n";
  time.restart();

  return 1;
}

