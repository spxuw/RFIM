#include "Parser.h"
#include "Find.h"

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//           The following functions are parsers, which convert different format files into
//           standard elist file, where node index starts from 0, and the first two numbers 
//           represent N=n and E=m, i.e. 
//            N    E
//            i -> j
//            p -> q
//            ......
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////////////////////////////
// Parser 1: convert Elist (2-col, in String format) to Elist (2-col, in Number format)
//
// we have to determine how many unique nodes we have from the Elist string data 
// Here, the input Elist data could contain arbitary string, e.g. gene's name, eprint's ID, etc.
// the output Elist data would be the standard elist format: i --> j
// with i, j are numbers which are the starting and ending nodes' indices of the arc i-->j
// the output Elist Number data is saved in both .txt and .bin formats
//
// Suppose we have a two-column data
// e.g. 
// acrR  acrA 
// acrR  acrB 
// acrR  acrR 
// acrR  micF 
// for this network, we just have 5 nodes and 4 edges
// we have the mapping 
// acrR  -- 0
// acrA  -- 1
// acrB  -- 2
// micF  -- 3
// the standard elist output would be
// 0  1
// 0  2
// 0  0 
// 0  3
// 
// Note that sometimes we have more columns in the string format, 
// in this case we have to prepare a two-column data!!!!!!! This is very important.

void Parser_s2n(char* fname1)
{
  ifstream fin(fname1, ios_base::in);
  if(!fin) {cout << "Cannot open " << fname1 << " for read."; exit(0);}

  char fname2[256];
  sprintf(fname2, "%s.nodemap", fname1);
  ofstream fout(fname2, ios_base::out);
  if(!fout) {cout << "Cannot open " << fname2 << " for write."; exit(0);}

  char fname[256];
  sprintf(fname, "%s.elist.t", fname1);
  ofstream tout(fname, ios_base::out);
  if(!tout) {cout << "Cannot open " << fname << " for write."; exit(0);}

  sprintf(fname, "%s.elist.b", fname1);
  ofstream bout(fname, ios_base::out|ios::binary);
  if(!bout) {cout << "Cannot open " << fname << " for write."; exit(0);}

  // do the index mapping first
  vector<string> nodenamelist; //  name is a string (could be words, not just number)
  string source; 
  string target;
  while (fin >> source >> target) {
    if(!find(nodenamelist, source)) {
      nodenamelist.push_back(source);
    }
    if(!find(nodenamelist, target)) {
      nodenamelist.push_back(target);
    }
  }
  fin.close();

  // save the nodemap file
  int n = nodenamelist.size();
  map<string,int> MAP; // save the map between nodename and the index
  for(int i=0; i<n; i++) {
    fout << i << ' ' << nodenamelist[i] << endl;
    MAP[nodenamelist[i]] = i;
  }

  // read the data again
  ifstream fin2(fname1, ios_base::in);
  if(!fin2) {cout << "Cannot open " << fname1 << " for read."; exit(0);}

  vector< set<int> > aoutlist(n);
  while (fin2 >> source >> target) {
    int i = MAP[source];
    int j = MAP[target];
    aoutlist[i].insert(j); // using set<int> will avoid getting repeated links
  }
  fin2.close();
  
  int m = 0;
  for(int i=0; i<n; i++) {
    m += aoutlist[i].size();
  }
  

  // save the standard edgelist file
  // Note that the first two lines store the information of N and E !!
  cout << fname1 << "\n# N= " << n << ", E= " << m << endl;
  tout << "#N= " << n << endl;
  tout << "#E= " << m << endl;

  int* nm = new int [2];
  nm[0]=n; nm[1]=m;
  bout.write(reinterpret_cast<char *>(nm), 2*sizeof(int));
  delete [] nm;
  int* endpoints = new int [2*m];
  int e=0;
  for(int i=0; i<n; i++) {
    for(set<int>::iterator it=aoutlist[i].begin(); it!=aoutlist[i].end(); it++) {
      int j = (*it);
      endpoints[e*2+0]=i;
      endpoints[e*2+1]=j;
      e++;
      tout << i << ' ' << j << endl;
    }
  }
  bout.write(reinterpret_cast<char *>(endpoints), 2*m*sizeof(int));
  delete [] endpoints;

  fout.close();
  tout.close();
  bout.close();
}
/////////////////////////////////////////////////////////////////////////////////////////////////





/////////////////////////////////////////////////////////////////////////////////////////////////
// same as above, just output the nodename map
void Parser_s2n(char* fname1, map<string,int>& MAP)
{
  ifstream fin(fname1, ios_base::in);
  if(!fin) {cout << "Cannot open " << fname1 << " for read."; exit(0);}

  char fname2[256];
  sprintf(fname2, "%s.nodemap", fname1);
  ofstream fout(fname2, ios_base::out);
  if(!fout) {cout << "Cannot open " << fname2 << " for write."; exit(0);}

  char fname[256];
  sprintf(fname, "%s.elist.t", fname1);
  ofstream tout(fname, ios_base::out);
  if(!tout) {cout << "Cannot open " << fname << " for write."; exit(0);}

  sprintf(fname, "%s.elist.b", fname1);
  ofstream bout(fname, ios_base::out|ios::binary);
  if(!bout) {cout << "Cannot open " << fname << " for write."; exit(0);}

  // do the index mapping first
  vector<string> nodenamelist; //  name is a string (could be words, not just number)
  string source; 
  string target;
  while (fin >> source >> target) {
    if(!find(nodenamelist, source)) {
      nodenamelist.push_back(source);
    }
    if(!find(nodenamelist, target)) {
      nodenamelist.push_back(target);
    }
  }
  fin.close();

  // save the nodemap file
  int n = nodenamelist.size();
  //map<string,int> MAP; // save the map between nodename and the index
  for(int i=0; i<n; i++) {
    fout << i << ' ' << nodenamelist[i] << endl;
    MAP[nodenamelist[i]] = i;
  }

  // read the data again
  ifstream fin2(fname1, ios_base::in);
  if(!fin2) {cout << "Cannot open " << fname1 << " for read."; exit(0);}

  vector< set<int> > aoutlist(n);
  while (fin2 >> source >> target) {
    int i = MAP[source];
    int j = MAP[target];
    aoutlist[i].insert(j); // using set<int> will avoid getting repeated links
  }
  fin2.close();
  
  int m = 0;
  for(int i=0; i<n; i++) {
    m += aoutlist[i].size();
  }
  

  // save the standard edgelist file
  // Note that the first two lines store the information of N and E !!
  cout << fname1 << "\n# N= " << n << ", E= " << m << endl;
  tout << "#N= " << n << endl;
  tout << "#E= " << m << endl;

  int* nm = new int [2];
  nm[0]=n; nm[1]=m;
  bout.write(reinterpret_cast<char *>(nm), 2*sizeof(int));
  delete [] nm;

  int* endpoints = new int [2*m];
  int e=0;
  for(int i=0; i<n; i++) {
    for(set<int>::iterator it=aoutlist[i].begin(); it!=aoutlist[i].end(); it++) {
      int j = (*it);
      endpoints[e*2+0]=i;
      endpoints[e*2+1]=j;
      e++;
      tout << i << ' ' << j << endl;
    }
  }
  bout.write(reinterpret_cast<char *>(endpoints), 2*m*sizeof(int));
  delete [] endpoints;

  fout.close();
  tout.close();
  bout.close();
}
/////////////////////////////////////////////////////////////////////////////////////////////////





/////////////////////////////////////////////////////////////////////////////////////////////////
// deal with the three-column interbank data: 
//e.g.  source     target     money-loan
void Parser_s2n_Weighted(char* fname1)
{
  ifstream fin(fname1, ios_base::in);
  if(!fin) {cout << "Cannot open " << fname1 << " for read."; exit(0);}

  char fname2[256];
  sprintf(fname2, "%s.nodemap", fname1);
  ofstream fout(fname2, ios_base::out);
  if(!fout) {cout << "Cannot open " << fname2 << " for write."; exit(0);}

  char fname[256];
  sprintf(fname, "%s.elist.t", fname1);
  ofstream tout(fname, ios_base::out);
  if(!tout) {cout << "Cannot open " << fname << " for write."; exit(0);}

  sprintf(fname, "%s.elist.b", fname1);
  ofstream bout(fname, ios_base::out|ios::binary);
  if(!bout) {cout << "Cannot open " << fname << " for write."; exit(0);}

  // do the index mapping first
  cout << "Do the index mapping first.\n";
  vector<string> nodenamelist; //  name is a string (could be words, not just number)
  string source; 
  string target;
  double weight;    
  while (fin >> source >> target >> weight) {
    if(!find(nodenamelist, source)) {
      nodenamelist.push_back(source);
    }
    if(!find(nodenamelist, target)) {
      nodenamelist.push_back(target);
    }
  }
  fin.close();

  // save the nodemap file
  cout << "Save the nodemap file.\n";
  int n = nodenamelist.size();
  map<string,int> MAP; // save the map between nodename and the index
  for(int i=0; i<n; i++) {
    fout << i << ' ' << nodenamelist[i] << endl;
    MAP[nodenamelist[i]] = i;
  }


  // read the data again
  cout << "Read the data again.\n";
  ifstream fin2(fname1, ios_base::in);
  if(!fin2) {cout << "Cannot open " << fname1 << " for read."; exit(0);}


  vector< vector<double> > W(n);
  for(int i=0; i<n; i++) {
    W[i].clear();
    W[i].resize(n,0);
  }

  vector< set<int> > aoutlist(n);
  while (fin2 >> source >> target >> weight) {
    int i = MAP[source];
    int j = MAP[target];
    aoutlist[i].insert(j); // using set<int> will avoid getting repeated links

    W[i][j] = weight;
  }
  fin2.close();
  
  int m = 0;
  for(int i=0; i<n; i++) {
    m += aoutlist[i].size();
  }
  

  // save the standard edgelist file
  cout << "Save the standard edgelist file.\n";
  // Note that the first two lines store the information of N and E !!
  cout << fname1 << "\n# N= " << n << ", E= " << m << endl;
  tout << "#N= " << n << endl;
  tout << "#E= " << m << endl;

  int* nm = new int [2];
  nm[0]=n; nm[1]=m;
  bout.write(reinterpret_cast<char *>(nm), 2*sizeof(int));
  delete [] nm;

  //int* endpoints = new int [2*m];
  double* endpoints = new double [3*m];
  int e=0;
  for(int i=0; i<n; i++) {
    for(set<int>::iterator it=aoutlist[i].begin(); it!=aoutlist[i].end(); it++) {
      int j = (*it);
      endpoints[e*3+0]=i;
      endpoints[e*3+1]=j;
      endpoints[e*3+2]=W[i][j];

      e++;
      tout << i << ' ' << j << ' ' << W[i][j] << endl;
    }
  }
  bout.write(reinterpret_cast<char *>(endpoints), 3*m*sizeof(double));
  delete [] endpoints;

  fout.close();
  tout.close();
  bout.close();
}
/////////////////////////////////////////////////////////////////////////////////////////////////




/////////////////////////////////////////////////////////////////////////////////////////////////
// deal with the three-column flux-coupling data: 
//e.g.  source-reaction (vi)     target-reation (vj)     edge-type
/*
The edge type corresponds to the type of coupling: 
1 = fully coupled:     vi is proportional to vj and vice versa:   vi <===> vj
2 = partially coupled: if vi=0 then vj=0 and vice versa:          vi <---> vj 
4 = directionally coupled: if vi=0, then vj=0                     vi  ---> vj       
*/

void Parser_s2n_FC(char* fname1)
{
  ifstream fin(fname1, ios_base::in);
  if(!fin) {cout << "Cannot open " << fname1 << " for read."; exit(0);}

  char fname2[256];
  sprintf(fname2, "%s.nodemap", fname1);
  ofstream fout(fname2, ios_base::out);
  if(!fout) {cout << "Cannot open " << fname2 << " for write."; exit(0);}

  char fname[256];
  sprintf(fname, "%s.elist.t", fname1);
  ofstream tout(fname, ios_base::out);
  if(!tout) {cout << "Cannot open " << fname << " for write."; exit(0);}

  sprintf(fname, "%s.elist.b", fname1);
  ofstream bout(fname, ios_base::out|ios::binary);
  if(!bout) {cout << "Cannot open " << fname << " for write."; exit(0);}

  // do the index mapping first
  cout << "Do the index mapping first.\n";
  vector<string> nodenamelist; //  name is a string (could be words, not just number)
  string source; 
  string target;
  int edgetype;    
  while (fin >> source >> target >> edgetype) {
    if(!find(nodenamelist, source)) {
      nodenamelist.push_back(source);
    }
    if(!find(nodenamelist, target)) {
      nodenamelist.push_back(target);
    }
  }
  fin.close();

  // save the nodemap file
  cout << "Save the nodemap file.\n";
  int n = nodenamelist.size();
  map<string,int> MAP; // save the map between nodename and the index
  for(int i=0; i<n; i++) {
    fout << i << ' ' << nodenamelist[i] << endl;
    MAP[nodenamelist[i]] = i;
  }


  // read the data again
  cout << "Read the data again.\n";
  ifstream fin2(fname1, ios_base::in);
  if(!fin2) {cout << "Cannot open " << fname1 << " for read."; exit(0);}

  vector< vector<int> > T(n);
  for(int i=0; i<n; i++) {
    T[i].clear();
    T[i].resize(n,0);
  }


  vector< set<int> > aoutlist(n);// using set<int> will avoid getting repeated links   
  while (fin2 >> source >> target >> edgetype) {
    int i = MAP[source];
    int j = MAP[target];
    T[i][j] = edgetype;

    switch(edgetype) {
    case 1 :     
    case 2 :
      aoutlist[i].insert(j);
      aoutlist[j].insert(i);
      break;
    case 4 :
      aoutlist[i].insert(j);
      break;
    default:
      cout << "edgetype should be 1,2, or 4.\n";
      exit(0);
    }
   
  }
  fin2.close();
  
  int m = 0;
  for(int i=0; i<n; i++) {
    m += aoutlist[i].size();
  }
  

  // save the standard edgelist file
  cout << "Save the standard edgelist file.\n";
  // Note that the first two lines store the information of N and E !!
  cout << fname1 << "\n# N= " << n << ", E= " << m << endl;
  tout << "#N= " << n << endl;
  tout << "#E= " << m << endl;

  int* nm = new int [2];
  nm[0]=n; nm[1]=m;
  bout.write(reinterpret_cast<char *>(nm), 2*sizeof(int));
  delete [] nm;

  //int* endpoints = new int [2*m];
  int* endpoints = new int [3*m];
  int e=0;
  for(int i=0; i<n; i++) {
    for(set<int>::iterator it=aoutlist[i].begin(); it!=aoutlist[i].end(); it++) {
      int j = (*it);
      endpoints[e*3+0]=i;
      endpoints[e*3+1]=j;
      endpoints[e*3+2]=T[i][j];

      e++;
      tout << i << ' ' << j << ' ' << T[i][j] << endl;
    }
  }
  bout.write(reinterpret_cast<char *>(endpoints), 3*m*sizeof(int));
  delete [] endpoints;

  fout.close();
  tout.close();
  bout.close();
}
/////////////////////////////////////////////////////////////////////////////////////////////////





/////////////////////////////////////////////////////////////////////////////////////////////////
// deal with the three-column weighted data: 
//      inventor1  inventor2  # of collaborations
// In case the network is too big, but the nodes are already represented by positive integers
// with some gaps in the indices, we may wanna skip the mapping procedure, which is quite time consuming.
// We just use the original index, and respect the gaps, which will be understood as isoated nodes.
void Parser_n2n_Undirected_Weighted(char* fname1)
{
  ifstream fin(fname1, ios_base::in);
  if(!fin) {cout << "Cannot open " << fname1 << " for read."; exit(0);}

  char fname2[256];
  sprintf(fname2, "%s.node.KW", fname1);
  ofstream fout(fname2, ios_base::out);
  if(!fout) {cout << "Cannot open " << fname2 << " for write."; exit(0);}

  char fname[256];
  sprintf(fname, "%s.elist.t", fname1);
  ofstream tout(fname, ios_base::out);
  if(!tout) {cout << "Cannot open " << fname << " for write."; exit(0);}

  sprintf(fname, "%s.elist.b", fname1);
  ofstream bout(fname, ios_base::out|ios::binary);
  if(!bout) {cout << "Cannot open " << fname << " for write."; exit(0);}

  // find the min and max of indices first 
  cout << "Find the min and max of the indices first.\n";
  int source; 
  int target;
  double weight;    

  int max=-1;
  int min=100000000;
  while (fin >> source >> target >> weight) {
    if(source>max) max = source;
    if(target>max) max = target;

    if(source<min) min = source;
    if(target<min) min = target;
  }
  fin.close();
  cout << " index min= " << min << " index max= " << max << endl;
       
  int n = max+1;

  // read the data again
  cout << "Read the data again to calculate the degree and weight for each node.\n";
  ifstream fin2(fname1, ios_base::in);
  if(!fin2) {cout << "Cannot open " << fname1 << " for read."; exit(0);}

  // weight of each node
  vector<double> W(n,0);
  // degree of each node
  vector<double> K(n,0);
  vector< set<int> > aoutlist(n);
  int E=0;
  while (fin2 >> source >> target >> weight) {
    int i = source;
    int j = target;
    aoutlist[i].insert(j); // using set<int> will avoid getting repeated links

    W[i] += weight;
    W[j] += weight;

    K[i]++;
    K[j]++;

    E++;
  }
  fin2.close();


  cout << "Read the data again to calculate the degree correlation, weight correlation.\n";
  double Ksum1 = 0;
  double Ksum2 = 0;
  double Ksum3 = 0;

  double Wsum1 = 0;
  double Wsum2 = 0;
  double Wsum3 = 0;

  ifstream fin3(fname1, ios_base::in);
  if(!fin3) {cout << "Cannot open " << fname1 << " for read."; exit(0);}
  while (fin3 >> source >> target >> weight) {
    int je = K[source]-1;
    int ke = K[target]-1;

    Ksum1 += (je*ke);
    Ksum2 += (je+ke);
    Ksum3 += (je*je + ke*ke);


    double jw = W[source];
    double kw = W[target];

    Wsum1 += (jw*kw);
    Wsum2 += (jw+kw);
    Wsum3 += (jw*jw + kw*kw);
  }

  Ksum1 /= E;
  Ksum2 /= (2*E);
  Ksum3 /= (2*E);
  Ksum2 *= Ksum2;
  double rd = (Ksum1-Ksum2)/(Ksum3-Ksum2);
  cout << "\n (undirected) degree correlation coefficient r_d = " << rd << endl;


  Wsum1 /= E;
  Wsum2 /= (2*E);
  Wsum3 /= (2*E);
  Wsum2 *= Wsum2;
  double rw = (Wsum1-Wsum2)/(Wsum3-Wsum2);
  cout << "\n (undirected) weight correlation coefficient r_w = " << rw << endl;
 




  // save the node weight file
  cout << "Save the node-degree-weight file.\n";
  int Nni = 0; // # of non-isolated nodes 
  int m = 0;

  double Kmean_ni = 0;
  double Wmean_ni = 0;

  for(int i=0; i<n; i++) {
    m += aoutlist[i].size();

    if(W[i]>0) {
      Nni++;
      Kmean_ni += K[i];
      Wmean_ni += W[i];
      fout << i << ' ' << K[i] << ' ' << W[i] << endl;
    }
  }
  Kmean_ni /= Nni;
  Wmean_ni /= Nni;
  cout << "--------------------------------------------------------\n";
  cout << "There are actually " << Nni << " nodes and " << m << " edges.\n";
  cout << "mean degree = " << Kmean_ni << endl;
  cout << "mean weight = " << Wmean_ni << endl;
  cout << "--------------------------------------------------------\n";
  
  // save the standard edgelist file
  cout << "Save the standard edgelist file.\n";
  // Note that the first two lines store the information of N and E !!
  cout << fname1 << "\n# N= " << n << ", E= " << m << endl;
  tout << "#N= " << n << endl;
  tout << "#E= " << m << endl;

  int* nm = new int [2];
  nm[0]=n; nm[1]=m;
  bout.write(reinterpret_cast<char *>(nm), 2*sizeof(int));
  delete [] nm;

  int* endpoints = new int [2*m];
  int e=0;
  for(int i=0; i<n; i++) {
    for(set<int>::iterator it=aoutlist[i].begin(); it!=aoutlist[i].end(); it++) {
      int j = (*it);
      endpoints[e*2+0]=i;
      endpoints[e*2+1]=j;
      e++;
      tout << i << ' ' << j << endl;
    }
  }
  bout.write(reinterpret_cast<char *>(endpoints), 2*m*sizeof(double));
  delete [] endpoints;

  fout.close();
  tout.close();
  bout.close();
       
}
/////////////////////////////////////////////////////////////////////////////////////////////////








/////////////////////////////////////////////////////////////////////////////////////////////////
// Parser 2: convert gml to Elist (in Number format)
/*
 edge
  [
    source 0
    target 1
    value 1
  ]
  edge
  [
    source 0
    target 2
    value 2
  ]
*/


/////////////////////////////////////////////////////////////////////////////////////////////////
const string numbers="0123456789";
// get the number from a string
void GetNumberFromaString(string str, int& x) 
{
  // x = atoi((str.substr(11)).c_str()); // this just works if we know exactly the posisition of the first number
  // in general, we have to do the following:
  int first = str.find_first_of(numbers);
  if(first == string::npos) 
    cout<<"find no numbers"<<endl;

  int last = str.find_last_of(numbers);
  if(last == string::npos) 
    cout<<"find no numbers"<<endl;

  // the following statements have the same function: convert string to int
  // x = atoi((str.substr(first, last-first+1)).c_str());
  istringstream buffer(str.substr(first, last-first+1)); 
  buffer >> x;
  
}
/////////////////////////////////////////////////////////////////////////////////////////////////


void Parser_gml(char* fname1)
{
  ifstream fin(fname1, ios_base::in);
  if(!fin) {cout << "Cannot open " << fname1 << " for read."; exit(0);}

  char fname2[256];
  sprintf(fname2, "%s.nodemap", fname1);
  ofstream fout(fname2, ios_base::out);
  if(!fout) {cout << "Cannot open " << fname2 << " for write."; exit(0);}

  char fname[256];
  sprintf(fname, "%s.elist.t", fname1);
  ofstream tout(fname, ios_base::out);
  if(!tout) {cout << "Cannot open " << fname << " for write."; exit(0);}

  sprintf(fname, "%s.elist.b", fname1);
  ofstream bout(fname, ios_base::out|ios::binary);
  if(!bout) {cout << "Cannot open " << fname << " for write."; exit(0);}


  //vector<string> nodenamelist; //  name is a string 
  //vector < vector<string> > edgelist;
  vector<int> nodelabellist; //  label is a number 
  vector < vector<int> > edgelist;

  int source; 
  int target;
  int value;
  string str;
  string::size_type pos;

  int m=0;
  while (!fin.eof()) {
    getline(fin, str);
    pos = str.find("source");

    if(pos!=string::npos && str.find_first_of(numbers)!=string::npos) { 
      GetNumberFromaString(str, source);
      
      while (!fin.eof()) {
	getline(fin, str);     
	pos = str.find("target");
	if(pos!=string::npos) { 
	  GetNumberFromaString(str, target);
	  //cout << source << "--->" << target << endl;

	  vector<int> edge(2);
	  edge[0]=source;
	  edge[1]=target;
	  edgelist.push_back(edge);
	  m++;

	  if(!find(nodelabellist, source)) {
	    nodelabellist.push_back(source);
	  }
	  if(!find(nodelabellist, target)) {
	    nodelabellist.push_back(target);
	  }
	  
	  break;
	}
      } 
    } // end of if finding a line containing "source"
  }// end of reading
  fin.close();
  
  int n = nodelabellist.size();
  for(int i=0; i<n; i++) {
    fout << i << ' ' << nodelabellist[i] << endl;
  }


  // save the standard edgelist file
  // Note that the first two lines stores the information of N and E !!
  cout << fname1 << "\n# N= " << n << ", E= " << m << endl;
  tout << "#N= " << n << endl;
  tout << "#E= " << m << endl;

  int* nm = new int [2];
  nm[0]=n; nm[1]=m;
  bout.write(reinterpret_cast<char *>(nm), 2*sizeof(int));
  delete [] nm;

  int* endpoints = new int [2*m];
  for(int e=0;e<m;e++) {
    int i= findpos(nodelabellist, edgelist[e][0]); // source
    int j= findpos(nodelabellist, edgelist[e][1]); // target
    endpoints[e*2+0]=i;
    endpoints[e*2+1]=j;
    tout << i << ' ' << j << endl;
  }
  bout.write(reinterpret_cast<char *>(endpoints), 2*m*sizeof(int));
  delete [] endpoints;

  fout.close();
  tout.close();
  bout.close();
}
/////////////////////////////////////////////////////////////////////////////////////////////////




/////////////////////////////////////////////////////////////////////////////////////////////////
// Parser 3: convert paj to Elist (in Number format)
// Note that in .paj file: 1. the node index starts from 1
//                         2. the node map has already be built 
/*
  ......
  20 "Suspended POC"
  21 "Sedimented POC"
  22 "Input"
  23 "Output"
  24 "Respiration"
*arcs
  22   1 552615.0
  22   2 552615.0
  ...... 
*/

void Parser_paj(char* fname1)
{
  ifstream fin(fname1, ios_base::in);
  if(!fin) {cout << "Cannot open " << fname1 << " for read."; exit(0);}

  char fname2[256];
  sprintf(fname2, "%s.temp", fname1);
  ofstream fout(fname2, ios_base::out);
  if(!fout) {cout << "Cannot open " << fname2 << " for write."; exit(0);}

  char fname[256];
  sprintf(fname, "%s.elist.t", fname1);
  ofstream tout(fname, ios_base::out);
  if(!tout) {cout << "Cannot open " << fname << " for write."; exit(0);}

  sprintf(fname, "%s.elist.b", fname1);
  ofstream bout(fname, ios_base::out|ios::binary);
  if(!bout) {cout << "Cannot open " << fname << " for write."; exit(0);}

  string str;
  string::size_type pos;
  string::size_type pos1;

  bool findarc = false;
  int m = 0;
  int n = 0;
  while (!fin.eof() && !findarc) {
    getline(fin, str);
    pos1 = str.find("*vertices");
    pos = str.find("arcs");

    if(pos1!=string::npos) { 
      GetNumberFromaString(str, n);
    }

    if(pos!=string::npos) { 
      //cout << str << endl;
      findarc = true;

      while (!fin.eof()) {
	getline(fin, str);     
	//cout << str.length() << endl; // debug
	if(str.length()>2) { // we find empty line sometime will str.length()==1, why?
	  m++;
	  //cout << str << endl;
	  fout << str << endl;
	}
	else {
	  break;
	}
      }
    } // end of if find arcs
  }// end of reading

  fin.close();
  fout.close();

  
  // get the standard Elist file (.txt)
  // Note that the first two lines store the information of N and E !!
  cout << fname1 << "\n# N= " << n << ", E= " << m << endl;
  tout << "#N= " << n << endl;
  tout << "#E= " << m << endl;
  char cmd [256];
  sprintf(cmd, "awk '{print $1-1, $2-1}' %s >> %s.elist.t", fname2,fname1);
  // here, note that "$1-1, $2-1" is due to that in .paj file, node index starts from 1
  system(cmd);
  tout.close();

  // get the standard Elist file (.bin) from the .txt file
  sprintf(fname, "%s.elist.t", fname1);
  ifstream tin(fname, ios_base::in);
  if(!tin) {cout << "Cannot open " << fname << " for read."; exit(0);}

  int* nm = new int [2];
  nm[0]=n; nm[1]=m;
  bout.write(reinterpret_cast<char *>(nm), 2*sizeof(int));
  delete [] nm;

  int source;
  int target;
  int* endpoints = new int [2*m];
  int e = 0;
  getline(tin, str); cout << str << endl;
  getline(tin, str); cout << str << endl;
  while (tin>> source >> target) {
    endpoints[e*2+0]=source;
    endpoints[e*2+1]=target;
    e++;
  }
  tin.close();

  bout.write(reinterpret_cast<char *>(endpoints), 2*m*sizeof(int));
  delete [] endpoints;

  bout.close();
}
/////////////////////////////////////////////////////////////////////////////////////////////////




/////////////////////////////////////////////////////////////////////////////////////////////////
// Parser 4: convert Elist.t to Elist.b
void Parser_t2b(char* fname1)
{
  // get the standard Elist file (.bin) from the .txt file
  char fname[256];
  sprintf(fname, "%s", fname1);
  ifstream tin(fname, ios_base::in);
  if(!tin) {cout << "Cannot open " << fname << " for read."; exit(0);}

  string str;
  string::size_type pos;

  int n = 0;
  getline(tin, str);
  pos = str.find("N");
  if(pos!=string::npos) { 
    GetNumberFromaString(str, n);
  }

  int m = 0;
  getline(tin, str);
  pos = str.find("E");
  if(pos!=string::npos) { 
    GetNumberFromaString(str, m);
  }

  cout << "#N= " << n << endl;
  cout << "#E= " << m << endl;

  sprintf(fname, "%s.elist.b", fname1);
  ofstream bout(fname, ios_base::out|ios::binary);
  if(!bout) {cout << "Cannot open " << fname << " for write."; exit(0);}

  int* nm = new int [2];
  nm[0]=n; nm[1]=m;
  bout.write(reinterpret_cast<char *>(nm), 2*sizeof(int));
  delete [] nm;

  int source;
  int target;
  int* endpoints = new int [2*m];
  int e = 0;
  while (tin>> source >> target) {
    endpoints[e*2+0]=source;
    endpoints[e*2+1]=target;
    e++;
  }
  tin.close();

  bout.write(reinterpret_cast<char *>(endpoints), 2*m*sizeof(int));
  delete [] endpoints;

  bout.close();
}
/////////////////////////////////////////////////////////////////////////////////////////////////




