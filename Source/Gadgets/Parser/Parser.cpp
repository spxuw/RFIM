#include "Parser.h"
#include "Find.h"

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

  
  vector<string> nodenamelist; 
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

  
  int n = nodenamelist.size();
  map<string,int> MAP; 
  for(int i=0; i<n; i++) {
    fout << i << ' ' << nodenamelist[i] << endl;
    MAP[nodenamelist[i]] = i;
  }

  
  ifstream fin2(fname1, ios_base::in);
  if(!fin2) {cout << "Cannot open " << fname1 << " for read."; exit(0);}

  vector< set<int> > aoutlist(n);
  while (fin2 >> source >> target) {
    int i = MAP[source];
    int j = MAP[target];
    aoutlist[i].insert(j); 
  }
  fin2.close();
  
  int m = 0;
  for(int i=0; i<n; i++) {
    m += aoutlist[i].size();
  }
  

  
  
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

  
  vector<string> nodenamelist; 
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

  
  int n = nodenamelist.size();
  
  for(int i=0; i<n; i++) {
    fout << i << ' ' << nodenamelist[i] << endl;
    MAP[nodenamelist[i]] = i;
  }

  
  ifstream fin2(fname1, ios_base::in);
  if(!fin2) {cout << "Cannot open " << fname1 << " for read."; exit(0);}

  vector< set<int> > aoutlist(n);
  while (fin2 >> source >> target) {
    int i = MAP[source];
    int j = MAP[target];
    aoutlist[i].insert(j); 
  }
  fin2.close();
  
  int m = 0;
  for(int i=0; i<n; i++) {
    m += aoutlist[i].size();
  }
  

  
  
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

  
  cout << "Do the index mapping first.\n";
  vector<string> nodenamelist; 
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

  
  cout << "Save the nodemap file.\n";
  int n = nodenamelist.size();
  map<string,int> MAP; 
  for(int i=0; i<n; i++) {
    fout << i << ' ' << nodenamelist[i] << endl;
    MAP[nodenamelist[i]] = i;
  }


  
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
    aoutlist[i].insert(j); 

    W[i][j] = weight;
  }
  fin2.close();
  
  int m = 0;
  for(int i=0; i<n; i++) {
    m += aoutlist[i].size();
  }
  

  
  cout << "Save the standard edgelist file.\n";
  
  cout << fname1 << "\n# N= " << n << ", E= " << m << endl;
  tout << "#N= " << n << endl;
  tout << "#E= " << m << endl;

  int* nm = new int [2];
  nm[0]=n; nm[1]=m;
  bout.write(reinterpret_cast<char *>(nm), 2*sizeof(int));
  delete [] nm;

  
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

  
  cout << "Do the index mapping first.\n";
  vector<string> nodenamelist; 
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

  
  cout << "Save the nodemap file.\n";
  int n = nodenamelist.size();
  map<string,int> MAP; 
  for(int i=0; i<n; i++) {
    fout << i << ' ' << nodenamelist[i] << endl;
    MAP[nodenamelist[i]] = i;
  }


  
  cout << "Read the data again.\n";
  ifstream fin2(fname1, ios_base::in);
  if(!fin2) {cout << "Cannot open " << fname1 << " for read."; exit(0);}

  vector< vector<int> > T(n);
  for(int i=0; i<n; i++) {
    T[i].clear();
    T[i].resize(n,0);
  }


  vector< set<int> > aoutlist(n);
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
  

  
  cout << "Save the standard edgelist file.\n";
  
  cout << fname1 << "\n# N= " << n << ", E= " << m << endl;
  tout << "#N= " << n << endl;
  tout << "#E= " << m << endl;

  int* nm = new int [2];
  nm[0]=n; nm[1]=m;
  bout.write(reinterpret_cast<char *>(nm), 2*sizeof(int));
  delete [] nm;

  
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

  
  cout << "Read the data again to calculate the degree and weight for each node.\n";
  ifstream fin2(fname1, ios_base::in);
  if(!fin2) {cout << "Cannot open " << fname1 << " for read."; exit(0);}

  
  vector<double> W(n,0);
  
  vector<double> K(n,0);
  vector< set<int> > aoutlist(n);
  int E=0;
  while (fin2 >> source >> target >> weight) {
    int i = source;
    int j = target;
    aoutlist[i].insert(j); 

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
 




  
  cout << "Save the node-degree-weight file.\n";
  int Nni = 0; 
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
  
  
  cout << "Save the standard edgelist file.\n";
  
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















const string numbers="0123456789";

void GetNumberFromaString(string str, int& x) 
{
  
  
  int first = str.find_first_of(numbers);
  if(first == string::npos) 
    cout<<"find no numbers"<<endl;

  int last = str.find_last_of(numbers);
  if(last == string::npos) 
    cout<<"find no numbers"<<endl;

  
  
  istringstream buffer(str.substr(first, last-first+1)); 
  buffer >> x;
  
}



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


  
  
  vector<int> nodelabellist; 
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
    } 
  }
  fin.close();
  
  int n = nodelabellist.size();
  for(int i=0; i<n; i++) {
    fout << i << ' ' << nodelabellist[i] << endl;
  }


  
  
  cout << fname1 << "\n# N= " << n << ", E= " << m << endl;
  tout << "#N= " << n << endl;
  tout << "#E= " << m << endl;

  int* nm = new int [2];
  nm[0]=n; nm[1]=m;
  bout.write(reinterpret_cast<char *>(nm), 2*sizeof(int));
  delete [] nm;

  int* endpoints = new int [2*m];
  for(int e=0;e<m;e++) {
    int i= findpos(nodelabellist, edgelist[e][0]); 
    int j= findpos(nodelabellist, edgelist[e][1]); 
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
      
      findarc = true;

      while (!fin.eof()) {
	getline(fin, str);     
	
	if(str.length()>2) { 
	  m++;
	  
	  fout << str << endl;
	}
	else {
	  break;
	}
      }
    } 
  }

  fin.close();
  fout.close();

  
  
  
  cout << fname1 << "\n# N= " << n << ", E= " << m << endl;
  tout << "#N= " << n << endl;
  tout << "#E= " << m << endl;
  char cmd [256];
  sprintf(cmd, "awk '{print $1-1, $2-1}' %s >> %s.elist.t", fname2,fname1);
  
  system(cmd);
  tout.close();

  
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







void Parser_t2b(char* fname1)
{
  
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





