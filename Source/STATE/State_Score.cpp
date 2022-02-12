#include "State.h"

using namespace std;

///////////////////////////////////////////////////////////////////////////////
void State::Get_BackgroundScore(string file)
{
 Get_BackgroundScore(1000, 2000, 1, file);
}

void State::Get_BackgroundScore(int Q, int kmax, int seed, string file)
{
  char fname [256];
  sprintf(fname,"%s-Q%d-kmax%d-seed%d.bg", file.c_str(), Q, kmax, seed);

  ifstream fin(fname, ios_base::in);
  if(fin) 
    Read_BackgroundScore(fname);
  else 
    Cal_BackgroundScore(Q, kmax, seed, file);
}
///////////////////////////////////////////////////////////////////////////////////////////////// 



///////////////////////////////////////////////////////////////////////////////////////////////// 
void State::Cal_BackgroundScore(int Q, int kmax, int seed, string file)
{
  Rand rand;
  rand.seed(seed);

  cout << "For 1 <= k <= " << kmax << ", we generate " << Q << " random subsets of size k to calculate the background module z-score distribution.\n";
   
  char fname [256];
  sprintf(fname,"%s-Q%d-kmax%d-seed%d.bg", file.c_str(), Q, kmax, seed);
  ofstream fout1(fname, ios_base::out);
  if(!fout1) {cout << "Cannot open " << fname << " for write."; exit(0);}

  MU1.clear();    
  SIGMA1.clear(); 
  MU2.clear();    
  SIGMA2.clear(); 

  MU1.resize(kmax+1,0);    
  SIGMA1.resize(kmax+1,0); 
  MU2.resize(kmax+1,0);    
  SIGMA2.resize(kmax+1,0); 
  
  for(int k=1; k<=kmax; k++)    {
    cout << "k= " << k << '\r';
    // for this particular k-value, we generate Q random subsets of size k, 
    // determine the mean and variance of ZA_1 and ZA_2
    vector<double> ZA1;
    vector<double> ZA2;
    for(int t=0;t<Q;t++) {
      vector<int> S;
      rand.randomsubset_shuf(N, k, S); 
      double sum=0;
      for(int j=0;j<k;j++)
	sum += h[S[j]];
      double ZA_1 = sum/sqrt(k);
      double ZA_2 = sqrt(k) * (sum/(double)k - havg); // correction.
      
      ZA1.push_back(ZA_1);
      ZA2.push_back(ZA_2);
    }
	
    // calculate its average and standard deviation
    double ZA1ave, ZA1var;
    Statistics(ZA1, ZA1ave, ZA1var);
    double ZA2ave, ZA2var;
    Statistics(ZA2, ZA2ave, ZA2var);
    
    MU1[k] = ZA1ave;
    SIGMA1[k] = sqrt(ZA1var);
    MU2[k] = ZA2ave;
    SIGMA2[k] = sqrt(ZA2var);

    fout1 << k << ' ' << MU1[k] << ' ' << SIGMA1[k] << ' ' << MU2[k] << ' ' << SIGMA2[k] << endl;
  }
  fout1.close();

}
///////////////////////////////////////////////////////////////////////////////




///////////////////////////////////////////////////////////////////////////////
void State::Read_BackgroundScore(string file)
{
  cout << "Read the background module z-score distribution.\n";
  char fname [256];
  sprintf(fname,"%s", file.c_str());
  ifstream fin(fname, ios_base::out);
  if(!fin) {cout << "Cannot open " << fname << " for read."; exit(0);}

  MU1.clear();    
  SIGMA1.clear(); 
  MU2.clear();    
  SIGMA2.clear(); 

  MU1.push_back(0);    
  SIGMA1.push_back(0);
  MU2.push_back(0);
  SIGMA2.push_back(0);
  
  int k;
  double mu1, sigma1, mu2, sigma2;
  
  while(fin >> k >> mu1 >> sigma1 >> mu2 >> sigma2) {
    MU1.push_back(mu1);
    SIGMA1.push_back(sigma1);
    MU2.push_back(mu2);
    SIGMA2.push_back(sigma2);
    //cout << k << ' ' << MU1[k] << ' ' << SIGMA1[k] << ' ' << MU2[k] << ' ' << SIGMA2[k] << endl;
  }
  fin.close();

}
///////////////////////////////////////////////////////////////////////////////





///////////////////////////////////////////////////////////////////////////////
void State::DFS(int u, vector<bool>& visited, Nbl& Component, char STATUS)
{
  visited[u] = true;
  Component.push_back(u);
  for(Nbl_itr p = A[u].begin(); p!= A[u].end(); p++) {
    int v = (*p);
    if(!visited[v] && spin[v]==STATUS)
      DFS(v, visited, Component, STATUS);
  }    

}
///////////////////////////////////////////////////////////////////////////////




///////////////////////////////////////////////////////////////////////////////
// calculate the (normalized) scores, size, etc of the module 
double State::CalModuleScore(Nbl& NL) 
{
  double ZA1, ZA2, SA1, SA2;
  int k = NL.size();
  int kmax = MU1.size();

  // Note that here if k is so large that MU[k] and SIGMA[k] were not calculated at all, then we simply ignore those modules. 
  if(k>kmax) {
    //cout << k << " > kmax=" << kmax << ". We will set the module score to be 0!\n";  
    ZA1 = ZA2 = SA1 = SA2 = 0;
  }
  else {
    double sum =0;
    for(Nbl_itr p = NL.begin(); p!= NL.end(); p++) 
      sum += h[*p];
    
    ZA1 = sum/sqrt((double)k);              // biased
    ZA2 = sqrt(k) * (sum/(double)k - havg); // corrected
    SA1 = (ZA1 - MU1[k])/SIGMA1[k];
    SA2 = (ZA2 - MU2[k])/SIGMA2[k];
  }
  //cout << k << ' ' << ZA1 << ' ' << ZA2 << ' ' << SA1 << ' ' << SA2 << endl; //debug

  return SA2;
}
///////////////////////////////////////////////////////////////////////////////



///////////////////////////////////////////////////////////////////////////////
// calculate the (normalized) scores of the module 
void State::CalModuleScore(Nbl& NL, double& SA1, double& SA2) 
{
  double sum =0;
  for(Nbl_itr p = NL.begin(); p!= NL.end(); p++) {
    sum += h[*p];
  }
  int k = NL.size();

  // Note that here if k is so large that MU[k] and SIGMA[k] were not calculated at all, then we will have inf error.
  int kmax = MU1.size();
  if(k>=kmax) {
    cout << k << " >=  kmax=" << kmax << endl;
    exit(0);
  }

  double ZA1 = sum/sqrt((double)k);                      // biased
  double ZA2 = sqrt(k) * (sum/(double)k - havg); // corrected

  // calculate the normalized ZA values
  SA1 = (ZA1 - MU1[k])/SIGMA1[k];
  SA2 = (ZA2 - MU2[k])/SIGMA2[k];
}
///////////////////////////////////////////////////////////////////////////////



///////////////////////////////////////////////////////////////////////////////
// calculate the (normalized) scores, size, etc of the module 
double State::CalModuleScore(Nbl& NL, module& Module) 
{
  double ZA1, ZA2, SA1, SA2;
  int k = NL.size();
  int kmax = MU1.size();

  // Note that here if k is so large that MU[k] and SIGMA[k] were not calculated at all, then we simply ignore those modules. 
  if(k>kmax) {
    //cout << k << " > kmax=" << kmax << ". We will set the module score to be 0!\n";  
    ZA1 = ZA2 = SA1 = SA2 = 0;
  }
  else {
    double sum =0;
    for(Nbl_itr p = NL.begin(); p!= NL.end(); p++) 
      sum += h[*p];
    
    ZA1 = sum/sqrt((double)k);              // biased
    ZA2 = sqrt(k) * (sum/(double)k - havg); // corrected
    SA1 = (ZA1 - MU1[k])/SIGMA1[k];
    SA2 = (ZA2 - MU2[k])/SIGMA2[k];
  }

  Module.score0 = SA1; 
  Module.score  = SA2; 
  Module.size = k;
  //cout << k << ' ' << ZA1 << ' ' << ZA2 << ' ' << SA1 << ' ' << SA2 << endl; //debug

  return SA2;
}
///////////////////////////////////////////////////////////////////////////////



///////////////////////////////////////////////////////////////////////////////
double State::GetScore()
{
  vector<bool> visited(N, false);
  int count=0;

  //Modules.clear();

  double Smax = 0;
  int    Imax = 0; //  index of the module with highest score 
  for(int i=0; i<N; i++) {
    if(!visited[i] && spin[i]==ACTIVE) {
      Nbl NL;
      DFS(i, visited, NL, ACTIVE);
      //NL.sort(); // sort the nodes
      //module Module;
      //double S = CalModuleScore(NL, Module);

      double S = CalModuleScore(NL);
      if(S>Smax) {
	Smax = S;
	Imax = count;
      }

      //Module.index =  count;
      //Module.nodelist = NL;
      //Modules.push_back(Module);
      
      count++;
    }
  }
  
  //cout << "Smax =" << Smax << endl;
  return Smax;
}
///////////////////////////////////////////////////////////////////////////////





///////////////////////////////////////////////////////////////////////////////
int State::CalModule(module& M)
{
  M.index= -1;
  M.score= 0;

  vector<bool> visited(N, false);
  int count=0;
  Modules.clear();

  moduleindex.clear();
  moduleindex.resize(N,-1);

  double Scoremax = -1e300;
  int    Imax = -1; // index of the module with highest score 

  int    Sizemax = 0;
  int    Jmax = -1; // index of the module with largest size, i.e., the LCC index 

  //vector<int> modulesize; // to plot the size distribution of active components

  for(int i=0; i<N; i++) {
    if(!visited[i] && spin[i]==ACTIVE) {
      Nbl NL;
      DFS(i, visited, NL, ACTIVE);
      NL.sort(); // sort the nodes
      module Module;

      double Score = CalModuleScore(NL, Module);
      if(Score>Scoremax) {
	Scoremax = Score;
	Imax = count;
      }

      int Size = NL.size();
      if(Size>Sizemax) {
	Sizemax = Size;
	Jmax = count;
      }

      //modulesize.push_back(Size);

      Module.index =  count;
      Module.nodelist = NL;
      Modules.push_back(Module);
      
      for(Nbl_itr p = NL.begin(); p!= NL.end(); p++) 
	moduleindex[*p] = count; 
      
      count++;
    }
  }

  if(Imax >= 0) 
    M = Modules[Imax];

  if(Jmax >= 0) {
    LCCindex = Jmax;
    cout << "LCCindex = " << LCCindex << endl; 
  }

  /*
  // calculate the histogram of the size distribution of the active components
  if(modulesize.size()>0) {
    char label[256];
    sprintf(label,"H=%e; m=%e; LCC=%d", Hext, m, Sizemax);
    GetHistogram_Accumulatively(modulesize, 1.2, label, "./data/DistSize.dat");
  }
  */

  return Sizemax;
}
///////////////////////////////////////////////////////////////////////////////



///////////////////////////////////////////////////////////////////////////////
// here we only care about the LCC of the active spins, without considering its score. YYL 04/11/2016
int State::CalLCC(module& M)
{
  M.index= -1;
  M.score= 0;

  vector<bool> visited(N, false);
  int count=0;
  Modules.clear();

  moduleindex.clear();
  moduleindex.resize(N,-1);

  int    Sizemax = 0;
  int    Jmax = -1; // index of the module with largest size, i.e., the LCC index 

  //vector<int> modulesize; // to plot the size distribution of active components

  for(int i=0; i<N; i++) {
    if(!visited[i] && spin[i]==ACTIVE) {
      Nbl NL;
      DFS(i, visited, NL, ACTIVE);
      NL.sort(); // sort the nodes
      module Module;
      int Size = NL.size();
      if(Size>Sizemax) {
	Sizemax = Size;
	Jmax = count;
      }

      //modulesize.push_back(Size);

      Module.index =  count;
      Module.nodelist = NL;
      Modules.push_back(Module);
      
      for(Nbl_itr p = NL.begin(); p!= NL.end(); p++) 
	moduleindex[*p] = count; 
      
      count++;
    }
  }

  
  if(Jmax >= 0) {
    LCCindex = Jmax;
    //cout << "LCCindex = " << LCCindex << endl; 
  }

  /*
  // calculate the histogram of the size distribution of the active components
  if(modulesize.size()>0) {
    char label[256];
    sprintf(label,"H=%e; m=%e; LCC=%d", Hext, m, Sizemax);
    GetHistogram_Accumulatively(modulesize, 1.2, label, "./data/DistSize.dat");
  }
  */

  return Sizemax;
}
///////////////////////////////////////////////////////////////////////////////




///////////////////////////////////////////////////////////////
void State::SaveModules(char* file)
{
 char fname [256];
  sprintf(fname,"%s.module_info", file);
  ofstream fout(fname, ios_base::out);
  if(!fout) {cout << "Cannot open " << fname << " for write."; exit(0);}
 
  //sprintf(fname,"%s.module_nodelist", file);
  //ofstream fout2(fname, ios_base::out);
  //if(!fout2) {cout << "Cannot open " << fname << " for write."; exit(0);}


  sort(Modules.begin(), Modules.end(), ModuleSortPredicate); 

  int Nm = Modules.size(); 
  cout << Normal << "\nFinally, we have " << Nm << " modules! Here we show the top-10 modules:\n"; 
  cout << Normal << "index," << Cyan <<" score," << Purple << " score0," << Green << " size," << Red << " nodelist\n";
  for (int a=0; a<Nm && a<10; a++) {
    cout << Normal << a << ',' 
	 << Cyan << Modules[a].score << ',' 
	 << Purple << Modules[a].score0 << ',' 
	 << Green << Modules[a].size 
	 << Red << " : {";
    for(Nbl_itr p = Modules[a].nodelist.begin(); p!= Modules[a].nodelist.end(); p++) 
      //cout << *p << ',';
      cout << NodeName[*p] << ',';
    cout << "}\n" << Normal;
  }
  
  for (int a=0; a<Nm; a++) {
    fout << a << ',' << Modules[a].score << ',' << Modules[a].score0 << ',' << Modules[a].size << ',';
    for(Nbl_itr p = Modules[a].nodelist.begin(); p!= Modules[a].nodelist.end(); p++) 
      //fout2 << *p << ' '; 
      //fout2 << NodeName[*p] << ',';
      fout << NodeName[*p] << ',';
    fout << endl;

    //fout2 << "\n";
  }
  fout.close();
  //fout2.close();

}
///////////////////////////////////////////////////////////////

