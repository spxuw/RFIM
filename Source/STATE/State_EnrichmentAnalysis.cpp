#include "State.h"

using namespace std;


int State::ReadGeneSetsFile(char* genesetsfile)
{
  cout << "\n Read the benchmark gene sets file...\n";
  benchmarkgenesets.clear();

  benchmarkgeneset_name.clear();
  benchmarkgeneset_pmin.clear();
  benchmarkgeneset_pmin_H.clear();
  benchmarkgeneset_pmin_LCC.clear();
  benchmarkgeneset_pmin_H_LCC.clear();

  Color.clear();
  Color.resize(N,0);

  int count=0;

  ifstream fin(genesetsfile, ios_base::in);
  string S;
  while(getline(fin,S)) {
    
    Splitter split(S, "\t"); 

    int k = (int)split.size();
    cout << "Gene set " << split[0] << " has " << k-2 << " genes." << endl;  
   
    gs geneset;
    for(int j=2; j<k; j++) { 
      string gene = split[j];
      if(!find(NodeName, gene)) {
	cout << "Gene " << gene << " is new. Added to the node map.\n"; 
	NodeName.push_back(gene);
	NodeMAP[gene] = NodeName.size()-1;
      }

      if(NodeMAP[gene]<N) 
	  Color[NodeMAP[gene]] = 1;
      

      geneset.push_back(NodeMAP[gene]);
    }

    sort(geneset.begin(), geneset.end());
    benchmarkgenesets.push_back(geneset);

    benchmarkgeneset_name.push_back(split[0]); 
    benchmarkgeneset_pmin.push_back(1.0);
    benchmarkgeneset_pmin_H.push_back(-100);
    benchmarkgeneset_pmin_LCC.push_back(1.0);
    benchmarkgeneset_pmin_H_LCC.push_back(-100);

    count++;
  }

  cout << "In total " << count << " benchmark gene sets.\n\n";  
  fin.close();

  


  return count;
}






double State::FisherExactTest(gs& X, gs& Y)
{
  vector<int>::iterator it;

  int Ntot = NodeName.size(); 
  int Nx   = X.size();
  int Ny   = Y.size();

  vector<int> v(Nx+Ny);
  it = set_intersection (X.begin(), X.end(), Y.begin(), Y.end(), v.begin());
  v.resize(it-v.begin());

  int a = v.size(); 
  int b = Ny - a;
  int c = Nx - a;
  int d = Ntot-Nx-b;


 
  double leftpval;
  double rightpval;
  double twopval;
  fishertest (a, b, c, d, leftpval, rightpval, twopval);

  cout << "        X    non-X  " << endl;
  cout << "Y       " <<    a   << "    " << b  << " : " << Ny << endl;
  cout << "non-Y   " <<    c   << "    " << d  << " : " << Ntot-Ny << endl;
  cout << "        " <<   Nx   << "    " << Ntot-Nx  << " : " << Ntot  << endl;
  cout << "Fisher's exact test: left sided pvalue = "<< leftpval << "; right sided pvalue = "<< rightpval << "; two sided pvlaue = "<< twopval << endl << endl;

  return twopval;
  
 
}






void State::FisherExactTest(vector<double>& P)
{
  P.clear();

  gs genelist_LCC;
  gs genelist;
  for(int i=0; i<N; i++) {
    if(spin[i]==UP) {
      genelist.push_back(i);
      if(moduleindex[i]==LCCindex)
	genelist_LCC.push_back(i);
    }
  }
  
  int Nbgs = benchmarkgenesets.size(); 
  for(int I=0; I<Nbgs; I++) {
    double pvalue_LCC = FisherExactTest(genelist_LCC, benchmarkgenesets[I]); 
    double pvalue     = FisherExactTest(genelist, benchmarkgenesets[I]);     
   
    P.push_back(pvalue_LCC);
    P.push_back(pvalue);

    if(pvalue < benchmarkgeneset_pmin[I]) {
      benchmarkgeneset_pmin[I] = pvalue;
      benchmarkgeneset_pmin_H[I] = Hext;
    }

    if(pvalue_LCC < benchmarkgeneset_pmin_LCC[I]) {
      benchmarkgeneset_pmin_LCC[I] = pvalue_LCC;
      benchmarkgeneset_pmin_H_LCC[I] = Hext;
    }

  }

}












void State::WriteGeneSetsFile(char* fname)
{
  cout << "\n Read the benchmark gene sets file...\n";

  ofstream gout(fname, ios::out);
  int Nbgs = benchmarkgenesets.size(); 
  for(int I=0; I<Nbgs; I++) {
    
    gout << benchmarkgeneset_name[I] << '\t'
	 << benchmarkgeneset_pmin[I] << '\t'
	 << benchmarkgeneset_pmin_H[I] << '\t'
	 << benchmarkgeneset_pmin_LCC[I] << '\t'
	 << benchmarkgeneset_pmin_H_LCC[I];
    
    int size = benchmarkgenesets[I].size();
    for(int i=0; i<size; i++)
      gout << '\t' << NodeName[benchmarkgenesets[I][i]];
    gout << endl;
  }
}






