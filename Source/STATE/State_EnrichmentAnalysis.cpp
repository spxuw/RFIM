#include "State.h"

using namespace std;


/////////////////////////////////////////////////////////////////////////////////////////////////
/* read benchmark gene sets from a multiple-row file (each row correspond to a particular benchmark gene set in the following GSEA format: 
EXTRINSIC_TO_PLASMA_MEMBRANE	http://www.broadinstitute.org/gsea/msigdb/cards/EXTRINSIC_TO_PLASMA_MEMBRANE	GNA14	APC2	GNAI1	SCUBE1	RGS19	EEA1	ARRB1	TDGF1	SYTL4	TGM3	SYTL2	GNAS	SYTL1
i.e., the first column is the geneset name, the second is its URL, the rest columns are the genes in the geneset. 

we then do the fisher exact test by considering each benchmark geneset and the disease module!
*/

/////////////////////////////////////////////////////////////////////////////////////////////////
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
    //Splitter split(S, " ");
    Splitter split(S, "\t"); // Note that the GSEA data files are using tab for separator!

    int k = (int)split.size();
    cout << "Gene set " << split[0] << " has " << k-2 << " genes." << endl;  
   
    gs geneset;
    for(int j=2; j<k; j++) { // note that we start from j=2, not j=0
      string gene = split[j];
      if(!find(NodeName, gene)) {
	cout << "Gene " << gene << " is new. Added to the node map.\n"; 
	NodeName.push_back(gene);
	NodeMAP[gene] = NodeName.size()-1;
      }

      if(NodeMAP[gene]<N) // if we find an existing gene in the network to be a benchmark gene, we color it
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
/////////////////////////////////////////////////////////////////////////////////////////////////





/////////////////////////////////////////////////////////////////////////////////////////////////
/*
                 X       non-X
         Y       a         b       :   Ny
     non-Y       c         d       : Ntot-Ny
                Nx      Ntot-Nx       Ntot
*/

double State::FisherExactTest(gs& X, gs& Y)
{
  vector<int>::iterator it;

  int Ntot = NodeName.size(); // total number of genes (including genes in PPI and new genes found in the benchmark gene sets.
  int Nx   = X.size();
  int Ny   = Y.size();

  vector<int> v(Nx+Ny);
  it = set_intersection (X.begin(), X.end(), Y.begin(), Y.end(), v.begin());
  v.resize(it-v.begin());

  int a = v.size(); // the size of the intersection of X and Y 
  int b = Ny - a;
  int c = Nx - a;
  int d = Ntot-Nx-b;

  /*
  cout << "X:";
  for (it=X.begin(); it!=X.end(); ++it)
    cout << ' ' << NodeName[*it];
  cout << '\n';

  cout << "Y:";
  for (it=Y.begin(); it!=Y.end(); ++it)
    cout << ' ' << NodeName[*it];
  cout << '\n';
  
  cout << "The intersection has " << (v.size()) << " elements:\n";
  for (it=v.begin(); it!=v.end(); ++it)
    cout << ' ' << NodeName[*it];
  cout << '\n';
  */
 
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
  //return rightpval;
 
}
/////////////////////////////////////////////////////////////////////////////////////////////////




/////////////////////////////////////////////////////////////////////////////////////////////////
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
  
  int Nbgs = benchmarkgenesets.size(); // # of benchmark gene sets. 
  for(int I=0; I<Nbgs; I++) {
    double pvalue_LCC = FisherExactTest(genelist_LCC, benchmarkgenesets[I]); // do the Fisher test for the LCC only 
    double pvalue     = FisherExactTest(genelist, benchmarkgenesets[I]);     // do the Fisher test for all the active genes
   
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
/////////////////////////////////////////////////////////////////////////////////////////////////










/////////////////////////////////////////////////////////////////////////////////////////////////
void State::WriteGeneSetsFile(char* fname)
{
  cout << "\n Read the benchmark gene sets file...\n";

  ofstream gout(fname, ios::out);
  int Nbgs = benchmarkgenesets.size(); // # of benchmark gene sets. 
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
/////////////////////////////////////////////////////////////////////////////////////////////////





