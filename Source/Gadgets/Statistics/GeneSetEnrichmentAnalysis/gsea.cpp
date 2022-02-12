#include "gsea.h"

bool mycmp_int (int i,int j) { return (i<j); }
bool mycmp_double (double i,double j) { return (i>j); }


/////////////////////////////////////////////////////////////////////////////////////////////////
// find whether the element j is in the vector nbv
bool find(vector<int>& nbv, int j)
{
  int n = nbv.size();
  for(int i=0; i<n; i++)
    if(j==nbv[i]) 
      return true;

  return false;
}
/////////////////////////////////////////////////////////////////////////////////////////////////


///////////////////////////////////////////////////////////////////////////
/* calculate the enrichment score (ES) 
Input:
r: r[j] : the correlation of gene j with the phenotype 
S: independently derived gene set S 

 */
void CalES(const vector<double>& r0, const vector<int>& S0)
{
  int N  = r0.size(); // total # of genes
  int M  = S0.size(); // # of genes in gene set S
    
  // sort r0
  vector<double> r(r0.begin(), r0.end());
  sort(r.begin(), r.end(), mycmp_double);

  // for(int i=0; i<N; i++) 
  //cout << r0[i] << ' ' << r[i] << endl;

  // sort S0
  vector<int> S(S0.begin(), S0.end());
  sort(S.begin(), S.end(), mycmp_int);

  //for(int j=0; j<M; j++) 
  //cout << S0[j] << ' ' << S[j] << endl;

  vector<int> Sc; // := the complementary set of S
  for(int i=0; i<N; i++) {
    if(!find(S, i))
      Sc.push_back(i);
  }

  //for(int k=0; k<N-M; k++) 
  //cout << Sc[k] << endl;

  
  // evaluate the fraction of genes in S ("hits") weighted by their 
  // correlation and the fraction of genes not in S ("misses") present
  // up to a given position i in L

  double N_R = 0;
  double p = 1.0;
  for(int j=0; j<M; j++) {
    N_R += pow(fabs(r[S[j]]), p);
  }
  cout << "N_R=" << N_R << endl;

  vector<double> Phit(N,0);
  vector<double> Pmiss(N,0);
  
  ofstream fout("test", ios::out);         

  for(int i=0; i<N; i++) {
    Phit[i] = 0;
    for(int j=0; j<M; j++) {
      if(S[j] <= i)
	Phit[i] += pow(fabs(r[S[j]]), p);
    }
    Phit[i] /= N_R;

    Pmiss[i] = 0;
    for(int j=0; j<M; j++) {
      if(Sc[j] <= i)
	Pmiss[i] ++;
    }
    Pmiss[i] /= (double)(N-M);

    cout << i << ' ' << Phit[i]-Pmiss[i] << endl;
    fout << i << ' ' << Phit[i]-Pmiss[i] << endl;
  }

}
///////////////////////////////////////////////////////////////////////////




