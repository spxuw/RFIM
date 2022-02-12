#include <iostream>
#include <fstream>
#include <vector>
#include "hypergeometric.h"

using namespace std;

///////////////////////////////////////////////////////////////////////////////////////
int main0(int argc, char** argv)
{
  if(argc!=5) {
    cout << "\nNumber of arguments is incorrect.\n";
    cout << "\nCommand format : Fisher exact test: a b c d\n ";
    exit(0);
  }
   
  int    a    = atoi(argv[1]);
  int    b    = atoi(argv[2]);
  int    c    = atoi(argv[3]);
  int    d    = atoi(argv[4]);

 pvalue_HypergeometricDistribution(a, b, c, d);
  
 
 exit(1);

}
///////////////////////////////////////////////////////////////////////////////////////



///////////////////////////////////////////////////////////////////////////////////////
int main(int argc, char** argv)
{
  if(argc!=5) {
    cout << "\nNumber of arguments is incorrect.\n";
    cout << "\nCommand format : Fisher exact test: a b c d\n ";
    exit(0);
  }
   
  int    a    = atoi(argv[1]);
  int    b    = atoi(argv[2]);
  int    c    = atoi(argv[3]);
  int    d    = atoi(argv[4]);

  double n11 = a;
  double n12 = b;
  double n21 = c;
  double n22 = d;

  double leftpval, rightpval, twopval;
 
  if( n11 < 0 || n12 < 0 || n21 < 0 || n22 < 0 ){
    cerr << "an input was negative: " << endl
	 << "n11: " << n11 << endl
	 << "n12: " << n12 << endl
	 << "n21: " << n21 << endl
          << "n22: " << n22 << endl;
    exit(1);
  }
  if(  n11 != floor(n11) || n12 != floor(n12)
       || n21 != floor(n21) || n22 != floor(n22) ){
    cerr << "an input was not an integer: " << endl
	 << "n11: " << n11 << endl
	 << "n12: " << n12 << endl
	 << "n21: " << n21 << endl
	 << "n22: " << n22 << endl;
     exit(1);
  }

  fishertest( n11, n12, n21, n22, leftpval, rightpval, twopval);
  cout   << "n11: " << n11 << endl
	 << "n12: " << n12 << endl
	 << "n21: " << n21 << endl
	 << "n22: " << n22 << endl;
  cout << "left sided test pvalue:  "<< leftpval << endl
       << "right sided test pvalue: "<< rightpval << endl
       << "two sided test pvlaue:   "<< twopval << endl;
 
 
  return 0;
}
///////////////////////////////////////////////////////////////////////////////////////




///////////////////////////////////////////////////////////////////////////////////////
int main3(int argc, char** argv)
{
  int first[] = {5,10,15,20,25};
  int second[] = {50,40,30,20,10};

  std::sort (first,first+5);     //  5 10 15 20 25
  std::sort (second,second+5);   // 10 20 30 40 50
 
  vector<int> X(5);
  vector<int> Y(5);
  copy(first, first+5, X.begin());
  copy(second, second+5, Y.begin());

  int Nx   = X.size();
  int Ny   = Y.size();
  int Ntot = Nx + Ny; 
  
  vector<int> v(Nx+Ny);
  vector<int>::iterator it = set_intersection (X.begin(), X.end(), Y.begin(), Y.end(), v.begin());
  v.resize(it-v.begin());

  int a = v.size(); // the size of the intersection of X and Y 
  int b = Ny - a;
  int c = Nx - a;
  int d = Ntot-Nx-b;

  cout << "        X    non-X  " << endl;
  cout << "Y       " <<    a   << "    " << b  << " : " << Ny << endl;
  cout << "non-Y   " <<    c   << "    " << d  << " : " << Ntot-Ny << endl;
  cout << "        " <<   Nx   << "    " << Ntot-Nx  << " : " << Ntot << endl << endl;
 
  double leftpval;
  double rightpval;
  double twopval;
  cout << "Fisher's exact test: "; //pvalue_HypergeometricDistribution(a, b, c, d);
  fishertest (a, b, c, d, leftpval, rightpval, twopval);
  cout << "left sided pvalue = "<< leftpval << "; right sided pvalue = "<< rightpval << "; two sided pvlaue = "<< twopval << endl << endl;

  return twopval;
 
}
/////////////////////////////////////////////////////////////////////////////////////////////////
