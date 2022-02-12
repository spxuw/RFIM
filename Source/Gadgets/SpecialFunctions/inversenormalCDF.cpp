/*  Computing the inverse of the normal (Gaussian) CDF. That is, given a probability p, 
    compute z such that Prob(Z < z) = p where Z is a standard normal (Gaussian) random variable. 
 */

#include "inversenormalCDF.h"
using namespace std;

const double EPS = numeric_limits<double>::epsilon();

//Method-I: http://www.johndcook.com/normal_cdf_inverse.html
double RationalApproximation(double t)
{
  // Abramowitz and Stegun formula 26.2.23.
  // The absolute value of the error should be less than 4.5 e-4.
  double c[] = {2.515517, 0.802853, 0.010328};
  double d[] = {1.432788, 0.189269, 0.001308};
  return t - ((c[2]*t + c[1])*t + c[0]) / 
    (((d[2]*t + d[1])*t + d[0])*t + 1.0);
}


double NormalCDFInverse(double p)
{

/*if (p <= 0.0 || p >= 1.0) {
    stringstream os;
    os << "Invalid input argument (" << p << "); must be larger than 0 but less than 1.";
    throw  invalid_argument( os.str() );
    }*/

  // in reading the gene-pvalue file, we do find some genes have p-value=1.
  // If this happens, we just replace it with a slightly small p-value.
  if(p<=0.0)
    p = EPS;
  if(p>=1.0)
    p = 1- EPS;

  // See article above for explanation of this section.
  if (p < 0.5) // F^-1(p) = - G^-1(p)
      return -RationalApproximation( sqrt(-2.0*log(p)) );
  else   // F^-1(p) = G^-1(1-p)
    return RationalApproximation( sqrt(-2.0*log(1-p)) );
}


void demo()
{
  cout << "\nShow that the NormalCDFInverse function is accurate at \n"
       << "0.05, 0.15, 0.25, ..., 0.95 and at a few extreme values.\n\n";

  double p[] =
    {
      0.0000001,
      0.00001,
      0.001,
      0.05,
      0.15,
      0.25,
      0.35,
      0.45,
      0.55,
      0.65,
      0.75,
      0.85,
      0.95,
      0.999,
      0.99999,
      0.9999999
    };

  // Exact values computed by Mathematica.
  double exact[] =
    {
      -5.199337582187471,
      -4.264890793922602,
      -3.090232306167813,
      -1.6448536269514729,
      -1.0364333894937896,
      -0.6744897501960817,
      -0.38532046640756773,
      -0.12566134685507402,
      0.12566134685507402,
      0.38532046640756773,
      0.6744897501960817,
      1.0364333894937896,
      1.6448536269514729,
      3.090232306167813,
      4.264890793922602,
      5.199337582187471
    };

  double maxerror = 0.0;
  int numValues = sizeof(p)/sizeof(double);
  cout << "p, exact CDF inverse, computed CDF inverse, diff\n\n";
  cout <<  setprecision(7);
  for (int i = 0; i < numValues; ++i)
    {
      double computed = NormalCDFInverse(p[i]);
      double error = exact[i] - computed;
      cout << p[i] << ", " << exact[i] << ", " 
	   << computed << ", " << error << "\n";
      if (fabs(error) > maxerror)
	maxerror = fabs(error);
    }

  cout << "\nMaximum error: " << maxerror << "\n\n";
}


/* 
Method-II: http://en.wikipedia.org/wiki/Error_function
   We can relate the inverse of normal CDF with the inverse-error function. And the inverse error function can be defined in terms of the Maclaurin series:
   erf^(-1) [p] = \sum_k=0^\infty  c_k/(2k+1) (\sqrt(pi)/2 p)^{2k+1}
   with c_0 = 1
   and  c_k = \sum_{m=0}^{k-1} c_m c_{k-1-m} / [(m+1)(2m+1)] 

Comment: It is slow and not accurate at all! Why????
Perhaps this is due to the fact that the above series summation converges very slow? One has to introduce a Lekner-like efficient summation. Method-I is an example.

*/

const double SQR2 = sqrt(2.0);

double NormalCDFInverse_YYL(double p) {
  if(p<0.5)
    return -SQR2 * inverseerrorfunction(1-2*p);
  else
    return +SQR2 * inverseerrorfunction(2*p-1);
}

//   erf^(-1) [p] = \sum_k=0^\infty  c_k/(2k+1) (\sqrt(pi)/2 p)^{2k+1}
//   with c_0 = 1
//   and  c_k = \sum_{m=0}^{k-1} c_m c_{k-1-m} / [(m+1)(2m+1)] 

const double PI=3.141592653589793238462643383279;
const double SQRPI = 0.5*sqrt(PI);

double inverseerrorfunction(double p) 
{
  int N=1000;
  vector<double> c(N);
  c[0]=1;
  for(int k=1; k<N; k++) {
    c[k] =0;
    for(int m=0;m<k;m++) 
      c[k] += c[m]*c[k-1-m]/(double)((m+1)*(m+2));

    //cout << k << ' ' << c[k] << endl;
  }

  
  double p2 = 0.25*PI*p*p;
  double pk = SQRPI * p;
  double sum =0;
  double temp = 1;
  for(int k=0; k<N && temp>1e-16; k++) {
    temp = c[k]/(double)(2*k+1) * pk;
    sum += temp;
    pk *= p2;
    //cout << k << ' ' << temp << endl; 
  }
  
  return sum;

}
