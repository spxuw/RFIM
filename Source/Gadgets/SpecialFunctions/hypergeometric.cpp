#include "hypergeometric.h"

/*  
Given a contigency table 
               Men     Women  	Row Total
Dieting	        a	b 	a + b = r
Non-dieting  	c	d 	c + d = s
Column Total	a + c	b + d	a + b + c + d (=N)
                  ||     ||   
                  m      n

Fisher showed that the probability of obtaining any such set of values was given by the hypergeometric distribution:

         (a+b, a) (c+d, c) 
p(a) = --------------------
         (N, a+c)

     = (a+b)! (c+d)! (a+c)! (b+d)! / [a! b! c! d! N!]

which can be used to calculate the p-value.
*/

// Method-I: Does not work for large numbers!! 
 
//////////////////////////////////////////////////////////////////////////////////////////////////
// return the binomial coefficient: C(m,k)
double binom(int m, int k)
{
  double c = 1.;
  for(int i=0; i<k; i++) { 
    double x= (double)(m-i)/(double)(i+1);
    c *= x;
    //cout << x << ',' << c << endl; // Note that for large m and k, very soon we will overflow problem.
  }

  //cout << "binom(" << m << ',' << k << ")=" << c << endl;  //debug

  return c;
}
//////////////////////////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////////////////////////
double Px_HypergeometricDistribution(int a, int b, int c, int d)
{
  // This is not good for large numbers.
  return binom(a+b, a) * binom(c+d, c) / binom(a+b+c+d, a+c);


/*One can easily prove the recursive relationship:????
p(a+1) = p(a) * ???
so, if we calculate p(a=0), we can then calculate all p(a) for a>0
*/

  /* The following is problematic!!!
  double P0=1.;
  for(int i=1; i<=b; i++) {
    P0 *= (double)(d+i)/(double)(c+d+i);
  }
  if(a==0) return P0;
  else {
    double Pa = P0;
    double x = b*c/(double)(d+1);
    for(int j=0; j<a; j++)
      Pa *= x/(double)(j+1); 
    return Pa;
  }
  */
}
//////////////////////////////////////////////////////////////////////////////////////////////////



//////////////////////////////////////////////////////////////////////////////////////////////////
/*
Thus, we can calculate the exact probability of any arrangement of the n teenagers into the 
four cells of the table, but Fisher showed that to generate a significance level, we need consider 
only the cases where the marginal totals are the same as in the observed table, and among those, 
only the cases where the arrangement is as extreme as the observed arrangement, or more so. 

An approach used by the fisher.test function in R is to compute the p-value by summing the probabilities 
for all tables with probabilities less than or equal to that of the observed table. 
*/

double pvalue_HypergeometricDistribution(int a, int b, int c, int d)
{
  double p0 = Px_HypergeometricDistribution(a, b, c, d);
  cout << "p0=" << p0 << endl;
  exit(0);

  double pvalue = 0;//p0;
  for(int x=0; x<=a+b; x++) 
  {
    double p = Px_HypergeometricDistribution(x, a+b-x, a+c-x, d-a+x);
    cout << x << ' ' << p << ' ' << p0 << endl;
    if(p<=p0)
      pvalue += p;
  }
  
  cout << "pvalue = " << pvalue << endl;
  return pvalue;
}
//////////////////////////////////////////////////////////////////////////////////////////////////







// Method-II: 
// The following code is taken from: http://lists.mcgill.ca/scripts/wa.exe?A2=ind9808b&L=stat-l&P=5654

//from numerical recipes in c
double gammln(double xx)
{
        double x,tmp,ser;
        static double cof[6]={76.18009173,-86.50532033,24.01409822,
                -1.231739516,0.120858003e-2,-0.536382e-5};
        int j;
 
        x=xx-1.0;
        tmp=x+5.5;
        tmp -= (x+0.5)*log(tmp);
        ser=1.0;
        for (j=0;j<=5;j++) {
                x += 1.0;
                ser += cof[j]/x;
        }
        return -tmp+log(2.50662827465*ser);
}
 
 
double logcomb ( double n, double x)
{   double lt = gammln( n+1.0 )
              - gammln(x+1.0)
              - gammln(n-x+1.0);
         return lt;
}

double hypergeometricprob(double x,
                          double n1d,
                          double n2d,
                          double nd1,
                          double nd2)
{   double n3 = nd1-x;
         double ndd = nd1 + nd2;
         double sum = logcomb( n1d, x)
              + logcomb( n2d, n3)
              - logcomb( ndd, nd1);
         return exp(sum);
}


// fisher's exact test
void fishertest (double a, double b, double c, double d, double &leftpval, double &rightpval, double &twopval)
{
  double r = a + b;
  double s = c + d;
  double m = a + c;
  double n = b + d;
  
  /* get range of variation */
  double lm = ( 0 > m-s ) ? 0 : m-s;
  double um = ( m < r   ) ? m : r;
  if( !((int)(um-lm+2)) ){
    leftpval= rightpval= twopval=1.0;
    return;
  }
  
  /* Fisher's exact test pval */
  double crit = hypergeometricprob( a, r, s, m, n );
  
  leftpval= rightpval= twopval=0;
  for( double x=lm; x<=um; x++){
    double prob = hypergeometricprob( x, r, s, m, n );
    if ( x <= a )      leftpval += prob;
    if ( x >= a )      rightpval += prob;
    if ( prob <= crit) twopval += prob;
  }
}
 
