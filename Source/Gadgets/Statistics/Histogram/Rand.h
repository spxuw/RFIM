#if !defined (RAND_H)
#define RAND_H

#include <cstdlib>
#include <iostream>
#include <cmath>
#include <vector>
#include <list>
#include <set>
#include <limits>
#include <gsl/gsl_sf_zeta.h>

//#include "mynrutil.h"
//#include "nr.h"

using namespace std;

const double TOLERANCE = numeric_limits<double>::epsilon();
const int    INTMAX    = numeric_limits<int>::max();

const int BIGPRIME = 899999963;
const double PI    = 3.1415926535897932384626433832795028841971693993751;

typedef double (*FUNC)(const double, const double, const double);
inline double fx(const double htilde, const double y, const double u) // here htilde=x is the argument of the function
{
  return (pow(htilde,y+1) - (y+1)*htilde + fabs(2*u-1)*y);
}    


inline int Getnode(list<int>& X, int index)
{
  int pos=0;
  for(list<int>::iterator p = X.begin(); p!= X.end(); p++, pos++) {
    if(pos==index) 
      return (*p);
  }

  return 0;
}

/*
inline double gammln(const double xx)
{
  int j;
  double x,y,tmp,ser;
  static const double cof[6]={76.18009172947146,-86.50532032941677,
			      24.01409824083091,-1.231739572450155,0.1208650973866179e-2,
			      -0.5395239384953e-5};

  y=x=xx;
  tmp=x+5.5;
  tmp -= (x+0.5)*log(tmp);
  ser=1.000000000190015;
  for (j=0;j<6;j++) ser += cof[j]/++y;
  return -tmp+log(2.5066282746310005*ser/x);
}
*/

class Rand {
  double c,cd,cm,u[97];
  int i97,j97 ;
  bool outputReady;
  double output;

  int idum; // used in ran1()

 public:
  Rand(int seed) {
    ranmarin(fabs(seed%BIGPRIME));
    outputReady=false;
    idum=-1*seed; // very important , cannot just set to be -1, must depend on seed

    flag_setFvector = false;
    flag_setPvector = false;
  }

  Rand() {
    ranmarin(1);
    outputReady=false;
    idum=-1;

    flag_setFvector = false;
    flag_setPvector = false;
  }

  void seed(int seed) {
    ranmarin(fabs(seed%BIGPRIME));
    outputReady=false;
    idum=-1*seed; // very important  , cannot just set to be -1, must depend on seed

    flag_setFvector = false;
    flag_setPvector = false;
  }


  // the worst uniform RNG
  double ran() { return rand()/(double)RAND_MAX;}
    
  double uniform()
  {
    double uni;
	
    uni=u[i97]-u[j97];
    if (uni<0.0) uni+=1.0;
    u[i97]=uni;
    if (--i97<0) i97=96;
    if (--j97<0) j97=96;
    c-=cd;
    if (c<0.0) c+=cm;
    uni-=c;
    if (uni<0.0) uni+=1.0;
    return(uni);
  }

  // "Minimal" random number generator of Park and Miller with Bays-Durham shuffle and added safeguards.
  // adopted from NR. YL 10/05/08
  double ran1()
  {
    const int IA=16807,IM=2147483647,IQ=127773,IR=2836,NTAB=32;
    const int NDIV=(1+(IM-1)/NTAB);
    const double EPS=3.0e-16,AM=1.0/IM,RNMX=(1.0-EPS);
    static int iy=0;
    static vector<int> iv(NTAB);
    int j,k;
    double temp;

    if (idum <= 0 || !iy) {
      if (-idum < 1) idum=1;
      else idum = -idum;
      for (j=NTAB+7;j>=0;j--) {
	k=idum/IQ;
	idum=IA*(idum-k*IQ)-IR*k;
	if (idum < 0) idum += IM;
	if (j < NTAB) iv[j] = idum;
      }
      iy=iv[0];
    }
    k=idum/IQ;
    idum=IA*(idum-k*IQ)-IR*k;
    if (idum < 0) idum += IM;
    j=iy/NDIV;
    iy=iv[j];
    iv[j] = idum;
    if ((temp=AM*iy) > RNMX) return RNMX;
    else return temp;
  }



  // Gives all numbers weighted equally, min and max inclusive
  int discrete(int min, int max)
  {
    //int ret= (int)(floor(uniform()*(max+1-min)+min));
    int ret= (int)(floor(ran1()*(max+1-min)+min));
    //		assert(ret<=max && ret>=min);
    return ret;
  }


  // Gives all numbers weighted equally, min and max inclusive
  // generate two different random intergers
  void discrete(int min, int max, int& a, int& b)
  {
    a = (int)(floor(ran1()*(max+1-min)+min));
    
    do {
      b = (int)(floor(ran1()*(max+1-min)+min));
    } while (a==b);
    
  }

  // Gives all numbers weighted equally, min and max inclusive
  // generate three different random intergers
  void discrete(int min, int max, int& a, int& b, int& c)
  {
    a = (int)(floor(ran1()*(max+1-min)+min));
    
    do {
      b = (int)(floor(ran1()*(max+1-min)+min));
    } while (b==a);

    do {
      c = (int)(floor(ran1()*(max+1-min)+min));
    } while (c==b || c==a);
    
  }

  // Gives all numbers weighted equally, min and max inclusive
  // generate 4 different random intergers
  void discrete(int min, int max, int& a, int& b, int& c, int& d)
  {
    a = (int)(floor(ran1()*(max+1-min)+min));
    
    do {
      b = (int)(floor(ran1()*(max+1-min)+min));
    } while (b==a);

    do {
      c = (int)(floor(ran1()*(max+1-min)+min));
    } while (c==b || c==a);

    do {
      d = (int)(floor(ran1()*(max+1-min)+min));
    } while (d==c || d==b || d==a);
   
  }




  // Gives all numbers weighted equally, min and max inclusive
  // generate n different random intergers
  void discrete(int min, int max, int n, set<int>& Set)
  {
    Set.clear();
    for(int i=0; i<n; ) {
      int a = (int)(floor(ran1()*(max+1-min)+min));
      if(Set.find(a)==Set.end()) {
	Set.insert(a);
	i++;
      }
    }
  }





  // min and max inclusive, gives all numbers in log-bin equally : [1,2), [2,4), [4,8), .... 
  int discrete_logbin(int xmin, int xmax)
  {
    double ratio = 1.15;
    if (xmin==0)
      xmin=1;

    int Nbins = (int)(floor(log(xmax/xmin)/log(ratio))) + 1; 
    // choose a bin uniformly
    int i = discrete(0, Nbins-1);
    int a = xmin * pow(ratio, (double)i);
    int b = a*ratio;
    if(b>xmax)
      b = xmax;
    
    int x = discrete(a,b);
    return x;
  }



  ///////////////////////////////////////////////////////////////////////////////
  void randomsubset(int N, int Nc, vector<int>& S)
  {
    list<int> Lbig;
    for(int i=0; i<N; i++) {
      Lbig.push_back(i);
    }

    int N0= N;
    int M = (2*Nc>N0) ? (N0-Nc) : Nc;
    
    list<int> Lsmall;
    for(int i=0; i<M; i++) {
      int index = discrete(0, N-1);
      //cout << index << endl; //test
      int j = Getnode(Lbig, index);
      Lsmall.push_back(j);
      Lbig.remove(j);
      N--;
    }

    S.clear();
    S.resize(Nc);

    if(2*Nc>N0) 
      copy(Lbig.begin(), Lbig.end(), S.begin());
    else
      copy(Lsmall.begin(), Lsmall.end(), S.begin());

  }
  ///////////////////////////////////////////////////////////////////////////////



  ///////////////////////////////////////////////////////////////////////////////
  // http://stackoverflow.com/questions/136474/best-way-to-pick-a-random-subset-from-a-collection/2564196#2564196
  // This is much faster!
  void randomsubset_knuth(int n, int m, vector<int>& S)
  {    
    //cout << m << ':';
    for (int i = 0; i < n; i++)
      /* select m of remaining n-i */
      if ((discrete(1,n) % (n-i)) < m) {
	//cout << i << ",";
	S.push_back(i);
	m--;
      }
    //cout << endl;
  }
  ///////////////////////////////////////////////////////////////////////////////

 
  /*
  // Poisson distribution
  int poisson(double xm)
  {
    static double sq,alxm,g,oldm=(-1.0);
    double em,t,y;

    if (xm < 12.0) {
      if (xm != oldm) {
	oldm=xm;
	g=exp(-xm);
      }
      em = -1;
      t=1.0;
      do {
	++em;
	t *= ran1();
      } while (t > g);
    } else {
      if (xm != oldm) {
	oldm=xm;
	sq=sqrt(2.0*xm);
	alxm=log(xm);
	g=xm*alxm-gammln(xm+1.0);
      }
      do {
	do {
	  y=tan(PI*ran1());
	  em=sq*y+xm;
	} while (em < 0.0);
	em=floor(em);
	t=0.9*(1.0+y*y)*exp(em*alxm-gammln(em+1.0)-g);
      } while (ran1() > t);
    }
    return (int)em;
  }
  */


  double gaussian(double sd)
  {
    double v1,v2,Rsq,z;
    if(outputReady)
      {
	outputReady = false;
	return output;
      }

    do
      {
	v1=(ran1()*2.0)-1.0; //v1=(uniform()*2.0)-1.0;
	v2=(ran1()*2.0)-1.0; //v2=(uniform()*2.0)-1.0;
	Rsq=v1*v1+v2*v2;
      } while (Rsq>=1.0);

    z=sqrt(-2.0*log(Rsq)/Rsq);
    output = v1*z*sd;
    outputReady=true;
    return v2*z*sd;
  }

  /* John Carpenter added */
  double lorentzian(double sd) 
  {
    double v1,v2,z;
    do {
      v2 = uniform();
      v1 = 2.0 * uniform() - 1.0;
      //v2 = ran1();            
      //v1 = 2.0 * ran1() - 1.0;
    } while ((v1*v1+v2*v2)>1.0);
    z = v2/v1;
    return sd*z/2.0; // in this way, sd=w=2 gamga = width
    //return sd*z;//Yang Liu, in this way, sd = gama, i.e. the half-width at half-maximum
  }


  /* Yang Liu added */
  //////////////////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////   power-law distribution //////////////////////////////////////////////
  // continuous distribution
  double powerlaw_continuous(double alpha, double xmin)
  {
    return ( xmin*pow(1-ran1(), 1./(1-alpha) ) );
  }

  // continuous powerlaw distribution within an interval (a, 1]
  // p(theta) \sim theta^{-beta}  for  theta \in (a, 1]
  //               0                    otherwise
  // here a = alpha/n
  double powerlaw_continuous_interval(double a, double beta)
  {
    return pow( (1- (1-ran1())*(1-pow(a, 1-beta))), 1/(1-beta) );
  }



  // discrete distribution approximated from the continuous disitribution
  // The approximation involved in this approach is largest
  // for the smallest value of x, which is by definition xmin.
  int powerlaw_discrete_approximate(double alpha, double xmin)
  {
    double r = ran1();
    int k = (int)( (xmin-0.5)*pow(1-r, 1./(1-alpha) ) + 0.5 );
    if(k<0) {
      cout << 1-r << ',' << k << endl;
      k =  INTMAX;
    }
    return k;

    //return (int)( (xmin-0.5)*pow(1-ran1(), 1./(1-alpha) ) + 0.5 );
  }


  // discrete distribution obtained by ``doubling up" and binary search
  int powerlaw_discrete(double alpha, int xmin)
  {
    // Set up the look-up-table for P(x)=Pr(X>=x)
    if(!flag_setPvector)
      SetPvector(alpha,xmin);

    double r = ran1();

    /*
    // For a given random number r, we first bracket a solution x to the equation
    // Pr(X>=x) = zeta(alpha,x)/zeta(alpha,xmin) = 1-r
    int x1,x2;
    x2 = xmin;
    do {
      x1 = x2;
      x2 = 2*x1;
    } while(P[x2]>1-r); // this will cause trouble when r=1. endless loop!
    // Actually this ``doubling up" is not necessary at all. Since we will do the binary search
    // any way.
    */


    // Then we pinpoint the solution within the range x1 to x2
    // by binary search. We need only continue the binary
    // search until the value of x is narrowed down to j ≤ x <
    // j+1 for some integer j: then we discard the integer part
    // and the result is a power-law distributed random integer.
    int j = Findj(P, 1-r);
    return j;

  }
  //////////////////////////////////////////////////////////////////////////////////////////////////





  // discrete distribution obtained by ``doubling up" and binary search
  int powerlaw_exponentialcutoff_discrete(double gamma, double kappa)
  {
    // Set up the look-up-table for P(x)=Pr(X>=x)
    if(!flag_setPvector)
      SetPvector_SFEC(gamma, kappa);

    double r = ran1();


    // Then we pinpoint the solution within the range x1 to x2
    // by binary search. We need only continue the binary
    // search until the value of x is narrowed down to j ≤ x <
    // j+1 for some integer j: then we discard the integer part
    // and the result is a power-law distributed random integer.
    int j = Findj(P, 1-r);
    return j;

  }
  //////////////////////////////////////////////////////////////////////////////////////////////////




  //////////////////////////////////////////////////////////////////////////////////////////////////

  ////////////////////     For arbitrary discrete distribution with given cumulative distribution function (CDF) /////
  //int arbitrary_discrete(vector<double>& CDF, int kmin)
  //{
  //double r = ran1();
  //int j = Findj(CDF, r);
  //return j+1+kmin;
  //}
  int arbitrary_discrete(vector<double>& CDF)
  {
    double r = ran1();
    // Given a monotonic array xx[0...n-1] and given a value x, returns a value j such that 
    // xx[j] < x < xx[j+1]. 
    int j = Findj(CDF, r);
    return j+1;
  }

  int arbitrary_discrete(vector<long double>& CDF)
  {
    double r = ran1();
    // Given a monotonic array xx[0...n-1] and given a value x, returns a value j such that 
    // xx[j] < x < xx[j+1]. 
    int j = Findj(CDF, r);
    return j+1;
  }

  //////////////////////////////////////////////////////////////////////////////////////////////////


  /* Yang Liu added */
  double lognormal(double sd) 
  {
    double v1,v2,Rsq,z;
    if(outputReady)
      {
	outputReady = false;
	return output;
      }

    do
      {
	v1=(ran1()*2.0)-1.0; //v1=(uniform()*2.0)-1.0;
	v2=(ran1()*2.0)-1.0; //v2=(uniform()*2.0)-1.0;
	Rsq=v1*v1+v2*v2;
      } while (Rsq>=1.0);

    z=sqrt(-2.0*log(Rsq)/Rsq);
    output = exp(v1*z*sd);
    outputReady=true;
    return exp(v2*z*sd);
  }


  /* Yang Liu added */
  double bimodal(double R)
  {
    //return (uniform()>0.5) ? (R) : (-R);
    return (ran1()>0.5) ? (R) : (-R);
  }

  /* Yang Liu added */
  double rectangular(double R)
  {
    //return ((uniform()*2-1)*R); // Bad ! Period is short!
    return ((ran1()*2-1)*R); 
  }


  /* Yang Liu added */
  // This parabolic distribution is obtained by using the transformation method. See NR:Ch7.2, Pg.291.
  // normalized parabolic distribution:  P(h) = (R^2-h^2)/(4R^3/3) with h \in [-R,R]
  // F(h) = \int_(-R)^h p(h) dh = 1/2 - h^3/(4R^3) + 3h/(4R)
  // given x \in [0,1), we want to find the root(s) of F(h) = x, i.e. 
  //              h^3 - 3R^2 h + 4 R^3(x-1/2) == 0

  // Basically we are sovling a cubic equation.
  // Note that here the Discriminant (Disc) for this simple cubic equation is negative (crucially due to x\in[0,1))
  // which means there are three real roots.
  // Now it is tricky to apply the analytic formula, e.g see http://en.wikipedia.org/wiki/Cubic_equation
  // Because Discriminant is negative, we have to explicitly write the three real roots rather than 
  // directly applying those formulas. Note that the C++ cannot handle sqrt(negative number). It will give you "nan"
  double parabolic(double R)
  {
    //double x = uniform();
    double x = ran1();

    int sign =  (x>0.5) ? (-1) : (+1);

    double R2 = R*R;
    double R3 = R2*R;

    double q     = 4*R3*(x-0.5);
    double p     = -3*R2;
    double Disc  = q*q/4. + p*p*p/27.;
    double r     = R3;//sqrt(q*q/4. - Disc);
    double theta = atan(sqrt(-Disc)/(-q/2.));
    double ronethird = pow(r,1/3.);
	  
    double h1 = sign * 2 * ronethird * cos(theta/3.);
    double h2 = sign * 2 * ronethird * cos((theta+2*PI)/3.);
    double h3 = sign * 2 * ronethird * cos((theta-2*PI)/3.);

    double y = 0;
    if(h1<R && h1>(-R)) 
      y = h1;
    else if(h2<R && h2>(-R)) 
      y = h2;
    else if(h3<R && h3>(-R)) 
      y = h3;
    else
      cout << "range is wrong!";
	
    //cout << x << ' ' << h1 << ' ' << h2 << ' ' << h3 << ' ' << y << endl;

    return y;
  }



  // This anti-parabolic distribution is obtained by using the transformation method. See NR:Ch7.2, Pg.291.
  // normalized parabolic distribution:  P(h) = h^2/(2R^3/3) with h \in [-R,R]
  // F(h) = \int_(-R)^h p(h) dh = h^3/(2R^3) + 1/2
  // given x \in [0,1), we want to find the root(s) of F(h) = x, i.e. 
  //              h^3 = (2x-1) R^3
  // This is easy: if x>0.5, then h =   (2x-1)^(1/3) * R
  //                    else      h = - (1-2x)^(1/3) * R

  double antiparabolic(double R)
  {
    //double x = uniform();
    double x = ran1();

    int sign =  (x>0.5) ? (+1) : (-1);

    double h = sign * pow(sign*(2*x-1),1./3.) * R;
    return h;
  }




  // Duxbury-Meinke type-I distribution:
  //
  //          y+1        
  // P(h) =  ----- [1 - (|h|/R)^y]    for -R <= h <= R
  //          2yR       
  // 
  // This can also be obtained by using the transformation method. See NR:Ch7.2, Pg.291.

  // Method 1: we can directly use a look-up-table. The idea follows:
  // We integrate 
  // F(h) = \int_(-R)^h p(h) dh  for a set of h values: h = -R, -R+dh, -R+2dh, ......, +R
  // and store those F(h) values in an array.
  // Then for given x \in [0,1), we just look up the F(h)-vector, the index of the array element with F(h)=x will
  // give the random variable h.
  double DMI_lookuptable(double y, double R)
  //double DMI(double y, double R)
  {
    if(!flag_setFvector)
      SetFvector(y,R);

    double u = ran1();
    int j = Findj(F, u);
    double h = -R + j*dh;

    return h;
  }


  // Method 2: we can also numerically solve the equation: F(h) == x to get h for a given x
  //double DMI_findroot(double y, double R)
  double DMI(double y, double R)
  {
    double u = ran1();
    int sign =  (u>0.5) ? (+1) : (-1);

    // Depending on u>0.5 or not, there are two cases (h>0 or h<0).
    // we need to find the root of the equation F(h) == u for a given u. 
    // It can be derived as: ~h^(y+1) - (y+1)*~h + sign* (2y*u - y) == 0 
    // with ~h := sign*h/R. 

    const double hacc = 1e-12;
    double htilde = rtbis(&fx, y, u, 0, 1, hacc);
    double h = sign*htilde*R;
    return h;
  }
    


  ////////// for discrete exponential distribution 
  // P(k) = [1-e^(-1/kappa)] e^(-k/kappa)
  // one can solve for the cumulative distribution function exactly:
  // CDF(K) = sum_{0}^K P(k)  = 1-e^(-(K+1)/kappa) = u with u in [0,1]
  // Inverse the CDF(K)=u we get
  // K(u) = -kappa log(1-u) - 1
  int exponential(double kappa)
  {
    double u = ran1();
    //int k = -kappa * log(1-u) -1;    // (1) this is not very good. But why?
    int k = -kappa * log(1-u);  // (2) this is better in the sense that the calculated nD match the analytical value
    return k;
  }
  
  /*
    (1) is obtained by inversing the CDF(K) = sum_{k=0}^K p(k)
    (2) is obtained by inversing the CCDF(K)= sum_{k=K}^\infty p(k), i.e. the
    complimentary cumulative distribution.
    See Clasuset's paper "Power-law distribution in empirical data", p39. Eq.D5.
    
    But why (2) is better?
  */



 private:

  /////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////////////


  ////////////////////     For DMI, method 1   ////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////////////
  // For DMI-distribution, we use the look-up-table.
  // Integrate F(h) = \int_(-R)^h p(h) dh  for a set of h values: h = -R, -R+dh, -R+2dh, ......, +R
  // and store those F(h) values in Fvector.

  bool flag_setFvector;
  vector<double> F;
  double dh;

  double P_DMI(double h, double y, double R)
  {
    if(h<-R || h>R) return 0;
    else            return (y+1)/(2*y*R)*(1-pow(fabs(h)/R,y));
  }

  void SetFvector(double y, double R)
  {
    int N = 10000001; 
    dh = R/(double)(N-1);
    F.resize(2*N-1,0);
	    
    //cout << dh << endl; // debug

    F[0]     = 0;
    F[2*N-2] = 1 - F[0];
    for(int i=1;i<N;i++)
      {
	double h1 = -R + i*dh;
	double h2 = h1 + dh;
		
	// use the trapezoidal rule to do the integration
	F[i]       = F[i-1] + 0.5* (P_DMI(h1,y,R)+P_DMI(h2,y,R))*dh;
	F[2*N-2-i] = 1- F[i];
      }
	    
    flag_setFvector = true;
  }

  // use bisection-method to find the corresponding h with F[h]=x
  // based on NR:locate(Vec_I_DP& xx, const DP x, int& j)
  // Given an array xx[0...n-1] and given a value x, returns a value j such that x is between
  // xx[j] and xx[j+1]. xx must be monotonic, either increasing or decreasing.
  // j=-1 or j=n-1 is returned to indicate that x is out of range
  int Findj(vector<double>& xx, const double x)
  {
    int j,jl,jm,ju;
    int n=xx.size();

    jl=-1; ju=n;  // initialize lower and upper limits

    bool ascnd = (xx[n-1] >= xx[0]);   // true if ascending order of table, false otherwise
    while (ju-jl > 1) 
      {
	jm=(ju+jl) >> 1;  // compute a midpoint
	if ((x >= xx[jm]) == ascnd) jl=jm; // replace either the lower limit
	else    	          ju=jm; // or the upper limit, as appropriate       
      }
	    
    if (x == xx[0]) j=0;
    else if (x == xx[n-1]) j=n-2;
    else j=jl;

    return j;
  }

  int Findj(vector<long double>& xx, const double x)
  {
    int j,jl,jm,ju;
    int n=xx.size();

    jl=-1; ju=n;  // initialize lower and upper limits

    bool ascnd = (xx[n-1] >= xx[0]);   // true if ascending order of table, false otherwise
    while (ju-jl > 1) 
      {
	jm=(ju+jl) >> 1;  // compute a midpoint
	if ((x >= xx[jm]) == ascnd) jl=jm; // replace either the lower limit
	else    	          ju=jm; // or the upper limit, as appropriate       
      }
	    
    if (x == xx[0]) j=0;
    else if (x == xx[n-1]) j=n-2;
    else j=jl;

    return j;
  }

  /////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////////////


  /////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////     For DMI, method 2   ////////////////////////////////////////////////////////


  // Using bisection, find the root of a function func known to lie between x1 and x2. The root, returned as
  // rtb, will be refined until its accuracy is +/- xacc
  double rtbis(FUNC func, const double y, const double u, const double x1, const double x2, const double xacc)
  {
    const int JMAX=40; // maximum allowed number of bisections.
    int j;
    double dx,f,fmid,xmid,rtb;
	    
    f=func(x1,y,u);
    fmid=func(x2,y,u);
    if (f*fmid >= 0.0) 
      cout << "Root must be bracketed for bisection in rtbis\n";
	    
    rtb = f < 0.0 ? (dx=x2-x1,x1) : (dx=x1-x2,x2);
    for (j=0;j<JMAX;j++) {
      fmid=func(xmid=rtb+(dx *= 0.5),y,u);
      if (fmid <= 0.0) rtb=xmid;
      if (fabs(dx) < xacc || fmid == 0.0) return rtb;
    }
    cout << "Too many bisections in rtbis\n";
    return 0.0;
  }

  /////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////////////


    


  void ranmarin(int ijkl)
  {
    int i,ii,j,jj,k,l,m ;
    double s,t ;

    int ij=ijkl/30082;
    int kl=ijkl-30082*ij;

    i=((ij/177)%177)+2 ;
    j=(ij%177)+2 ;
    k=((kl/169)%178)+1 ;
    l=kl%169 ;
    for (ii=0;ii<97;ii++) {
      s=0.0 ;
      t=0.5 ;
      for (jj=0;jj<24;jj++) {
	m=(((i*j)%179)*k)%179 ;
	i=j;
	j=k;
	k=m;
	l=(53*l+1)%169;
	if (((l*m)%64)>=32) s+=t;
	t*=0.5;
      }
      u[ii]=s;
    }
    c=362436.0/16777216.0;
    cd=7654321.0/16777216.0;
    cm=16777213.0/16777216.0;
    i97=96;
    j97=32;
  }
	


  ////////////////////     For powerlaw_discrete   ////////////////////////////////////////////////////////
  // The look-up-table of P(x)=Pr(X>=x) = zeta(alpha,x)/zeta(alpha,xmin) 
  // Store those P(x) values in Pvector.
  bool flag_setPvector;
  vector<double> P;

  //Hurwitz zeta function (brute-forcely)
  // But this convergence is extremely slow, especially for small s.
  double Hzeta(double s, double q)
  {
    double sum = 0;
    double nterm = 0;
    int n;
    for(n=0; ;n++)
      {
	nterm = pow((n+q),-s);
	sum += nterm;

	if(fabs(nterm/sum) < TOLERANCE)
	  break;
      }
    cout << "nmax= " << n << endl; //test
    return sum;
  }

  void SetPvector(double alpha, int xmin)
  {
    // calculate Hzeta(alpha,xmin) first
    double a = gsl_sf_hzeta(alpha, xmin);//Hzeta(alpha, xmin);
    //cout << "a= " << a << endl;//test

    int N = 10000000+xmin; 

    P.resize(N,2);
    // Here we'd better set all P(x)=2 (which is unphysical) for x<xmin
    // becauce in performing binary search we would like P[] to be monotonic
    // if we set all P(x)=0 for x<xmin, since P(xmin)=1, then P(x)<1 for x>xmin
    // then P[] is not monotonic

    P[xmin] = 1;  
    for(int i=xmin;i<N-1;i++)
      {
	// P(x+1) = P(x) - 1/(x^alpha * zeta(alpha,xmin))
	P[i+1] = P[i] - 1/(pow((double)i,alpha)*a);  
      }
	    
    //for(int i=0;i<100;i++)
    //cout << P[i] << ',';


    flag_setPvector = true;
  }





/////////////////////////////////////////////////////////////////////////////////////////////////
/**
 * cpx_polylog_sum -- compute the polylogarithm by direct summation
 *
 * Li_s(z) = sum_{n=1}^infty z^n/ n^s
 * 
 * The magnitude of z must be less than one in order for the 
 * summation to be carried out.
 */
double polylog(double s, double z)
{
  double sum = 0;
  double zn = z;
  int nmax = 1;
  for(int n=1; ; n++) {
    double term = zn/pow((double)n,s);
    sum += term;
    zn *= z;
    if(term < TOLERANCE) {
      nmax = n;
      break;
    }
  }
  //cout << "n= " << nmax << endl; //test

  return sum;
}
/////////////////////////////////////////////////////////////////////////////////////////////////


  // for   int powerlaw_exponentialcutoff_discrete(double kappa)
  void SetPvector_SFEC(double gamma, double kappa)
  {
    double enkappa= exp(-1/kappa);
    double C = polylog(gamma, enkappa);

    int N = 10000000; 

    P.resize(N,2);
    // Here we'd better set all P(x)=2 (which is unphysical) for x<xmin
    // becauce in performing binary search we would like P[] to be monotonic
    // if we set all P(x)=0 for x<xmin, since P(xmin)=1, then P(x)<1 for x>xmin
    // then P[] is not monotonic

    P[1] = 1;  
    double a = enkappa;
    for(int i=1;i<=N;i++) {
      P[i+1] = P[i] - pow((double)i,-gamma) * a/ C ;
      a *= enkappa;
    }
	    
    flag_setPvector = true;
  }






   
};	

#endif
