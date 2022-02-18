#if !defined (RAND_H)
#define RAND_H

#include <cstdlib>
#include <iostream>
#include <cmath>
#include <vector>
#include <list>
#include <set>
#include <limits>
#include <algorithm>

#include <gsl/gsl_sf_zeta.h>

using namespace std;

const double TOLERANCE = numeric_limits<double>::epsilon();
const int    INTMAX    = numeric_limits<int>::max();

const int BIGPRIME = 899999963;
const double PI    = 3.1415926535897932384626433832795028841971693993751;

typedef double (*FUNC)(const double, const double, const double);
inline double fx(const double htilde, const double y, const double u) 
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

  return -1;
}

class Rand {
  double c,cd,cm,u[97];
  int i97,j97 ;
  bool outputReady;
  double output;

  int idum; 

 public:
  Rand(int seed) {
    ranmarin(fabs(seed%BIGPRIME));
    outputReady=false;
    idum=-1*seed; 

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
    idum=-1*seed; 

    flag_setFvector = false;
    flag_setPvector = false;
  }


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



  int discrete(int min, int max)
  {
    int ret= (int)(floor(ran1()*(max+1-min)+min));
    return ret;
  }


  void discrete(int min, int max, int& a, int& b)
  {
    a = (int)(floor(ran1()*(max+1-min)+min));
    
    do {
      b = (int)(floor(ran1()*(max+1-min)+min));
    } while (a==b);
    
  }


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


  void discrete(int min, int max, int n, vector<int>& Vec)
  {
    set<int> Set;
 
    Vec.clear();
    for(int i=0; i<n; ) {
      int a = (int)(floor(ran1()*(max+1-min)+min));
      if(Set.find(a)==Set.end()) {
	Set.insert(a);
	Vec.push_back(a);
	i++;
      }
    }
  }



  int discrete_logbin(int xmin, int xmax)
  {
    double ratio = 1.15;
    if (xmin==0)
      xmin=1;

    int Nbins = (int)(floor(log(xmax/xmin)/log(ratio))) + 1; 
    int i = discrete(0, Nbins-1);
    int a = xmin * pow(ratio, (double)i);
    int b = a*ratio;
    if(b>xmax)
      b = xmax;
    
    int x = discrete(a,b);
    return x;
  }



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

  void randomsubset_knuth(int n, int m, vector<int>& S)
  {    
    for (int i = 0; i < n; i++)
      if ((discrete(1,n) % (n-i)) < m) {
	S.push_back(i);
	m--;
      }
  }


  void randomsubset_shuf(int n, int m, vector<int>& S)
  {   
    int i,j,t;
    vector<int> x(n);
    for (i = 0; i < n; i++)
      x[i] = i;
    for (i = 0; i < m; i++) {
      j = discrete(i, n-1);
      t = x[i]; 
      x[i] = x[j]; 
      x[j] = t;
    }
    sort(x.begin(), x.begin()+m);
    for (i = 0; i< m; i++)
      S.push_back(x[i]);
  }

  void randomsubset_shuffle(int n, int m, vector<int>& S)
  {   
    vector<int> x(n);
    for(int i=0; i<n; i++)
      x[i] = i;

    random_shuffle(x.begin(), x.end());

    for (int i=0; i<m; i++) 
      S.push_back(x[i]);
  }
  

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


	v1=(uniform()*2.0)-1.0;
	v2=(uniform()*2.0)-1.0;

	Rsq=v1*v1+v2*v2;
      } while (Rsq>=1.0);

    z=sqrt(-2.0*log(Rsq)/Rsq);
    output = v1*z*sd;
    outputReady=true;
    return v2*z*sd;
  }

  double lorentzian(double sd) 
  {
    double v1,v2,z;
    do {
      v2 = uniform();
      v1 = 2.0 * uniform() - 1.0;
    } while ((v1*v1+v2*v2)>1.0);
    z = v2/v1;
    return sd*z/2.0;
  }

  double powerlaw_continuous(double alpha, double xmin)
  {
    return ( xmin*pow(1-ran1(), 1./(1-alpha) ) );
  }

  double powerlaw_continuous_interval(double a, double beta)
  {
    return pow( (1- (1-ran1())*(1-pow(a, 1-beta))), 1/(1-beta) );
  }


  int powerlaw_discrete_approximate(double alpha, double xmin)
  {
    double r = ran1();
    int k = (int)( (xmin-0.5)*pow(1-r, 1./(1-alpha) ) + 0.5 );
    if(k<0) {
      cout << 1-r << ',' << k << endl;
      k =  INTMAX;
    }
    return k;

  }


  int powerlaw_discrete(double alpha, int xmin)
  {
    if(!flag_setPvector)
      SetPvector(alpha,xmin);

    double r = ran1();
    int j = Findj(P, 1-r);
    return j;

  }
  





  int powerlaw_exponentialcutoff_discrete(double gamma, double kappa)
  {
    if(!flag_setPvector)
      SetPvector_SFEC(gamma, kappa);

    double r = ran1();

    int j = Findj(P, 1-r);
    return j;

  }
  
  int arbitrary_discrete(vector<double>& CDF)
  {
    double r = ran1();
    int j = Findj(CDF, r);
    return j+1;
  }

  int arbitrary_discrete(vector<long double>& CDF)
  {
    double r = ran1();
    int j = Findj(CDF, r);
    return j+1;
  }



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
	v1=(ran1()*2.0)-1.0; 
	v2=(ran1()*2.0)-1.0; 
	Rsq=v1*v1+v2*v2;
      } while (Rsq>=1.0);

    z=sqrt(-2.0*log(Rsq)/Rsq);
    output = exp(v1*z*sd);
    outputReady=true;
    return exp(v2*z*sd);
  }


 
  double bimodal(double R)
  {
    
    return (ran1()>0.5) ? (R) : (-R);
  }

 
  double rectangular(double R)
  {
    
    return ((ran1()*2-1)*R); 
  }


 
  
  
  
  
  

  
  
  
  
  
  
  double parabolic(double R)
  {
    
    double x = ran1();

    int sign =  (x>0.5) ? (-1) : (+1);

    double R2 = R*R;
    double R3 = R2*R;

    double q     = 4*R3*(x-0.5);
    double p     = -3*R2;
    double Disc  = q*q/4. + p*p*p/27.;
    double r     = R3;
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
	
    

    return y;
  }



  
  
  
  
  
  
  

  double antiparabolic(double R)
  {
    
    double x = ran1();

    int sign =  (x>0.5) ? (+1) : (-1);

    double h = sign * pow(sign*(2*x-1),1./3.) * R;
    return h;
  }




  
  
  
  
  
  
  

  
  
  
  
  
  
  double DMI_lookuptable(double y, double R)
  
  {
    if(!flag_setFvector)
      SetFvector(y,R);

    double u = ran1();
    int j = Findj(F, u);
    double h = -R + j*dh;

    return h;
  }


  
  
  double DMI(double y, double R)
  {
    double u = ran1();
    int sign =  (u>0.5) ? (+1) : (-1);

    
    
    
    

    const double hacc = 1e-12;
    double htilde = rtbis(&fx, y, u, 0, 1, hacc);
    double h = sign*htilde*R;
    return h;
  }
    


  
  
  
  
  
  
  int exponential(double kappa)
  {
    double u = ran1();
    
    int k = -kappa * log(1-u);  
    return k;
  }
  
 



 private:

  
  


  
  
  
  
  

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
	    
    

    F[0]     = 0;
    F[2*N-2] = 1 - F[0];
    for(int i=1;i<N;i++)
      {
	double h1 = -R + i*dh;
	double h2 = h1 + dh;
		
	
	F[i]       = F[i-1] + 0.5* (P_DMI(h1,y,R)+P_DMI(h2,y,R))*dh;
	F[2*N-2-i] = 1- F[i];
      }
	    
    flag_setFvector = true;
  }

  
  
  
  
  
  int Findj(vector<double>& xx, const double x)
  {
    int j,jl,jm,ju;
    int n=xx.size();

    jl=-1; ju=n;  

    bool ascnd = (xx[n-1] >= xx[0]);   
    while (ju-jl > 1) 
      {
	jm=(ju+jl) >> 1;  
	if ((x >= xx[jm]) == ascnd) jl=jm; 
	else    	          ju=jm; 
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

    jl=-1; ju=n;  

    bool ascnd = (xx[n-1] >= xx[0]);   
    while (ju-jl > 1) 
      {
	jm=(ju+jl) >> 1;  
	if ((x >= xx[jm]) == ascnd) jl=jm; 
	else    	          ju=jm; 
      }
	    
    if (x == xx[0]) j=0;
    else if (x == xx[n-1]) j=n-2;
    else j=jl;

    return j;
  }

  
  


  
  
  


  
  
  double rtbis(FUNC func, const double y, const double u, const double x1, const double x2, const double xacc)
  {
    const int JMAX=40; 
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
	


  
  
  
  bool flag_setPvector;
  vector<double> P;

  
  
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
    cout << "nmax= " << n << endl; 
    return sum;
  }

  void SetPvector(double alpha, int xmin)
  {
    
    double a = gsl_sf_hzeta(alpha, xmin);
    

    int N = 10000000+xmin; 

    P.resize(N,2);
    
    
    
    

    P[xmin] = 1;  
    for(int i=xmin;i<N-1;i++)
      {
	
	P[i+1] = P[i] - 1/(pow((double)i,alpha)*a);  
      }
	    
    
    


    flag_setPvector = true;
  }


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
  

  return sum;
}



  
  void SetPvector_SFEC(double gamma, double kappa)
  {
    double enkappa= exp(-1/kappa);
    double C = polylog(gamma, enkappa);

    int N = 10000000; 

    P.resize(N,2);
    
    
    
    

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
