#if !defined (RAND_H)
#define RAND_H

#include <iostream>
#include <cmath>
using namespace std;



const int BIGPRIME = 899999963;
const double PI             = 3.1415926535897932384626433832795028841971693993751;

class Rand {
	double c,cd,cm,u[97];
	int i97,j97 ;
	bool outputReady;
	double output;

public:
	Rand(int seed) {
		ranmarin(abs(seed%BIGPRIME));
		outputReady=false;
	}

	Rand() {
		ranmarin(1);
		outputReady=false;
	}

	void seed(int seed) {
		ranmarin(abs(seed%BIGPRIME));
		outputReady=false;
	}

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

	int discrete(int min, int max)
	{
		int ret= (int)(floor(uniform()*(max+1-min)+min));
		return ret;
	}


	double gaussian(double sd)
	{
	  double x,y,r,z;
	  if(outputReady)
	    {
	      outputReady = false;
	      return output*sd;
	    }

	  do
	    {
	      x=(uniform()*2.0)-1.0;
	      y=(uniform()*2.0)-1.0;
	      r=x*x+y*y;
	    } while (r>=1.0);

	  z=sqrt(-2.0*log(r)/r);
	  output = x*z;
	  outputReady=true;
	  return y*z*sd;
	}

	double lorentzian(double sd) {
	  double x,y,z;
	  do {
	    x = uniform();
	    y = 2.0 * uniform() - 1.0;
	  } while ((x*x+y*y)>1.0);
	  z = y/x;
	  return sd*z/2.0;
	}


    double bimodal(double R)
	{
	    return (uniform()>0.5) ? (R) : (-R);
	}


    double rectangular(double R)
	{
	    return ((uniform()*2-1)*R); 
	}



    double parabolic(double R)
      {
	  double x = uniform();
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

	  double y;
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



private:
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
	    
};	

#endif
