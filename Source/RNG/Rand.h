/*
 * Copyright (C) 1998 Matthew C. Kuntz, James P. Sethna, Karin A. Dahmen
 * and John Carpenter
 *
 * This file is part of the Hysteresis program.
 *
 * The Hysteresis program is free software; you can redistribute it
 * and/or modify it under the terms of the GNU General Public License
 * version 2 as published by the Free Software Foundation.
 *
 * See the file COPYING for details.
 */


#if !defined (RAND_H)
#define RAND_H

#include <iostream>
#include <cmath>
using namespace std;

//#include <cstdlib>


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

	// Gives all numbers weighted equally, min and max inclusive
	int discrete(int min, int max)
	{
		int ret= (int)(floor(uniform()*(max+1-min)+min));
//		assert(ret<=max && ret>=min);
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

	/* John Carpenter added */
	double lorentzian(double sd) {
	  double x,y,z;
	  do {
	    x = uniform();
	    y = 2.0 * uniform() - 1.0;
	  } while ((x*x+y*y)>1.0);
	  z = y/x;
	  return sd*z/2.0;
	}


    /* Yang Liu added */
    double bimodal(double R)
	{
	    return (uniform()>0.5) ? (R) : (-R);
	}

    /* Yang Liu added */
    double rectangular(double R)
	{
	    return ((uniform()*2-1)*R); 
	}

    /* Yang Liu added */
    // This parabolic distribution is obtained by using the transformation method. See NR:Ch7.2, Pg.291.
    // normalized parabolic distribution:  P(h) = (R^2-h^2)/(4R^3/3) with h \in [-R,R]
    // F(h) = \int_0^h p(h) dh = 1/2 - h^3/(4R^3) + 3h/(4R)
    // given x \in [0,1), we want to find the root(s) of F(h) = x, i.e. 
    //              h^3 - 3R^2 h + 4 R^3(x-1/2) == 0

    // Basically we are sovling a cubic equation.
    // Note that here the Discriminant (Disc) for this simple cubic equation is negative (crucially due to x\in[0,1))
    // which means there are three real roots.
    // Now it is tricky to apply the analytic formula, e.g see http://en.wikipedia.org/wiki/Cubic_equation
    // Because Discriminant is negatice, we have to explicitly write the three real roots rather than 
    // directly applying those formulas. Note that the C++ cannot handle sqrt(negative number). It will give you "nan"

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
	
	  //cout << x << ' ' << h1 << ' ' << h2 << ' ' << h3 << ' ' << y << endl;

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
