
#include "nrutil.h"

float **y,*xx;

void rkdumb(vstart,nvar,x1,x2,nstep,derivs)
float vstart[],x1,x2;
int nstep,nvar;
void (*derivs)();
{
	void rk4();
	int i,k;
	float x,h;
	float *v,*vout,*dv;

	v=vector(1,nvar);
	vout=vector(1,nvar);
	dv=vector(1,nvar);
	for (i=1;i<=nvar;i++) {
		v[i]=vstart[i];
		y[i][1]=v[i];
	}
	xx[1]=x1;
	x=x1;
	h=(x2-x1)/nstep;
	for (k=1;k<=nstep;k++) {
		(*derivs)(x,v,dv);
		rk4(v,dv,nvar,x,h,vout,derivs);
		if ((float)(x+h) == x) nrerror("Step size too small in routine rkdumb");
		x += h;
		xx[k+1]=x;
		for (i=1;i<=nvar;i++) {
			v[i]=vout[i];
			y[i][k+1]=v[i];
		}
	}
	free_vector(dv,1,nvar);
	free_vector(vout,1,nvar);
	free_vector(v,1,nvar);
}
