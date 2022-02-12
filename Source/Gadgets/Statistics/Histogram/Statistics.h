#if !defined (STATISTICS_H)
#define STATISTICS_H

#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <vector>
#include <algorithm>
#include "Rand.h"
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_bessel.h>

using namespace std;

const int         XMAX = 600;
const long double EXPXMAX = exp((double)XMAX);
const long double EXP     = exp(1.0);

double Beta(double t, double x, double Itot); 
double Gaussian(double x, double mean, double var);
long double Poisson(int lambda, int k);

void Statistics(const vector<double>& data, double& ave, double& squareave, double& var);

void Statistics(const vector<double>& data, double& max, double& min, double& ave, double& var, double& error);
void Statistics(const vector<int>& data, double& max, double& min, double& ave, double& var, double& error);
void Statistics(const vector<int>& data, double& xmax, double& xmin, double& xave, double& xvar, double& xerror, double& x2ave);

//void Statistics(const vector<int>& data, double& xmax, double& xmin, double& xave, double& xvar, double& xerror, double& x2ave, double& x1inverseave);

void Statistics(const vector<int>& data, double& Xmin, double& Xmax, 
		double& X1inverseave, //  <1/(X+1)>
		double& logXave,      //  < logX >
		double& Xsrave,       //  < X^0.5 > 
		double& Xave,         //  < X > 
		double& X2ave );      //  < x^2 > 

void moment(const vector<int>& data, double &ave, double &adev, double &sdev, double &var, double &skew, double &curt);

void StatisticalDispersion(vector<int>& X, double& H7, double& H8, double& H9, double& H10);

double GetXsecondmin(const vector<double>& data, double Xmin);

void Statistics(const vector<double>& data, 
		double& min, double& Q1, double& Q2, double& Q3, double& max, 
		double& ave, double& var, double& error,
		ofstream& fout1, ofstream& fout2);

void Binning(const vector<double>& data, int B, vector<double>& Binsize, vector<double>& Q);
int GetBinIndex(const double& x, const vector<double>& Q);

#endif


