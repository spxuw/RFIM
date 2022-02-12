#if !defined (HISTOGRAM_H)
#define HISTOGRAM_H

#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <vector>
#include <set>
#include <algorithm>

#include "Rand.h"
//#include "nr.h"

#include "Statistics.h"

using namespace std;


class ProbofW
{
 public:
    double w;       // w 
    double p;    // Prob(W=w)
};



//typedef vector<double> Vec_DP;

// P(k) = C(N-1,k) p^k (1-p)^(N-1-k)
double Binomial(int N, double p, int k);

class Bin
{
 public:
    double below;    // the below boundary of this bin
    double above;    // the above boundary
    double count;    // # of data in this bin

    double normalized_count;       // # of data in this bin/total # of data
    double sum;      // the sum of all points in the bin
    double ave;      // the average
    double var;      // the variance
    double err;      // the error of the mean: error bar !!!!!!!!!

    double fpv;      // fraction of positive values
};
// log bin for discrete distribution, for GKK model
void GetHistogram(const vector<int>& X, double ratio, char* fname);
void GetHistogram(const vector<int>& X, double ratio, char* fname, double& Kmean, double& K2mean);
void GetHistogram(const vector<int>& X, double ratio, char* fname, double& Kmean, double& K2mean, double& Kmax);

// linear bin for discrete distribution
void GetHistogram(const vector<int>& X);
void GetHistogram(const vector<int>& X, char* fname);
void GetHistogram(const vector<int>& X, char* fname, double& Kmean, double& K2mean, double& Kmax);

// linear bin for discrete distribution, for Erdos-Renyi model. p is the connection probability
void GetHistogram_ER(const vector<int>& X, double p, char* fname);
void GetHistogram_ER(const vector<int>& X, double p, char* fname, double& Kmean, double& K2mean, double& Kmax);


// linear bin for continuous distribution
void GetHistogram(const vector<int>& X, int Nbins);


// log bin for continous distribution, e.g. the distribution of betweeness
void GetHistogram(const vector<double>& X, double ratio, char* fname);
void GetHistogram_Accumulatively(const vector<int>& X, double ratio, char* label, char* fname);


// log bin for continous distribution, e.g. the distribution of betweeness
void GetHistogram(const vector<double>& X, double ratio, char* fname, double& Xmin);

// linear bin for continuous distribution, e.g. the eigenvalue distribution
void GetHistogram(const vector<double>& X, int Nbins, char* fname);
void GetPQw(const vector<double>& W, int Nbins, vector<ProbofW>& P, vector<ProbofW>& Q);
void GetPw(const vector<double>& W, int Nbins, vector<ProbofW>& P);


//
//template <typename T> 
//void GetDiscreteDistribution(vector<T>& X, char* fname);

void GetDiscreteDistribution(vector<double>& X, char* fname);
void GetDiscreteDistribution(vector<string>& X, char* fname);

void GetPz(vector<int>& K, vector<double>& P, double& z);
void GetPQz(vector<int>& K, vector<double>& P, vector<double>& Q, double& z);
void GetPQz(vector<int>& K, vector<double>& P, vector<double>& Q, double& z, double& v);
void GetPQz(vector<int>& K, vector<int>& Count, vector<double>& P, vector<double>& Q, double& z, double& v);

void GetPQz(vector<int>& K, vector<long double>& P, vector<long double>& Q, double& z, double& v);
void GetPQz(vector<int>& K, vector<int>& Count, vector<long double>& P, vector<long double>& Q, double& z, double& v);

//void GetPQz(vector<int>& Kin, vector<int>& Kout, vector<double>& P, vector<double>& Q, double& z);


void GetHistogram(const vector<int>& X, char* fname, const vector<double>& Y);
void GetHistogram(const vector<double>& X, int Nbins, char* fname, const vector<double>& Y);
void GetHistogram(const vector<int>& X, double ratio, char* fname, const vector<double>& Y);
void GetHistogram(const vector<double>& X, double ratio, char* fname, const vector<double>& Y);
void GetDiscreteDistribution(const vector<double>& X, char* fname, const vector<double>& Y);


class Bin2D
{
 public:
    double xbelow;    // the below boundary of this bin
    double xabove;    // the above boundary

    double ybelow;    // the below boundary of this bin
    double yabove;    // the above boundary

    double count;    // # of data in this bin

    double normalized_count;       // # of data in this bin/total # of data
    double sum;      // the sum of all points in the bin
    double ave;      // the average
    double var;      // the variance
    double err;      // the error of the mean: error bar !!!!!!!!!

    vector<double> vec;
};


void Get2DHistogram(const vector<double>& X, const vector<double>& Y, double ratiox, double ratioy, char* fname, const vector<double>& Z);
void Get2DHistogram_u(const vector<double>& X, const vector<double>& Y, double ratiox, double intervaly, char* fname, const vector<double>& Z);

void Get2DHistogram_u2(const vector<double>& X, const vector<double>& Y, double ratiox, int Nbins_y, char* fname, const vector<double>& Z);


#endif
