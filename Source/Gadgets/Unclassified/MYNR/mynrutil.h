#ifndef _MYNR_H_
#define _MYNR_H_



#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <string>
#include <list>
#include <set>
#include <stack>
#include <deque>
#include <limits>








#include "Statistics.h"
#include "ULink.h"
#include "Link.h"

#include "nr.h"

#include "Rand.h"

using namespace std;

const double EPS = numeric_limits<double>::epsilon(); 
const double SMALL = numeric_limits<double>::min()/numeric_limits<double>::epsilon();


typedef list<int> Nbl; 

void mynrerror(const string error_text);



int* ivector(int nl, int nh);
double* dvector(int nl, int nh);
void free_ivector(int* v, int nl, int nh);
void free_dvector(double* v, int nl, int nh);



int** imatrix(int nrl, int nrh, int ncl, int nch);
double** dmatrix(int nrl, int nrh, int ncl, int nch);
void free_imatrix(int** m, int nrl,int nrh, int ncl, int nch);
void free_dmatrix(double** m, int nrl,int nrh, int ncl, int nch);



int*** imatrix(int xi, int xf, int yi, int yf, int zi, int zf);
double*** dmatrix(int xi, int xf, int yi, int yf, int zi, int zf);
void free_imatrix(int*** m, int xi, int xf, int yi, int yf, int zi, int zf);
void free_dmatrix(double*** m, int xi, int xf, int yi, int yf, int zi, int zf);



int** imatrix(int M, int* T);
double** dmatrix(int M, int* T);
void free_imatrix(int** X, int M);
void free_dmatrix(double** X, int M);



int*** imatrix(int M, int* T, int d);
double*** dmatrix(int M, int* T, int N);
void free_imatrix(int*** X, int M, int* T);
void free_dmatrix(double*** X, int M, int* T);




void GenerateCovarianceMatrix(double**& C, int d, int seed);
void CholeskyDecomposition(double**& C, double**& L, int d);

void PrintMatrix(double**& A, int d);
void MatrixCopy(double**& A, double**& B, int d);
void MatrixCopy(double**& A, double**& B, int n, int m);
void MatrixTranspose(double**& A, double**& B, int d);
void MxM(double**& A, double**& B, double**& C, int d);
void MxV(double**& A, double*& B, double*& C, int d);
void MxM(double**& A, double**& B, double**& C, int n, int m, int q);
void VaV(double*& A, double*& B, double*& C, int d);
void VsV(double*& A, double*& B, double*& C, int d);
double VdotV(double*& A, double*& B, int d);

void LUD(double**& a, int*& indx, double& s, int d);
void LUbksb(double**& a, int*& indx, double*& b, int d);
void MatrixInverse(double**& a, double**&y, double& det, int d);
void Separate(double**& a, double**& xl, double**& xu, int d);


int BisectionSearch(vector<double>& xx, const double x);
int BisectionSearch(double* xx, int n, const double x);

bool find(vector<set<int> >& Vset, set<int>& s);
bool find(vector< vector<char>  >& Vstate, vector<char>& state);
bool find(vector< vector<int> >& Vv, vector<int>& v);
bool find(vector< vector<string> >& Vs, vector<string>& v);
bool find(list<int>& X, int j);
bool find(vector<string>& X, string j);
bool find(vector<set<int> >& X, set<int>& j);
bool find(vector<list<int> >& X, list<int>& j);


bool find(set<int>& nbl, int j);

bool find(Nbl& nbl, int j);
int findpos(Nbl& nbl, int j);

bool find(vector<int>& nbv, int j);
int findpos(vector<int>& nbv, int j);

bool find(vector<string>& nbs, string s);
int findpos(vector<string>& nbs, string s);

bool find(deque<int>& S, int x);
void print(deque<int>& Q);

bool find(list<ULink>& L, ULink e);
bool find(list<Link>& L, Link e);


int getnodeinNbl(Nbl& nbl, int j);
double getnodeinlist(list<double>& X, int index);
int getnodeinlist(list<int>& X, int index);






template <typename Iter> 
void getNbl(Iter begin, Iter end); 

void GetNumberFromaString(string str, int& x); 

double polylog(double s, double z);
double polylog(double s, double z, int prec);

int FindRoots(double fx(const double), double X1, double X2, double& rootmin, double& rootmax, vector<double>& roots);
int FindRoots_long(long double fx(const long double), long double X1, long double X2, long double& rootmin, long double& rootmax, vector<long double>& roots);

int FindRoots_u(double fx(const double), double X1, double X2, vector<double>& roots);


void zbrak_rrg(double fx(const int r, const double x), const int r, const double x1, const double x2, const int n,
	       Vec_O_DP &xb1, Vec_O_DP &xb2, int &nroot);
double rtbis_rrg(DP func(const int, const DP), const int r, const DP x1, const DP x2, const DP xacc);
int FindRoots_RRG(double frx(const int, const double), const int r, double X1, double X2, vector<double>& roots);


int FindRoots(double fx(const double, const double), const double a, double X1, double X2, double& rootmin, double& rootmax, vector<double>& roots);


int FindRoots(double fx(const double, const double, const double), const double a, const double b, double X1, double X2, double& rootmin, double& rootmax, vector<double>& roots);




double binom(int m, int k);

int print(Nbl& nbl);





double NR_gammq(const double a, const double x);





double GetPoisson(double lambda, int k);
void   Get_ER_concavemean(double lambda, double& logk1, double& ksr);
double Get_ER_AbsoluteDeviation(double z);
double Get_ER_ShannonEntropy(double z);
double Get_ER_RMD(double z);
double Get_ER_QCD(double z);


double Get_Gammas_StaticModel(double z, double a, int k);
double Get_Gammas_StaticModel_Adaptive(double z, double a, int k);

void   Get_SF_SM_concavemean(double z, double a, double& logk1, double& ksr);
double Get_SF_SM_AbsoluteDeviation(double z, double a);
double Get_SF_SM_ShannonEntropy(double z, double a);
double Get_SF_SM_RMD(double z, double a);
double Get_SF_SM_QCD(double z, double a);

double Get_Pk_ChungLu(vector<double>& w, int k);
double Get_CL_RMD(vector<double>& P, double z);


double Get_SF_EC_RMD(double gamma, double kappa, double z);
double Get_SF_EC_D_RMD(double gamma, double kappa, double z);

double Get_SF_RMD(double gamma, double z);
double Get_SF_D_RMD(double gamma, double z);

double Get_EXP_D_RMD(double enkappa, double z);

void STP(double* X, int a, double* Y, int b, double* C);
void STP(double** A, int m, int n,  double** B, int p, int q, double** C);
void print(double** C, int a, int b); 



void MatrixFormIII_check_test(int n, int m);
void MatrixPrint(double** X, int a, int b);
void MatrixPrint(int** X, int a, int b);
void MxM(int**& A, int**& B, int**& C, int d);

void MatrixCopy(double**& A, double**& B, int d);
void MatrixCopy(double**& A, double**& B, int n, int m);
void MatrixCopy(int**& A, int**& B, int n, int m);

bool Transform_to_MatrixFormIII(int**& Q, int n, int m, int**& rowindex, int**& colindex);








double Norm(vector<double>& V);
double Sum(vector<double>& V);
double Error(vector<double>& V1, vector<double>& V2);



int GetOverlap(set<int>& S1, set<int>& S2);

void Char2String(char& c, string& s);
string Int2String(int i);


#endif
