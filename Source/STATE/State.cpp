#include "State.h"

using namespace std;
 








Vec State::h; 
double State::hmin=0;
double State::hmax=0;
double State::havg=0;
double State::hsum=0;
double State::abshmin=0;
double State::abshmax=0;
double State::dEmin=0;

double             State::c=0; 
 
vector<Nbl>        State::A;
vector<int>        State::K;
vector<Link>       State::LINK;
vector<string>     State::NodeName;
map<string,int>    State::NodeMAP; 
vector<double>     State::PVALUE;

vector<gs>       State::benchmarkgenesets;  
vector<string>   State::benchmarkgeneset_name;  
vector<double>   State::benchmarkgeneset_pmin;  
vector<double>   State::benchmarkgeneset_pmin_H;
vector<double>   State::benchmarkgeneset_pmin_LCC;
vector<double>   State::benchmarkgeneset_pmin_H_LCC;


vector<double>     State::MU1;
vector<double>     State::SIGMA1;
vector<double>     State::MU2;
vector<double>     State::SIGMA2;
map<string,int>    State::MAP;
map<string,double> State::Weight;
double             State::Wsum=0;


bool descending(double i, double j) {return (i>j);}
bool ModuleSortPredicate0(const module& M1, const module& M2) {return M1.score0 > M2.score0;}
bool ModuleSortPredicate(const module& M1, const module& M2) {return M1.score > M2.score;}
