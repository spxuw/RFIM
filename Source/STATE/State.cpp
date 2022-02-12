#include "State.h"

using namespace std;
 

//////////////////////////////////////////////////////////////////////////////////////////////////////
// Note that those static members must be defined in the .cpp file.
// 
// Without setting the size of the default vector to be large enough
// there will be Segmentation fault. 
// I have no idea what's going on here and what is the lower limit of the 
// size you have to set. So I just set it to be 256^3.
Vec State::h; //(16777216);  
double State::hmin=0;
double State::hmax=0;
double State::havg=0;
double State::hsum=0;
double State::abshmin=0;
double State::abshmax=0;
double State::dEmin=0;

double             State::c=0; // mean degree
 
vector<Nbl>        State::A;//(16777216);
vector<int>        State::K;//(16777216);
vector<Link>       State::LINK;//(16777216);
vector<string>     State::NodeName;//(16777216);
map<string,int>    State::NodeMAP; 
vector<double>     State::PVALUE;

vector<gs>       State::benchmarkgenesets;  
vector<string>   State::benchmarkgeneset_name;  // gene set name
vector<double>   State::benchmarkgeneset_pmin;  // the minimum p-value this gene set will obtain with increasing H 
vector<double>   State::benchmarkgeneset_pmin_H;// the corresponding H value 
vector<double>   State::benchmarkgeneset_pmin_LCC;// the minimum p-value this gene set will obtain with increasing H (consider only the active LCC) 
vector<double>   State::benchmarkgeneset_pmin_H_LCC;// the corresponding H value 


vector<double>     State::MU1;//(16777216);    
vector<double>     State::SIGMA1;//(16777216); 
vector<double>     State::MU2;//(16777216);    
vector<double>     State::SIGMA2;//(16777216); 
map<string,int>    State::MAP;
map<string,double> State::Weight;
double             State::Wsum=0;
//////////////////////////////////////////////////////////////////////////////////////////////////////

bool descending(double i, double j) {return (i>j);}
bool ModuleSortPredicate0(const module& M1, const module& M2) {return M1.score0 > M2.score0;}
bool ModuleSortPredicate(const module& M1, const module& M2) {return M1.score > M2.score;}
