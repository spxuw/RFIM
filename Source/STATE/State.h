#ifndef _STATE_H_
#define _STATE_H_

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cmath>
#include <vector>
#include <list>
#include <queue> 
#include <deque>
#include <map>
#include <stack>
#include <set>
#include <limits.h>
#include <stack> 
#include <time.h>       /* time_t, time, ctime */



#include "Rand.h"
#include "types.h"
#include "main.h"
#include "find.h"
#include "color.h"
#include "Link.h"
#include "hypergeometric.h"
#include "Splitter.h"
#include "Histogram.h"
#include "Statistics.h"



using namespace std;

////////////////////////////// Note that we'd better use ``signed char" instead of "char" to guarantee it is signed. Note that ``char" is the smallest addressable unit of the machine that can contain basic character set. It is an integer type. Actual type can be either signed or unsigned depending on the implementation. YYL 01/14/2015
const char UP=1;
const char DOWN=-1;
const char ACTIVE = 1;
const char INACTIVE = -1;

/*
const char UP=1;
const char DOWN=0;
const char ACTIVE = 1;
const char INACTIVE = 0;
*/

const char NONE=2;
extern const hType   RESOLUTION;
#define J      RESOLUTION          // magnify the coupling constant 
#define CC     (4*J)               // capacity constant  

typedef vector<int> gs; // gene set



//////////////////////////////////////////////////////////////////////////////////////////////////////////
class module
{
 public:
  int     index;
  double  score0;       // ZA1 = sum/sqrt(k)
  double  score;        // ZA2 = sqrt(k) * (sum/(double)k - zave) 
  int     size;
  Nbl     nodelist; 
};
bool descending(double i, double j);
bool ModuleSortPredicate0(const module& M1, const module& M2);
bool ModuleSortPredicate(const module& M1, const module& M2);
//////////////////////////////////////////////////////////////////////////////////////////////////////////



//////////////////////////////////////////////////////////////////////////////////////////////////////////
class State 
{
 public:
    // constructor and destructor; see State_constructor.cpp 
    State();    // default constructor 
    State(int dimension, int length,  double disorder, int rngseed);  // construct a State without setting spin values
    State(int dimension, int length, double disorder, int rngseed, char dir); // State with all spins=DIR (UP 1 or DOWN 0)
    State(int dimension, int length, double disorder, int rngseed, double H); // State all spins=UP(DN) if H>0(<=0)
    State(const State &X); 
    ~State(); 


    // initialization, see State_initialize.cpp
    void Initialize_withdir(char dir);  // initialize a State spins with all spins=DIR (UP 1 or DOWN 0)
    void Initialize_fromfield(double H);  // initialize a State spins with all spins=UP(1) if H>0 or DOWN(0) if H<=0
    void Initialize_fromfile(int k);     // initialize State from the GS file: ./data/GS-jobID/k


    // set a default Random Field for all States, see State_setRandomField.cpp  
    static void SetRandomField(int dimension, int length, double disorder, int rngseed);
    //static void Set_RandomField_Cal_dEmin(int dimension, int length, double disorder, int rngseed);    

    void SetEffectiveField(double H);          // set the external field and the effective field 
    void SetEffectiveField_typeII(double lambda);


    // The following structures and functions are used to solve the Max-Flow problem; 
    list<myarc> ArcList;     // a list of arcs used in the network

    // simple mappling,      // see State_mapping.cpp
    void Mapping();                                    // mapping to a DIMACS-format file             
    void Mapping(long &n, long &m, long &s, long &t);  // mapping into ArcList

    // advanced mapping, using eariler solution (frozen spins);  see State_advanced_mapping.cpp
    map<int,int> NodeMap;  // the map between the DOWN spin's location (key) and the index of the nodes (value)
    vector<int> LocVec;    // the vector stores the map between the index of the nodes and the DOWN spin's location 
    bool mapping;          // whether the NodeMap has been built 
    void Node_Map_Loc(char FZdir); // set up the map between the free spin's location and the node index, 
    void Mapping_Frozenspins(char FZdir);                                     // mapping to a DIMACS-format file          
    void Mapping_Frozenspins(char FZdir,long &n, long &m, long &s, long &t);  // mapping into ArcList


    // calculate the internal energy of the state. see State_energy.cpp
    void   CalE0();   


    // inline small functions
    int    GetM()     const {return M;}
    double Getm()     const {return M/(double)N;}
    double GetE()     const {return E;}
    double GetE0()    const {return E0;}
    int    Get_Nfs()  const {return Nfs;}
    int    Get_Nup()  const {return (N+M)/2;}
    int    Get_Ndn()  const {return (N-M)/2;}
    double Get_nup()  const {return 0.5*(1+m);}
    double Get_ndn()  const {return 0.5*(1-m);}
    double Get_dEmin() const {return dEmin;}        
    double Get_abshmin() const {return abshmin;}        
    double Get_hmin() const {return hmin;}        
    double Get_hmax() const {return hmax;}        
    double Get_havg() const {return havg;}        
    double Get_hsum() const {return hsum;}        
    double Get_Hext() const {return Hext;}        
    hType  Get_heff(int i) const {return heff[i];}

    // "Set" implies that spin[i] could have no initial value
    void SetUP(int loc) {spin[loc] = UP;}        
    void SetDN(int loc) {spin[loc] = DOWN;}   

    // "Flip" implies that spin[i] has been previously assigned an initial value
    void FlipUP(int loc) {spin[loc] = UP;}        
    void FlipDN(int loc) {spin[loc] = DOWN;}   
     
    void CalMbyNup(int Nup)  {M= 2*Nup-N; m= M/(double)N;} 
    void CalMbyNdn(int Ndn)  {M= N-2*Ndn; m= M/(double)N;} 

    void SetE0(char signed dir) {
      if(!network)  E0 = (dir==UP) ? (-N*D-hsum) : (-N*D+hsum);       // lattice with PBC=1
      else   	    E0 = (dir==UP) ? (-Wsum-hsum) : (-Wsum+hsum);     // network 
    } // set E0 according to m=+1 or -1

    void CalE(double H) {E=E0 - H*M;}

    void SetR(double disorder) {R= disorder;}

    // calculate the difference of two states, see State_avalanche.cpp
    friend void Cal_Avalanche(State& A, State& B, list<int>& avalanche);  // Attention: the symbol & is necessary. 
    friend void Show_Avalanche(State& A, State& B);  // Attention: the symbol & is necessary. 
    void Show_State();
    void Save_State(int k);

    // assign an index to a single avalanche 
    void DrawOneAvalanche(int index, list<int> & avalanche);

    // draw all the avalanches with different colors (according to their indeies) during the whole M(H) process
    void getRGBcolor(int avalancheindex, double& r,double& g,double& b);
    void DrawAllAvalanches(char* fname);



    // Most of the following tool-funcitons are adopted from the Hysteresis Code; see State_tools.cpp 
    int  GetL() const {return L;}
    int  GetN() const {return N;}
    int  GetZ() const {return Z;}
    int  GetD() const {return D;}
    int  GetNbLocs(int j) const {return neighborLocs[j];}
    char GetSpin(int i) const {return spin[i];}
    int  GetLoc(const int* coords) const;
    void GetCoords(int loc, int* newCoords) const;
 
    void GetNeighbors(int currentLoc) const;

    void Reset(); // Reset for the calculation of ground state at different external field



    //////////////////////////////  For Network Ony  ///////////////////////////////////////////////////// 
    State(int num, double disorder, int rngseed, char dir);  
    State(int num, double disorder, int rngseed);  
    State(int num, char dir);  
    State(int num, double disorder, char dir);  
    State(int num, double disorder);
    State(int num);  

    static void Set_RandomPvalue_RandomBond(int num, double kmean, double gamma, int rngseed);
  
    static void Set_RandomField_RandomBond(int num, double kmean, double gamma, double disorder, int rngseed);
    static void Set_RandomField_RandomBond_Lattice(int dim, int a, double disorder, int seed, bool PBC);
    static void Read_Fields_Bonds(string file);
    static void Randomize_Fields(int seed);
    static void Randomize_Bonds(int seed);


    void Cal_BackgroundScore(int Q, int kmax, int seed, string file);
    void Read_BackgroundScore(string file);
    void Get_BackgroundScore(int Q, int kmax, int seed, string file);
    void Get_BackgroundScore(string file);

    void   DFS(int u, vector<bool>& visited, Nbl& Component, char STATUS);
    double CalModuleScore(Nbl& NL, module& Module); 
    double CalModuleScore(Nbl& NL); 
    void   CalModuleScore(Nbl& NL, double& SA1, double& SA2); 


    double GetScore();
    int   CalModule(module& M);
    void   SaveModules(char* file);

    int   CalLCC(module& M);


    void Save_State(string file, int k);
    void Initialize_fromfile(string file, int k);  
    
    
    bool CheckNeighbors(char dir, int u); 
    void WriteGML(char* fname);
    void WriteGML(char dir, char* fname);
    void WriteGML(char dir, int smin, char* fname);
    void WriteGML_HightlightBenchmarkGenes(char dir, int smin, char* fname);
    void WriteGML_HightlightBenchmarkGenes_NearestNeighbors(char dir, int smin, char* fname);

    int    edgeindex(int i, int j);
    double edgeweight(int i, int j);

    void DFS();
    void DFS(int u, vector<bool>& visited, Nbl&  Component);

    int ReadGeneSetsFile(char* genesetsfile);
    void FisherExactTest(vector<double>& P);
    double FisherExactTest(gs& X, gs& Y);
    void WriteGeneSetsFile(char* fname);
    //////////////////////////////////////////////////////////////////////////////////////////////////////



    /////////////////////////    For single-spin flip dynamics     ///////////////////////////////////////
    void SingleSpinFlip_Initialize(int dist, double R, int seed);
    void FlipNext(char direction);
    int  FindNext(char direction);
    void FlipSpin(int loc, char direction);
    void HysteresisLoop(char* fname);
    void CalMetastableState(double H0);
    double GetH() const {return H;}
    int    GetAvalancheCount() const {return AvalancheCount;}
    //////////////////////////////////////////////////////////////////////////////////////////////////////


    /////////////////////////    For jActive                       ///////////////////////////////////////
    void HMS_Ideker(int itype, double Ti, double Tf, int seed, char* file);  
    //////////////////////////////////////////////////////////////////////////////////////////////////////


    /////////////////////////    For dmGWAS                        ///////////////////////////////////////
    int   GetHighestscoringCommonNeighbors(int i, int k);
    void  GreedyGrowth(Nbl& NL, double r);
    void  HMS_Jia(double r, char* file);
    void  CalOverlapping();
    bool  CheckFNSConnectivity(vector<int>& NV);
    void  CheckFNSConnectivity(char* elistfile, char* fnsfile);
    //////////////////////////////////////////////////////////////////////////////////////////////////////




 private:
    char* spin;                  // the ground state (GS)

    int* avalancheindex;         // the index (occuring time) of the avalanche during the M(H) process

    int D;                       // dimensionality
    int L;                       // the linear size of the hypercubic lattice
    int Z;	                 // Z=2D, The co-ordination number of the lattice   
    int N;                       // total number of spins in the system (lattice), added by Yang Liu 07/21/05   
    double R;                    // disorder parameter(standard deviation of the Gaussian) 
    int seed;                    // RNG seed

    int* size;                   
    int* stride;
    mutable int* neighborLocs;
    mutable int** neighborCoords;

    /////////////////////////   Static Members   /////////////////////////////////////////////////////////////
    static Vec h;                // The random fields are stored in a default vector
    static double hmin;          // the minimum of the random fields
    static double hmax;          // the maximum of the random fields
    static double havg;          // the average of the random fields,in N \to \infty limit, havg \to 0.
    static double hsum;          // the sum of the random fields = havg * N
    static double abshmin;       // the minimum of |h[i]|
    static double abshmax;       // the maximum of |h[i]|

    static double dEmin;         // the smallest energy difference between two ground states
    //////////////////////////////////////////////////////////////////////////////////////////////////////////


    double  Hext;                // uniform external field
    hType*  heff;                // rounded effective local field (hi+Hext) (integer) 

    double E0;                   // internal energy E0= - sum_<i,j> Si*Sj - sum_i h[i]*Si 
    double E;                    // total energy    E = E0 - Hext * M
    int    M;                    // magnetization   M = sum_i Si = nUP-nDN = 2*nUP-N = N-2*nDN
    double m;                    // order parameter m = M/N;
    int    Nfs;                  // number of free sites (=N - number of frozen spins)


    //////////////////////////////  For Network Only  ///////////////////////////////////////////////////// 
    bool             network; // tell if the system is a random network or a hypercubic lattice 
    vector<module>   Modules; // equivalent to active connected component defined in Ideker's algorithm
    vector<int>      moduleindex; // moduleindex[i] labels the index of the module that node i belongs to
    int              LCCindex;// the module index of the LCC (largest connected component of the active genes)

 
    static vector<gs>       benchmarkgenesets;  // vector of benchmark gene sets 

    static vector<string>   benchmarkgeneset_name;  // gene set name
    static vector<double>   benchmarkgeneset_pmin;  // the minimum p-value this gene set will obtain with increasing H 
    static vector<double>   benchmarkgeneset_pmin_H;// the corresponding H value 
    static vector<double>   benchmarkgeneset_pmin_LCC;// the minimum p-value this gene set will obtain with increasing H (consider only the active LCC) 
    static vector<double>   benchmarkgeneset_pmin_H_LCC;// the corresponding H value 


    static double           c; // mean degree
    static vector<string>   NodeName;  
    static map<string,int>  NodeMAP;
    static vector<double>   PVALUE;
    static vector<Nbl>         A; 
    static vector<int>         K; 
    static map<string,int>     MAP; 
    static map<string,double>  Weight; // The random bonds are stored in a default map
    static double              Wsum;   // the sum of the random bonds
    static vector<Link>        LINK;  

    static vector<double> MU1;    
    static vector<double> SIGMA1; 
    static vector<double> MU2;    
    static vector<double> SIGMA2; 

    vector<int>   Color; // color of a node indicates whether it is a benchmark gene (1) or not (0)



    int Ncc; // # of connected components
    vector<Nbl>  AllComponents; // connected component
    vector<int>  ccindex; // which connected component each node belongs to 
    list<int> LClist; // the vertex list of the largest component

    vector<bool> LC; // whether the node belongs to the largest component
    int          lcc_index; // the index of the largest component
    int          Nlc; // number of nodes in the largest connected component
    double       Nac; // average size of connected components

    double       nlc; // :=Nlc/N
    int          Elc; // number of edges in the largest connected component
    double       mlc; // edge density for the largest connected component 

    int          Kmaxlc; // the highest degree in the largest connected component
    int          leader; // we choose the leader to be the node with highest degree (within the largest connected component, of course)
    //////////////////////////////////////////////////////////////////////////////////////////////////////



    /////////////////////////    For single-spin flip dynamics     ///////////////////////////////////////
    int* nUp;                   // The number of neighbors pointing up for each spin is stored here
    queue<int> spinFlipQueue;   // A queue to store the spins that will flip in the next shell
    int time;                   // clock used in the propagation of avalanche
    int AvalancheCount;         // count the total # of avalanches
    double H;                   // external field 
    //////////////////////////////////////////////////////////////////////////////////////////////////////



    /////////////////////////    For jActive and dmGWAS            ///////////////////////////////////////
    Rand rand;
    //////////////////////////////////////////////////////////////////////////////////////////////////////




};
//////////////////////////////////////////////////////////////////////////////////////////////////////////



//////////////////////////////////////////////////////////////////////////////////////////////////////////
inline int State::edgeindex(int i, int j)
{
  stringstream sst; 
  sst << i << ">" << j;  
  return MAP[sst.str()];
}

//There are many different types of edge weights

// (1)  the edge weight is read from the weighted edge list 
inline double State::edgeweight(int i, int j) { stringstream sst; sst << i << ">" << j; return Weight[sst.str()];}
//inline double State::edgeweight(int i, int j) { return 1.0;}
//inline double State::edgeweight(int i, int j) { return ((h[i])+(h[j]))*(2*abshmax); }

// (2) Here the edge weight is defined to be related to node weights.
// This type of edge weight is used in this paper: 
// An integrative network algorithm identifies age-associated differential methylation interactome hotspots targeting stem-cell differentiation pathways 
// James West1,2, Stephan Beck3, Xiangdong Wang4 & Andrew E. Teschendorff1,2,4
//inline double State::edgeweight(int i, int j) { return (fabs(h[i])+fabs(h[j]))/(2*abshmax); }


// (3) Shall we try Jij = fabs(h[i] + h[j]) / (2*abshmax)?
// Then if h[i]>>0 and h[j]<<0, then Jij will be small.
//inline double State::edgeweight(int i, int j) { return fabs(h[i] + h[j])/(2*abshmax); }
// This case, Why Jij should be very small? node i is very likely to be a disease gene, node j is very unlikely to be a disease gene,
// their interaction is zero??? Does this make sense? This implies that only nodes with same sign will have interactions. This results in strong localization!!! 


// (4) Shall we try Jij = wij * fabs(h[i] + h[j]) / (2*abshmax)?
// inline double State::edgeweight(int i, int j) { stringstream sst; sst << i << ">" << j; return Weight[sst.str()]*fabs(h[i] + h[j])/(2*abshmax); }


// (5) Shall we try Jij = wij * [fabs(h[i]) + fabs(h[j])] / (2*abshmax)?
// inline double State::edgeweight(int i, int j) { stringstream sst; sst << i << ">" << j; return Weight[sst.str()]*(fabs(h[i]) + fabs(h[j]))/(2*abshmax); }




inline void State::FlipSpin(int loc, char direction)
{
  spin[loc] = direction;
  M += 2*direction;
}



//////////////////////////////////////////////////////////////////////////////////////////////////////////


void Statistics(const vector<double>& data, double& ave, double& var);
void Get_N_E_NameList(string file, int& num, int& numedge, vector<string>& namelist);  

void Char2String(char& c, string& s);
string Int2String(int i);


#endif /* ! _STATE_H_ */
