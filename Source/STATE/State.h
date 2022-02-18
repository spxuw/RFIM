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
#include <time.h>      



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

const char UP=1;
const char DOWN=-1;
const char ACTIVE = 1;
const char INACTIVE = -1;


const char NONE=2;
extern const hType   RESOLUTION;
#define J      RESOLUTION          
#define CC     (4*J)               

typedef vector<int> gs; 




class module
{
 public:
  int     index;
  double  score0;       
  double  score;        
  int     size;
  Nbl     nodelist; 
};
bool descending(double i, double j);
bool ModuleSortPredicate0(const module& M1, const module& M2);
bool ModuleSortPredicate(const module& M1, const module& M2);





class State 
{
 public:
    State();    
    State(int dimension, int length,  double disorder, int rngseed);  
    State(int dimension, int length, double disorder, int rngseed, char dir); 
    State(int dimension, int length, double disorder, int rngseed, double H); 
    State(const State &X); 
    ~State(); 


    void Initialize_withdir(char dir);  
    void Initialize_fromfield(double H);  
    void Initialize_fromfile(int k);     


    static void SetRandomField(int dimension, int length, double disorder, int rngseed);

    void SetEffectiveField(double H);          
    void SetEffectiveField_typeII(double lambda);


    list<myarc> ArcList;     

    void Mapping();                                    
    void Mapping(long &n, long &m, long &s, long &t);  

    map<int,int> NodeMap;  
    vector<int> LocVec;    
    bool mapping;          
    void Node_Map_Loc(char FZdir); 
    void Mapping_Frozenspins(char FZdir);                                     
    void Mapping_Frozenspins(char FZdir,long &n, long &m, long &s, long &t);  


    void   CalE0();   


    
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

    void SetUP(int loc) {spin[loc] = UP;}        
    void SetDN(int loc) {spin[loc] = DOWN;}   

    void FlipUP(int loc) {spin[loc] = UP;}        
    void FlipDN(int loc) {spin[loc] = DOWN;}   
     
    void CalMbyNup(int Nup)  {M= 2*Nup-N; m= M/(double)N;} 
    void CalMbyNdn(int Ndn)  {M= N-2*Ndn; m= M/(double)N;} 

    void SetE0(char signed dir) {
      if(!network)  E0 = (dir==UP) ? (-N*D-hsum) : (-N*D+hsum);       
      else   	    E0 = (dir==UP) ? (-Wsum-hsum) : (-Wsum+hsum);     
    } 

    void CalE(double H) {E=E0 - H*M;}

    void SetR(double disorder) {R= disorder;}

    friend void Cal_Avalanche(State& A, State& B, list<int>& avalanche);  
    friend void Show_Avalanche(State& A, State& B);  
    void Show_State();
    void Save_State(int k);

    void DrawOneAvalanche(int index, list<int> & avalanche);

    void getRGBcolor(int avalancheindex, double& r,double& g,double& b);
    void DrawAllAvalanches(char* fname);



    int  GetL() const {return L;}
    int  GetN() const {return N;}
    int  GetZ() const {return Z;}
    int  GetD() const {return D;}
    int  GetNbLocs(int j) const {return neighborLocs[j];}
    char GetSpin(int i) const {return spin[i];}
    int  GetLoc(const int* coords) const;
    void GetCoords(int loc, int* newCoords) const;
 
    void GetNeighbors(int currentLoc) const;

    void Reset(); 



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
    



    void SingleSpinFlip_Initialize(int dist, double R, int seed);
    void FlipNext(char direction);
    int  FindNext(char direction);
    void FlipSpin(int loc, char direction);
    void HysteresisLoop(char* fname);
    void CalMetastableState(double H0);
    double GetH() const {return H;}
    int    GetAvalancheCount() const {return AvalancheCount;}


    void HMS_Ideker(int itype, double Ti, double Tf, int seed, char* file);  
    


    int   GetHighestscoringCommonNeighbors(int i, int k);
    void  GreedyGrowth(Nbl& NL, double r);
    void  HMS_Jia(double r, char* file);
    void  CalOverlapping();
    bool  CheckFNSConnectivity(vector<int>& NV);
    void  CheckFNSConnectivity(char* elistfile, char* fnsfile);




 private:
    char* spin;                  

    int* avalancheindex;         

    int D;                       
    int L;                       
    int Z;	                 
    int N;                       
    double R;                    
    int seed;                    

    int* size;                   
    int* stride;
    mutable int* neighborLocs;
    mutable int** neighborCoords;

    
    static Vec h;                
    static double hmin;          
    static double hmax;          
    static double havg;          
    static double hsum;          
    static double abshmin;       
    static double abshmax;       

    static double dEmin;         


    double  Hext;                
    hType*  heff;                

    double E0;                   
    double E;                    
    int    M;                    
    double m;                    
    int    Nfs;                  


    bool             network; 
    vector<module>   Modules; 
    vector<int>      moduleindex; 
    int              LCCindex;

 
    static vector<gs>       benchmarkgenesets;  

    static vector<string>   benchmarkgeneset_name;  
    static vector<double>   benchmarkgeneset_pmin;  
    static vector<double>   benchmarkgeneset_pmin_H;
    static vector<double>   benchmarkgeneset_pmin_LCC;
    static vector<double>   benchmarkgeneset_pmin_H_LCC;


    static double           c; 
    static vector<string>   NodeName;  
    static map<string,int>  NodeMAP;
    static vector<double>   PVALUE;
    static vector<Nbl>         A; 
    static vector<int>         K; 
    static map<string,int>     MAP; 
    static map<string,double>  Weight; 
    static double              Wsum;   
    static vector<Link>        LINK;  

    static vector<double> MU1;    
    static vector<double> SIGMA1; 
    static vector<double> MU2;    
    static vector<double> SIGMA2; 

    vector<int>   Color; 



    int Ncc; 
    vector<Nbl>  AllComponents; 
    vector<int>  ccindex; 
    list<int> LClist; 

    vector<bool> LC; 
    int          lcc_index; 
    int          Nlc; 
    double       Nac; 

    double       nlc; 
    int          Elc; 
    double       mlc; 

    int          Kmaxlc; 
    int          leader; 
    



    
    int* nUp;                   
    queue<int> spinFlipQueue;   
    int time;                   
    int AvalancheCount;         
    double H;                   
    



    
    Rand rand;
    




};





inline int State::edgeindex(int i, int j)
{
  stringstream sst; 
  sst << i << ">" << j;  
  return MAP[sst.str()];
}




inline double State::edgeweight(int i, int j) { stringstream sst; sst << i << ">" << j; return Weight[sst.str()];}



























inline void State::FlipSpin(int loc, char direction)
{
  spin[loc] = direction;
  M += 2*direction;
}






void Statistics(const vector<double>& data, double& ave, double& var);
void Get_N_E_NameList(string file, int& num, int& numedge, vector<string>& namelist);  

void Char2String(char& c, string& s);
string Int2String(int i);


#endif
