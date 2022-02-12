#ifndef LINESEGMENT_H_
#define LINESEGMENT_H_
#include <iostream>
using namespace std;


// line segment in the M-H curve, representing a ground state and its H-boundaries
class LineSegment                
{
 public:
    // constructor
    LineSegment(){}

    LineSegment(double b, double B, double e0, double e, int mag, int num, int Index)
	{Hlower=b; Hupper=B; E0= e0; E= e; M = mag; N=num; index=Index;}

    LineSegment(State& C, int Index)
	{E0= C.GetE0(); E=C.GetE(); M= C.GetM(); N=C.GetN(); index=Index;}

    // for the sorting of the Hlower field
    bool operator() (const LineSegment& A, const LineSegment& B) const 
	{return A.Get_Hlower() < B.Get_Hlower() ;}

    // get the crossing field of two states
    friend double Get_CrossField(LineSegment& A, LineSegment& B)  // Attention: the symbol & is necessary. 
	{ return (A.GetE0() - B.GetE0())/(double)(A.GetM() - B.GetM());}

    double Get_Hlower() const {return Hlower;}        
    double Get_Hupper() const {return Hupper;}        
    double GetE0()    const {return E0;}
    int    GetM()     const {return M;}
    void   CalE(double H) {E=E0 - H*M;}
    double GetE()     const {return E;}
    int    Get_index()     const {return index;}

    void   Set_Hlower(double b) {Hlower=b;}
    void   Set_Hupper(double B) {Hupper=B;}
    void   Set_Hboundaries(double b, double B) {Hlower=b;Hupper=B;}

    double  Getm()     const {return M/(double)N;}
    double  Get_nup()  const {return 0.5*(1+M/(double)N);}
    double  Get_ndn()  const {return 0.5*(1-M/(double)N);}

 private:
    double Hlower;               // the lower boundary of the validity of the GS
    double Hupper;               // the upper boundary of the validity of the GS
    double E0;                   // internal energy E0= - sum_<i,j> Si*Sj - sum_i h[i]*Si 
    double E;                    // total energy    E = E0 - Hext * M
    int    M;                    // magnetization   M = sum_i Si = nUP-nDN = 2*nUP-N = N-2*nDN
    int    N;                    // number of spins
    int    index;                // index of this LineSegment  
};

#endif /* ! LINESEGMENT_H_ */
