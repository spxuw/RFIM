#ifndef LINESEGMENT_H_
#define LINESEGMENT_H_
#include <iostream>
using namespace std;



class LineSegment                
{
 public:
    
    LineSegment(){}

    LineSegment(double b, double B, double e0, double e, int mag, int num, int Index)
	{Hlower=b; Hupper=B; E0= e0; E= e; M = mag; N=num; index=Index;}

    LineSegment(State& C, int Index)
	{E0= C.GetE0(); E=C.GetE(); M= C.GetM(); N=C.GetN(); index=Index;}

    
    bool operator() (const LineSegment& A, const LineSegment& B) const 
	{return A.Get_Hlower() < B.Get_Hlower() ;}

    
    friend double Get_CrossField(LineSegment& A, LineSegment& B)  
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
    double Hlower;               
    double Hupper;               
    double E0;                   
    double E;                    
    int    M;                    
    int    N;                    
    int    index;                
};

#endif
