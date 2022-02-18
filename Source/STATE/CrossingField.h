#ifndef _CROSSINGFIELD_H_
#define _CROSSINGFIELD_H_
#include <cmath>
#include <iostream>
#include <queue>
using namespace std;

class CrossingField
{
 public:
    CrossingField(double hx, int i1, int i2) {Hx=hx;k=i1;l=i2;};
    double GetHx() const {return Hx;}
    int    Getk()  const {return k;}
    int    Getl()  const {return l;}

 private:
    double  Hx;   
    int     k;    
    int     l;    

};


#endif 
