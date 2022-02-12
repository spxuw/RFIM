#ifndef _CROSSINGFIELD_H_
#define _CROSSINGFIELD_H_
#include <cmath>
#include <iostream>
#include <queue>
using namespace std;

// To implement Vives' algorithm to calculate the evolution of the Ground States of the RFIM
// as sweeping the external field, we introduce the following class. Y.L. 06/06/06
// For Algorithm details, see Vives' paper JCP 160, 117-125 (2000)

class CrossingField
{
 public:
    CrossingField(double hx, int i1, int i2) {Hx=hx;k=i1;l=i2;};// constructor
    double GetHx() const {return Hx;}
    int    Getk()  const {return k;}
    int    Getl()  const {return l;}

 private:
    double  Hx;   // the crossing field
    int     k;    // index of the ground state
    int     l;    // index of the ground state

};


#endif // _CROSSINGFIELD_H_
