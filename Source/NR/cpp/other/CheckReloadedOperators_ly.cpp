/*************************************************************************
 *    Physics 598CPA Fall 2004                                           *
 *    Problem Set #8: PDE, 2D hydrodynamics                              *
 *    Author: Yang Liu (yangliu@uiuc.edu)                                *
 *    Date: 10/22/04                                                     *
 *************************************************************************/
#include <string>
#include <iostream>
#include <fstream>
#include <iomanip>
#include "nr.h"
using namespace std;

int main()
{

// check reloaded operator of NRMat3d_DP 10/23/04
    
    Mat3D_DP O;
    cout << "O= " << endl << O << endl;

    Mat3D_DP A(2,3,4);  
    cout << "A= " << endl << A << endl;
 
    Mat3D_DP B(2,3,4);
    B = 2.5;
    cout << "B= " << endl << B << endl;

    Mat3D_DP C(2,3,4);
    C = B;
    cout << "C= " << endl << C << endl;
  
    double d1[24] = {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24}; 
    Mat3D_DP D1(d1,2,3,4);
    cout << "D1= " << endl << D1 << endl;

    double d2[2][3][4] = { {{1,2,3,4},   {5,6,7,8},     {9,10,11,12}},
                           {{13,14,15,16},{17,18,19,20},{21,22,23,24}}}; 
   
//    Mat3D_DP D2(d2,2,3,4); //bad
   Mat3D_DP D2(**d2,2,3,4);// OK
//    Mat3D_DP D2(*d2[0],2,3,4); // OK
    cout << "D2= " << endl << D2 << endl;


/*
// check reloaded operator of NRMat_DP 10/23/04
    
    Mat_DP O;
    cout << "O= " << endl << O << endl;

    Mat_DP A(4,5);  
    cout << "A= " << endl << A << endl;
 
    Mat_DP B(2.0,4,5);
    cout << "B= " << endl << B << endl;

    const DP a = 1.0; 
    Mat_DP C(a,4,5);
    cout << "C= " << endl << C << endl;

    double d1[20] = {11,12,13,14,15,21,22,23,24,25,31,32,33,34,35,41,42,43,44,45}; 
    Mat_DP D1(d1,4,5);
    cout << "D1= " << endl << D1 << endl;

    double d2[4][5] = {{11,12,13,14,15},{21,22,23,24,25},{31,32,33,34,35},{41,42,43,44,45}}; 
//    Mat_DP D2(d2,4,5); //bad
//   Mat_DP D2(d2[0],4,5);// OK
    Mat_DP D2(*d2,4,5); // OK


    cout << "D2= " << endl << D2 << endl;


    Mat_DP E = D2;
    cout << "E= " << endl << E << endl;

    Mat_DP F(4,5);
    F = 2.0;
    cout << "F= " << endl << F << endl;

    F += E; 
    cout << "F += E " << endl << F << endl;

    F = 2.0;
    F -= E; 
    cout << "F -= E " << endl << F << endl;

    F = 2.0;
    F += a; 
    cout << "F += a " << endl << F << endl;

    F = 2.0;
    F -= a; 
    cout << "F -= a " << endl << F << endl;

    F =2.0;
    Mat_DP G;
    G = -F;
    cout << "G = -F " << endl << G << endl;

    G = +F;
    cout << "G = + F " << endl << G << endl;

    cout << "F= " << endl << F << endl;
    cout << "G= " << endl << G << endl;

    cout << "F + G " << endl << F+G << endl;
    cout << "F - G " << endl << F-G << endl;

    cout << "F + 3 " << endl << F+3.0 << endl;
    cout << "F - 3 " << endl << F-3.0 << endl;


    cout << " 1.4*F " << endl << 1.4*F << endl;
    cout << " F*1.5 " << endl << F*1.5 << endl;
    cout << " F/0.5 " << endl << F/0.5 << endl;
    cout << " F*G " << endl << F*G << endl;

    cout << "F= " << endl << F << endl;
    cout << "G= " << endl << G << endl;
    cout << "F + 3*0.7*G " << endl << F + 3.0*0.7*G << endl;
    
    F += 2*0.3*G;
    cout << "F + 2*0.3*G " << endl << F << endl;

    F =2.0;
    Vec_DP V(2.0,5);
    cout << "V= " << V << endl;
    cout << " F * V " << endl << F*V << endl;

*/

/* check reloaded operator of NRVec_DP 10/23/04
   
    Vec_DP O;
    cout << "O= " << O << endl;

    Vec_DP A(5);
    cout << "A= " << A << endl;

    Vec_DP B(2.0,5);
    cout << "B= " << B << endl;

    const DP a = 1.0;
    Vec_DP C(a,5);
    cout << "C= " << C << endl;

    double b[5] = {1,2,3,4,5}; 
    Vec_DP D(b,5);
    cout << "D= " << D << endl;

    Vec_DP E = D;
    cout << "E= " << E << endl;

    Vec_DP F(5);
    F = 2.0;
    cout << "F= " << F << endl;

    F += E; 
    cout << "F += E " << F << endl;

    F = 2.0;
    F -= E; 
    cout << "F -= E " << F << endl;

    F = 2.0;
    F += a; 
    cout << "F += a " << F << endl;

    F = 2.0;
    F -= a; 
    cout << "F -= a " << F << endl;

    F =2.0;
    Vec_DP G;
    G = -F;
    cout << "G = -F " << G << endl;

    G = +F;
    cout << "G = + F " << G << endl;

    cout << "F= " << F << endl;
    cout << "G = " << G << endl;

    cout << "F + G " << F+G << endl;
    cout << "F - G " << F-G << endl;

    cout << "F + 3 " << F+3.0 << endl;
    cout << "F - 3 " << F-3.0 << endl;

    cout << " 1.4*F " << 1.4*F << endl;
    cout << " F*1.5 " << F*1.5 << endl;
    cout << " F/0.5 " << F/0.5 << endl;
    cout << " F*G " << F*G << endl;
    cout << " F dor G " << dot(F,G) << endl;

*/
}





