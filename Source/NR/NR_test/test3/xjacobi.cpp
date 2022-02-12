#include <iostream>
#include <iomanip>
#include "nr.h"
using namespace std;

// Driver for routine jacobi


int main(void)
{
        DP a_d[3*3]=
          {1.0,2.0,3.0,
          2.0,2.0,3.0,
          3.0,3.0,3.0};

	DP b_d[3*3]=
          {2.0,-1.0,-1.0,
           -1.0,2.0,-1.0,
           -1.0,-1.0,2.0};

        int i,j,k,l,kk,ll,nrot;
        Mat_DP a(a_d,3,3);
	Mat_DP b(b_d,3,3);

        cout << fixed << setprecision(6);
          Vec_DP d(3),r(3);
          Mat_DP v(3,3);
          //NR::jacobi(a,d,v,nrot);
          NR::jacobi(b,d,v,nrot);

          cout << "number of JACOBI rotations: " << nrot << endl;
          cout << "eigenvalues: " << endl;
          for (j=0;j<3;j++) {
            cout << setw(12) << d[j];
            if ((j+1) % 5 == 0) cout << endl;
          }

          cout << endl << "eigenvectors:" << endl;
          for (j=0;j<3;j++) {
            cout << setw(9) << "number" << setw(4) << (j+1) << endl;
            for (k=0;k<3;k++) {
              cout << setw(12) << v[k][j];
              if ((k+1) % 5 == 0) cout << endl;
            }
            cout << endl;
          }

          cout << endl << "reorder eigenvalues and eigenvectors:" << endl;
	  NR::eigsrt(d,v);
          cout << "eigenvalues: " << endl;
          for (j=0;j<3;j++) {
            cout << setw(12) << d[j];
            if ((j+1) % 5 == 0) cout << endl;
          }

          cout << endl << "eigenvectors:" << endl;
          for (j=0;j<3;j++) {
            cout << setw(9) << "number" << setw(4) << (j+1) << endl;
            for (k=0;k<3;k++) {
              cout << setw(12) << v[k][j];
              if ((k+1) % 5 == 0) cout << endl;
            }
            cout << endl;
          }



        
        return 0;
}



