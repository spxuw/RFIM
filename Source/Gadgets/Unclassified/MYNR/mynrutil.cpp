// Note that here, we are using offset pointer:
// for example, 
//     v=(int *)calloc((unsigned) (nh-nl+1),sizeof(int));
//     v-=nl;
// if we set nl=0, nh=4, then as usual
//      v[0] means the first element
//      v[4] means the last element 
// if we set nl=1, nh=5
// then v[1] actually means the first element
//      v[5] ................   last element
// 
// Note that if we use offset pointer, when we do the "free" operation,
// we have to recover the pointer, and then do the "free"!!!
// Y.L. 01/07/09

/*
**      File:   nrutil.c
**      Purpose: Memory allocation routines borrowed from the
**		book "Numerical Recipes" by Press, Flannery, Teukolsky,
**		and Vetterling. 
*/

#include <iostream>
#include <string>
#include "mynrutil.h"

using namespace std;
const double DBL_MAX = numeric_limits<double>::max();
const double DBL_MIN = numeric_limits<double>::min();


inline void mynrerror(const string error_text)
{
  cout << "Numerical Recipes run-time error..." << endl;
  cout << error_text << endl;
  cout << "...now exiting to system..." << endl;
  exit(1);
}


///////////////////////////////////////////////////////////////////////
int* ivector(int nl, int nh)
{
  int *v;
  v = new int [nh-nl+1];
  if (!v) mynrerror("allocation failure in ivector()");
  return v-nl;
}

double* dvector(int nl, int nh)
{
  double* v = new double [nh-nl+1];
  if (!v) mynrerror("allocation failure in dvector()");
  return v-nl;
}

void free_ivector(int* v, int nl, int nh)
{
  v += nl;
  delete [] v;
}

void free_dvector(double* v, int nl, int nh)
{
  v += nl;
  delete [] v;
}
///////////////////////////////////////////////////////////////////////


///////////////////////////////////////////////////////////////////////
int** imatrix(int nrl, int nrh, int ncl, int nch)
{
  int** m = new int* [nrh-nrl+1];
  if (!m) mynrerror("allocation failure 1 in matrix()");
  m -= nrl;

  for(int i=nrl;i<=nrh;i++) {
    m[i] = new int [nch-ncl+1];
    if (!m[i]) mynrerror("allocation failure 2 in matrix()");
    m[i] -= ncl;
  }
  return m;
}


double** dmatrix(int nrl, int nrh, int ncl, int nch)
{
  double** m = new double* [nrh-nrl+1];
  if (!m) mynrerror("allocation failure 1 in matrix()");
  m -= nrl;

  for(int i=nrl;i<=nrh;i++) {
    m[i] = new double [nch-ncl+1];
    if (!m[i]) mynrerror("allocation failure 2 in matrix()");
    m[i] -= ncl;
  }
  return m;
}

void free_imatrix(int** m, int nrl,int nrh, int ncl, int nch)
{
  for(int i=nrh;i>=nrl;i--) 
    {
      m[i] += ncl;
      delete [] m[i];
    }

  m += nrl;
  delete [] m;
}

void free_dmatrix(double** m, int nrl,int nrh, int ncl, int nch)
{
  for(int i=nrh;i>=nrl;i--) 
    {
      m[i] += ncl;
      delete [] m[i];
    }

  m += nrl;
  delete [] m;
}
///////////////////////////////////////////////////////////////////////



///////////////////////////////////////////////////////////////////////
int*** imatrix(int xi, int xf, int yi, int yf, int zi, int zf)
{
  int*** m = new int** [xf-xi+1];
  if (!m) mynrerror("allocation failure 1 in matrix()");
  m -= xi;

  for(int x=xi;x<=xf;x++) 
    {
      m[x] = new int* [yf-yi+1];
      if (!m[x]) mynrerror("allocation failure 2 in matrix()");
      m[x] -= yi;

      for(int y=yi;y<=yf;y++) 
	{
	  m[x][y] = new int [zf-zi+1];
	  if (!m[x][y]) mynrerror("allocation failure 3 in matrix()");
	  m[x][y] -= zi;
	}
    }
  return m;
}


double*** dmatrix(int xi, int xf, int yi, int yf, int zi, int zf)
{
  double*** m = new double** [xf-xi+1];
  if (!m) mynrerror("allocation failure 1 in matrix()");
  m -= xi;

  for(int x=xi;x<=xf;x++) 
    {
      m[x] = new double* [yf-yi+1];
      if (!m[x]) mynrerror("allocation failure 2 in matrix()");
      m[x] -= yi;

      for(int y=yi;y<=yf;y++) 
	{
	  m[x][y] = new double [zf-zi+1];
	  if (!m[x][y]) mynrerror("allocation failure 3 in matrix()");
	  m[x][y] -= zi;
	}
    }
  return m;
}


void free_imatrix(int*** m, int xi, int xf, int yi, int yf, int zi, int zf)
{
  for(int x=xf; x>=xi; x--) 
    {
      for(int y=yf; y>=yi; y--) 
	{
	  m[x][y] += zi;
	  delete [] m[x][y];
	}
      m[x] += yi;
      delete [] m[x];
    }
  m += xi;
  delete [] m;
}

void free_dmatrix(double*** m, int xi, int xf, int yi, int yf, int zi, int zf)
{
  for(int x=xf; x>=xi; x--) 
    {
      for(int y=yf; y>=yi; y--) 
	{
	  m[x][y] += zi;
	  delete [] m[x][y];
	}
      m[x] += yi;
      delete [] m[x];
    }
  m += xi;
  delete [] m;
}
///////////////////////////////////////////////////////////////////////



//////////////////////////////////////////////////////////////////////
int** imatrix(int M, int* T)
{
  int** X = new int* [M];
  if (!X) mynrerror("allocation failure 1 in matrix()");
  X--;

  for(int m=1;m<=M;m++) {
    X[m] = new int [T[m]];
    if (!X[m]) mynrerror("allocation failure 2 in matrix()");
    X[m]--;
  }
  return X;
}


double** dmatrix(int M, int* T)
{
  double** X = new double* [M];
  if (!X) mynrerror("allocation failure 1 in matrix()");
  X--;

  for(int m=1;m<=M;m++) {
    X[m] = new double [T[m]];
    if (!X[m]) mynrerror("allocation failure 2 in matrix()");
    X[m]--;
  }
  return X;
}


void free_imatrix(int** X, int M)
{
  for (int m=1; m<=M; m++) {
    X[m]++;
    delete [] X[m];
  }

  X++;
  delete [] X;
}


void free_dmatrix(double** X, int M)
{
  for (int m=1; m<=M; m++) {
    X[m]++;
    delete [] X[m];
  }

  X++;
  delete [] X;
}
//////////////////////////////////////////////////////////////////////



//////////////////////////////////////////////////////////////////////
int*** imatrix(int M, int* T, int d)
{
  int*** X = new int** [M];
  if (!X) mynrerror("allocation failure 1 in matrix()");
  X--;

  for(int m=1;m<=M;m++) {
    X[m] = new int* [T[m]];
    if (!X[m]) mynrerror("allocation failure 2 in matrix()");
    X[m]--;

    for(int t=1;t<=T[m];t++)
      {
	X[m][t] = new int [d];
	if (!X[m][t]) mynrerror("allocation failure 3 in matrix()");
	X[m][t]--;
      }
  }
    
  return X;
}


double*** dmatrix(int M, int* T, int N)
{
  double*** X = new double** [M];
  if (!X) mynrerror("allocation failure 1 in matrix()");
  X--;

  for(int m=1;m<=M;m++) {
    X[m] = new double* [T[m]];
    if (!X[m]) mynrerror("allocation failure 2 in matrix()");
    X[m]--;

    for(int t=1;t<=T[m];t++)
      {
	X[m][t] = new double [N];
	if (!X[m][t]) mynrerror("allocation failure 3 in matrix()");
	X[m][t]--;
      }
  }
    
  return X;
}


void free_imatrix(int*** X, int M, int* T)
{
  for (int m=1; m<=M; m++) {
    for (int t=1; t<=T[m]; t++) {
      X[m][t]++;
      delete [] X[m][t];
    }
    X[m]++;
    delete [] X[m];
  }
    
  X++;
  delete [] X;
}


void free_dmatrix(double*** X, int M, int* T)
{
  for (int m=1; m<=M; m++) {
    for (int t=1; t<=T[m]; t++) {
      X[m][t]++;
      delete [] X[m][t];
    }
    X[m]++;
    delete [] X[m];
  }
    
  X++;
  delete [] X;
}
//////////////////////////////////////////////////////////////////////




/////////////////////////////////////////////////////////////////
/* generate a symmtric and positive-definite matrix randomly. */
void GenerateCovarianceMatrix(double**& C, int d, int seed)
{
  static Rand rand;
  rand.seed(seed);
    
  double** M = dmatrix(1,d,1,d);
  for(int i=1; i<=d; i++)
    for(int j=1; j<=d; j++)
      M[i][j] = rand.gaussian(1.0);

  for(int i=1; i<=d; i++) {
    for(int j=1; j<=d; j++) {
      C[i][j] = 0;
      for(int k=1;k<=d;k++)
	C[i][j] += M[i][k]*M[j][k];
    }
  }
    
  free_dmatrix(M,1,d,1,d);
}
/////////////////////////////////////////////////////////////////



/////////////////////////////////////////////////////////////////
void CholeskyDecomposition(double**& C, double**& L, int d)
{
  double** A = dmatrix(1,d,1,d);
  MatrixCopy(C,A,d);

  double* p = dvector(1,d);
    
  int i,j,k;
  double sum;
  for (i=1;i<=d;i++) {
    for (j=i;j<=d;j++) {
      for (sum=A[i][j],k=i-1;k>=1;k--) 
	sum -= A[i][k]*A[j][k];
      if (i == j) {
	if (sum <= 0.0)
	  mynrerror("choldc failed");
	p[i]=sqrt(sum);
      } 
      else A[j][i]=sum/p[i];
    }
  }

  for(i=1; i<=d; i++)
    {
      L[i][i] = p[i];
      for(j=1; j<i; j++)
	{
	  L[i][j] = A[i][j];
	  L[j][i] = 0;
	}
    }
}
/////////////////////////////////////////////////////////////////



///////////////////////////////////////////////////
void PrintMatrix(double**& A, int d)
{
  for(int i=1;i<=d;i++)
    {
      for(int j=1;j<=d;j++) {
	cout.width(8);
	cout << A[i][j] << ' ';
      }
      cout << endl;
    }
  cout << endl;
}
///////////////////////////////////////////////////


///////////////////////////////////////////////////
// copy A(dxd) to B(dxd)
void MatrixCopy(double**& A, double**& B, int d)
{
  for(int i=1;i<=d;i++)
    for(int j=1;j<=d;j++)
      B[i][j] = A[i][j];
}

// copy A(nxm) to B(nxm)
void MatrixCopy(double**& A, double**& B, int n, int m)
{
  for(int i=1;i<=n;i++)
    for(int j=1;j<=m;j++)
      B[i][j] = A[i][j];
}

// copy A(nxm) to B(nxm)
void MatrixCopy(int**& A, int**& B, int n, int m)
{
  for(int i=1;i<=n;i++)
    for(int j=1;j<=m;j++)
      B[i][j] = A[i][j];
}
///////////////////////////////////////////////////


///////////////////////////////////////////////////
// transpose A(dxd) to B(dxd)
void MatrixTranspose(double**& A, double**& B, int d)
{
  for(int i=1;i<=d;i++)
    for(int j=1;j<=d;j++)
      B[i][j] = A[j][i];
}
///////////////////////////////////////////////////


///////////////////////////////////////////////////
// multiply A(dxd) by B(dxd), save result to C
void MxM(double**& A, double**& B, double**& C, int d)
{
    
  for(int i=1;i<=d;i++){
    for(int j=1;j<=d;j++){
      C[i][j] = 0;
      for(int k=1;k<=d;k++)
	C[i][j] += A[i][k]*B[k][j];
    }
  }
}
///////////////////////////////////////////////////


///////////////////////////////////////////////////
// multiply A(nxm) by B(mxq), save result to C(nxq)
void MxM(double**& A, double**& B, double**& C, int n, int m, int q)
{
    
  for(int i=1;i<=n;i++){
    for(int j=1;j<=q;j++){
      C[i][j] = 0;
      for(int k=1;k<=m;k++)
	C[i][j] += A[i][k]*B[k][j];
    }
  }
}
///////////////////////////////////////////////////



///////////////////////////////////////////////////
// multiply A(dxd) by B(dx1), save result to C(dx1)
void MxV(double**& A, double*& B, double*& C, int d)
{
  for(int i=1;i<=d;i++){
    C[i] = 0;
    for(int j=1;j<=d;j++){
      C[i] += A[i][j]*B[j];
    }
  }
}
///////////////////////////////////////////////////


///////////////////////////////////////////////////
//  A(dx1) + B(dx1) = C(dx1)
void VaV(double*& A, double*& B, double*& C, int d)
{
  for(int i=1;i<=d;i++)
    C[i] = A[i] + B[i];
}
///////////////////////////////////////////////////


///////////////////////////////////////////////////
//  A(dx1) - B(dx1) = C(dx1)
void VsV(double*& A, double*& B, double*& C, int d)
{
  for(int i=1;i<=d;i++)
    C[i] = A[i] - B[i];
}
///////////////////////////////////////////////////


///////////////////////////////////////////////////
//  A(dx1) . B(dx1) = C
double VdotV(double*& A, double*& B, int d)
{
  double sum=0;
  for(int i=1;i<=d;i++)
    sum += A[i]*B[i];

  return sum;
}
///////////////////////////////////////////////////



//////////////////////////////////////////////////////////////////////
//LU Decomposition:
// Given a matrix a(dxd), this routine replaces it by the LU decomposition 
// of a rowwise permutation of itself. a is input. On output, it is arranged.
// index(dx1) is an output vector that records the row permutation effected 
// by the partial pivoting. s is output as +/-1 depending on whether the 
// number of row interchanges was even or odd, respectively.
void LUD(double**& a, int*& indx, double& s, int d)
{
  const double TINY=1.0e-20;
  int i,j,k;
  int imax=0;
  double big,dum,sum,temp;

  double* vv = dvector(1,d);
  s=1.0;
  for (i=1;i<=d;i++) {
    big=0.0;
    for (j=1;j<=d;j++){
      if ((temp=fabs(a[i][j])) > big) 
	big=temp;
    }
    if (big == 0.0) 
      {
	cout << "Singular matrix in routine LUD\n";
	exit(0);
      }
    vv[i]=1.0/big;
  }

  for (j=1;j<=d;j++) {
    for (i=1;i<j;i++) {
      sum=a[i][j];
      for (k=1;k<i;k++) 
	sum -= a[i][k]*a[k][j];
      a[i][j]=sum;
    }
    big=0.0;
    for (i=j;i<=d;i++) {
      sum=a[i][j];
      for (k=1;k<j;k++) 
	sum -= a[i][k]*a[k][j];
      a[i][j]=sum;
      if ((dum=vv[i]*fabs(sum)) >= big) {
	big=dum;
	imax=i;
      }
    }
    if (j != imax) {
      for (k=1;k<=d;k++) {
	dum=a[imax][k];
	a[imax][k]=a[j][k];
	a[j][k]=dum;
      }
      s = -s;
      vv[imax]=vv[j];
    }
    indx[j]=imax;
    if (a[j][j] == 0.0) a[j][j]=TINY;
    if (j != d) {
      dum=1.0/(a[j][j]);
      for (i=j+1;i<=d;i++) 
	a[i][j] *= dum;
    }
  }

  free_dvector(vv,1,d);
}
//////////////////////////////////////////////////////////////////////



//////////////////////////////////////////////////////////////////////
// Compose separately the lower and upper matrices
void Separate(double**& a, double**& xl, double**& xu, int d)
{
  int k,l;
  for (k=1;k<=d;k++) {
    for (l=1;l<=d;l++) {
      if (l > k) {
	xu[k][l]=a[k][l];
	xl[k][l]=0.0;
      } else if (l < k) {
	xu[k][l]=0.0;
	xl[k][l]=a[k][l];
      } else {
	xu[k][l]=a[k][l];
	xl[k][l]=1.0;
      }
    }
  }
}
//////////////////////////////////////////////////////////////////////



//////////////////////////////////////////////////////////////////////
// solves the set of d linear equations A*X=B
// Here, a(dxd) is input, not as the matrix A but rather as its LU decomposition,
// determined by the above routine. indx(dx1) is input as the permutation vector 
// returned by LUD. b(dx1) is input as the RHS vector B, and returns with the 
// solution vector X. a and indx are not modified by this routine and can be left
// in place for successive calls with different RHS b. This routine takes into account
// the possibility that b will begin with many zero elements, so it is efficient 
// for use in matrix inversion.
void LUbksb(double**& a, int*& indx, double*& b, int d)
{
  int i,ii=1,ip,j;
  double sum;

  for (i=1;i<=d;i++) {
    ip=indx[i];
    sum=b[ip];
    b[ip]=b[i];
    if (ii != 1){
      for (j=ii-1;j<i;j++) 
	sum -= a[i][j]*b[j];
    }
    else if (sum != 0.0)
      ii=i+1;
    b[i]=sum;
  }
  for (i=d;i>=1;i--) {
    sum=b[i];
    for (j=i+1;j<=d;j++) 
      sum -= a[i][j]*b[j];
    b[i]=sum/a[i][i];
  }
}
//////////////////////////////////////////////////////////////////////



//////////////////////////////////////////////////////////////////////
// Use the above LU decomposition and backsubstitution routines, it is 
// completely straightforward to find the inverse of a matrix column by column.
// Also the determinant of matrix a will be retured.
void MatrixInverse(double**& a, double**&y, double& det, int d)
{
  // we don't want to destory matrix a, so we make a copy first
  double** a0 = dmatrix(1,d,1,d);
  MatrixCopy(a,a0,d);

  double* col = dvector(1,d);
  int* indx   = ivector(1,d);
  double s;

  LUD(a0,indx,s,d); // note that here a0 will be modified.

  int i,j;
  for(j=1;j<=d;j++)
    s *= a0[j][j];
  det = s; 

  for(j=1;j<=d;j++){
    for(i=1;i<=d;i++)
      col[i]=0.0;
    col[j]=1.0;
    LUbksb(a0, indx, col, d);
    for(i=1;i<=d;i++)
      y[i][j]=col[i];
  }

  free_dmatrix(a0,1,d,1,d);
  free_dvector(col,1,d);
  free_ivector(indx,1,d);
}
//////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////
// Generate permutation matrix P
void PermutationMatrix(int**& P, int n)
{
  for(int i=1; i<=n; i++) 
    for(int j=1; j<=n; j++)
      P[i][j] = 0;

  vector<int> x(n,0);
  for(int i=0; i<n; i++) 
    x[i] = i+1;
  random_shuffle(x.begin(), x.end());

  for(int i=1; i<=n; i++) 
    P[i][x[i-1]] = 1;
}
//////////////////////////////////////////////////////////////////////


/////////////////////////////////////////////////////////////////////
// multiply A(axb) by B(bxc), save result to C(axc)
void MxM(int**& A, int**& B, int**& C, int a, int b, int c)
{
    
  for(int i=1;i<=a;i++){
    for(int j=1;j<=c;j++){
      C[i][j] = 0;
      for(int k=1;k<=b;k++)
	C[i][j] += A[i][k]*B[k][j];
    }
  }
}
/////////////////////////////////////////////////////////////////////



///////////////////////////////////////////////////
// multiply P1(axa) * Q(axb) * P2(bxb) with P1 P2 two permutation matries, save result to Y(axb)
void P1QP2(int**& Q, int**& Y, int a, int b)
{
  int** P1 = imatrix(1,a, 1,a);
  PermutationMatrix(P1, a);
  //MatrixPrint(P1, a, a);
  int** P2 = imatrix(1,b, 1,b);
  PermutationMatrix(P2, b);
  //MatrixPrint(P2, b, b);

  int** P1Q = imatrix(1,a, 1,b);
  MxM(P1, Q, P1Q, a, a, b);
  MxM(P1Q, P2, Y, a, b, b);

  free_imatrix(P1Q,1,a,1,b);
  free_imatrix(P1,1,a,1,a);
  free_imatrix(P2,1,b,1,b);
}
///////////////////////////////////////////////////



///////////////////////////////////////////////////
// multiply P1(axa) * Q0(axb) * P2(bxb) with P1 P2 two permutation matries, save result to Q(axb) and also keep tracking the permutation orders 
void P1QP2(int**& Q0, int**& Q, int**& rowindex, int**& colindex, int a, int b)
{
  // copy Q0 to Q
  MatrixCopy(Q0,Q,a,b);

  int** rowindex0 = imatrix(1,a, 1,b);
  MatrixCopy(rowindex,rowindex0,a,b);

  int** colindex0 = imatrix(1,a, 1,b);
  MatrixCopy(colindex,colindex0,a,b);

  // P1 . Q0 ==> permutate rows of Q0
  vector<int> x(a,0);
  for(int i=0; i<a; i++) 
    x[i] = i+1;
  random_shuffle(x.begin(), x.end());

  for(int i=1; i<=a; i++) {
    for(int j=1; j<=b; j++) { 
      Q[i][j] = Q0[x[i-1]][j];
      rowindex[i][j] = rowindex0[x[i-1]][j];
    }
  }

  int** Q1 = imatrix(1,a,1,b);
  MatrixCopy(Q,Q1,a,b);

  // Q . P2 ==> permutate cols of Q
  vector<int> y(b,0);
  for(int j=0; j<b; j++) 
    y[j] = j+1;
  random_shuffle(y.begin(), y.end());

  for(int j=1; j<=b; j++) {
    for(int i=1; i<=a; i++) { 
      Q[i][j] = Q1[i][y[j-1]];
      colindex[i][j] = colindex0[i][y[j-1]];
    }
  }  

  free_imatrix(Q1,1,a,1,b);
  free_imatrix(rowindex0,1,a,1,b);
  free_imatrix(colindex0,1,a,1,b);

}
///////////////////////////////////////////////////



/////////////////////////////////////////////////////////
// permute row-i1 and row-i2 for the matrix Q(axb)
void Permute_Row(int**& Q, int**& rowindex, int a, int b, int i1, int i2)
{
  for(int j=1;j<=b;j++) {
    swap(Q[i1][j], Q[i2][j]);
    swap(rowindex[i1][j], rowindex[i2][j]);
  }
}
///////////////////////////////////////////////////


/////////////////////////////////////////////////////////
// permute col-j1 and col-j2 for the matrix Q(axb)
void Permute_Col(int**& Q, int**& colindex, int a, int b, int j1, int j2)
{
  for(int i=1;i<=a;i++) {
    swap(Q[i][j1], Q[i][j2]);
    swap(colindex[i][j1], colindex[i][j2]);
  }
}
///////////////////////////////////////////////////



///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/* Check whether a Matrix Q is in the Form-III, i.e. there exists two permutation matries P1, P2 such that 

             [ @ ... @ x 0 ..........0 ]
             | .     . @ x 0 .       . |
             | .     . . @ x 0 .     . |
   P1 Q P2 = | .     . .   @ x 0 .   . |
             | .     . .     @ x 0 . . |
             | .     . .       @ x 0 . |
             | .     . .         @ x 0 |
             [ @ ... @ @ ..........@ x ]

   where the x-elements must be non-zero. The @-elements can be either zero or non-zero. 

   This is used to determine whether a structured system [A,B] is strongly structural controllable.
   For details, see the paper entitled "On strong structural controllability of linear systems" by K.J. Reinschke et al.
*/
 
void MatrixFormIII_check_test(int n, int m)
{
  Rand rand;
  rand.seed(1);

  /*//make a matrix in the exact form-III
               1     m m+1           m+n
             [ @ ... @ x 0 ..........0 ]1
             | .     . @ x 0 .       . |2
             | .     . . @ x 0 .     . |
   i.e. Q0 = | .     . .   @ x 0 .   . |
             | .     . .     @ x 0 . . |
             | .     . .       @ x 0 . |
             | .     . .         @ x 0 |
             [ @ ... @ @ ..........@ x ]n
   */
  /*
  int** Q0 = imatrix(1,n,1,n+m);

  for(int j=1; j<=m; j++) 
    for(int i=1; i<=n; i++) 
      Q0[i][j] = (rand.ran1()>0.5) ? 1 : 0;

  for(int j=m+1; j<=n+m; j++) { 
    for(int i=1; i<=n; i++){
      if(j-i==m)
	Q0[i][j] = 2;
      else if(j-i>m)
	Q0[i][j] = 0;
      else
	Q0[i][j] = (rand.ran1()>0.5) ? 1 : 0;
    }
  }
  //Q0[1][n+m]=5;
  //Q0[2][2]=0;
  //Q0[n][2]=0;
  */
  


   // make a matrix in such a form of [A,B] that it can be transformed into form-III
   /*
               1 2           n n+1   n+m
             [ x 0 ..........0 @ ... @ ]1
             | @ x 0 .       . @ ... @ |2
             | . @ x 0 .     . @ ... @ |
   i.e. Q0 = | .   @ x 0 .   . @ ... @ |
             | .     @ x 0 . . @ ... @ |
             | .       @ x 0 . @ ... @ |
             | .         @ x 0 @ ... @ |
             [ @ ..........@ x @ ... @ ]n
   */
  
  int** Q0 = imatrix(1,n,1,n+m);

  for(int j=1; j<=n; j++) {
    for(int i=1; i<=n; i++) {
      if(j-i==0)
	Q0[i][j] = 10+i;
      else if(j-i>0)
	Q0[i][j] = 0;
      else
	Q0[i][j] = (rand.ran1()>0.5) ? 1 : 0;// 1;
    }
  }

  for(int j=n+1; j<=n+m; j++) 
    for(int i=1; i<=n; i++)
      Q0[i][j] = (rand.ran1()>0.5) ? 1 : 0; // 1;
  


  /* test
  //This form of Q0 can be used to test the tracking of row, col index after the permutations 
  for(int i=1; i<=n; i++) 
    for(int j=1; j<=n+m; j++) 
      Q0[i][j] = 10*i + j;
  */

  cout << " Q0="; MatrixPrint(Q0, n, n+m);

  //store the original row index, col index of each element in Q
  int** rowindex = imatrix(1,n,1,n+m);
  int** colindex = imatrix(1,n,1,n+m);
  for(int i=1; i<=n; i++) { 
    for(int j=1; j<=n+m; j++) {
      rowindex[i][j] = i;
      colindex[i][j] = j;
    }
  } 
  
  cout << "Original row-col index:\n";
  for(int i=1; i<=n; i++) { 
    for(int j=1; j<=n+m; j++) {
      cout << rowindex[i][j]  << colindex[i][j] << " ";
    }
    cout << endl; 
  }

  // permutate the matrix Q0, using two permutation matries P1 (nxn) and P2 ((n+m)x(n+m))
  int** Q = imatrix(1,n, 1,n+m);
  P1QP2(Q0, Q, rowindex, colindex, n, n+m); // Q = P1 Q0 P2
  //MatrixCopy(Q0,Q,n,n+m);
  cout << " Q="; MatrixPrint(Q, n, n+m);

  cout << "row-col index:\n";
  for(int i=1; i<=n; i++) { 
    for(int j=1; j<=n+m; j++) {
      cout << rowindex[i][j]  << colindex[i][j] << " ";
    }
    cout << endl; 
  }


  
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  cout << "\n\nNow check if the matrix Q is in the form-III using Reinschke algorithm.\n";
  // (1) initialization: I=n; J=n+m; V=0;
  int I = n;
  int J = n+m;
  int V = 0;

  while(I>0) {
    cout << "I = " << I << " J = " << J; // 

    //cout << " Q="; MatrixPrint(Q, n, n+m);

    // (2) set V equal to the number of non-zeros of that column I_s of the IxJ submatrix von Q=[A,B], which fulfills
    //     the following conditions:
    //     * The column I_s is different from zero; 
    //     * The column I_s contains a minimum number of non-zeros.
    vector<int> Nnzero(J+1,0); // Nnzero[j] = # of non-zeros in column j
    for(int j=1; j<=J; j++) {
      for(int i=1; i<=I; i++) {
	if(Q[i][j]!=0)
	  Nnzero[j]++;
      }
    }
    
    int I_s = 1;
    int Nnz_min = I;
    for(int j=1; j<=J; j++) {
      if(Nnzero[j]>0 && Nnzero[j]<=Nnz_min) {
	Nnz_min = Nnzero[j];
	I_s = j;
      }
    }
    V = Nnz_min;
    //cout << "I_s = " << I_s << "; J = " << J << "; V = " << V << endl; //test
    cout << "; V = " << V << endl; //test
 
    // (3) if V!=1 then the system is not strongly structurally controllable.
    if(V!=1) {
      cout << "\nThe system Q=[A,B] is not strongly structurally controllable.\n";
      break;
    }

    // (4) set I_z equal to the row index of the last non-zero in the column I_s of the IxJ submatrix von Q=[A,B]
    int I_z = 0;
    for(int i=1; i<=I; i++) {
      if(Q[i][I_s]!=0)
	I_z = i;
    }
    //cout << "I_z = " << I_z << "; I = " << I << endl; //test

    
    // (5) if I_z != I then permute the row I_z and I of the matrix Q=[A,B]
    if(I_z != I) 
      Permute_Row(Q, rowindex, n,n+m,I_z,I);

    
    // (6) if I_s != J then permute the column I_s and J of the matrix Q=[A,B]
    if(I_s != J)
      Permute_Col(Q, colindex, n,n+m,I_s,J);
    
    
    // (7) I=I-1; J=J-1; V=0;
    I--;
    J--;
    V=0;

  }

  cout << " Q="; MatrixPrint(Q, n, n+m);

  cout << "row-col index:\n";
  for(int i=1; i<=n; i++) { 
    for(int j=1; j<=n+m; j++) {
      cout << rowindex[i][j] << colindex[i][j] << " ";
    }
    cout << endl; 
  }


  // (8) if I=0; then term rank of Q=[A,B]=n
  if(I==0) {
    cout << "term rank Q = " << n <<"; Q is in the form-III\n";

    // Now consider the second condition of Reinschke's theorem
    int n_aii = 0;
    for(int i=1; i<=n; i++) {
      if(rowindex[i][m+i]==i && colindex[i][m+i]==i) {
	n_aii++;
      }
    }

    if(n_aii>0) 
      cout << "\nThe x-elements of the form-III  of [A-lambda I, B] contain a_ii element \n"
	   << "==> [A,B] is NOT strongly structurally controllable.\n\n";
    else
      cout << "\nThe x-elements of the form-III  of [A-lambda I, B] contain no a_ii element \n"
	   << "==> [A,B] is  strongly structurally controllable.\n\n";
  }
   

  free_imatrix(Q0,1,n,1,n+m);
  free_imatrix(Q,1,n,1,n+m);
  free_imatrix(rowindex,1,n,1,n+m);
  free_imatrix(colindex,1,n,1,n+m);

}
//////////////////////////////////////////////////////////////////////






///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/* Check whether a Matrix Q is in the Form-III, i.e. there exists two permutation matries P1, P2 such that 

             [ @ ... @ x 0 ..........0 ]
             | .     . @ x 0 .       . |
             | .     . . @ x 0 .     . |
   P1 Q P2 = | .     . .   @ x 0 .   . |
             | .     . .     @ x 0 . . |
             | .     . .       @ x 0 . |
             | .     . .         @ x 0 |
             [ @ ... @ @ ..........@ x ]

   where the x-elements must be non-zero. The @-elements can be either zero or non-zero. 

   This is used to determine whether a structured system [A,B] is strongly structural controllable.
   For details, see the paper entitled "On strong structural controllability of linear systems" by K.J. Reinschke et al.

   Note that here, Q is given by the structured system [A,B].
*/
 
bool Transform_to_MatrixFormIII(int**& Q, int n, int m, int**& rowindex, int**& colindex)
{
  Rand rand;
  rand.seed(1);

  //store the original row index, col index of each element in Q
  for(int i=1; i<=n; i++) { 
    for(int j=1; j<=n+m; j++) {
      rowindex[i][j] = i;
      colindex[i][j] = j;
    }
  } 
  
  //store the original Q matrix
  int** Q0 = imatrix(1,n,1,n+m);
  MatrixCopy(Q, Q0, n, n+m);

  
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  cout << "\n\nNow check if the matrix Q is in the form-III using Reinschke algorithm.\n";
  // (1) initialization: I=n; J=n+m; V=0;
  int I = n;
  int J = n+m;
  int V = 0;

  while(I>0) {
    cout << "I = " << I << " J = " << J; // 

    // (2) set V equal to the number of non-zeros of that column I_s of the IxJ submatrix von Q=[A,B], which fulfills
    //     the following conditions:
    //     * The column I_s is different from zero; 
    //     * The column I_s contains a minimum number of non-zeros.
    vector<int> Nnzero(J+1,0); // Nnzero[j] = # of non-zeros in column j
    for(int j=1; j<=J; j++) {
      for(int i=1; i<=I; i++) {
	if(Q[i][j]!=0)
	  Nnzero[j]++;
      }
    }
    
    int I_s = 1;
    int Nnz_min = I;
    for(int j=1; j<=J; j++) {
      if(Nnzero[j]>0 && Nnzero[j]<=Nnz_min) {
	Nnz_min = Nnzero[j];
	I_s = j;
      }
    }
    V = Nnz_min;
    //cout << "I_s = " << I_s << "; J = " << J << "; V = " << V << endl; //test
    cout << "; V = " << V << endl; //test
 
    // (3) if V!=1 then the system is not strongly structurally controllable.
    if(V!=1) {
      cout << "\nThe system Q=[A,B] is not strongly structurally controllable.\n";
      break;
    }

    // (4) set I_z equal to the row index of the last non-zero in the column I_s of the IxJ submatrix von Q=[A,B]
    int I_z = 0;
    for(int i=1; i<=I; i++) {
      if(Q[i][I_s]!=0)
	I_z = i;
    }
    //cout << "I_z = " << I_z << "; I = " << I << endl; //test

    
    // (5) if I_z != I then permute the row I_z and I of the matrix Q=[A,B]
    if(I_z != I) 
      Permute_Row(Q, rowindex, n,n+m,I_z,I);

    
    // (6) if I_s != J then permute the column I_s and J of the matrix Q=[A,B]
    if(I_s != J)
      Permute_Col(Q, colindex, n,n+m,I_s,J);
    
    
    // (7) I=I-1; J=J-1; V=0;
    I--;
    J--;
    V=0;

  }

  //cout << " Q0="; MatrixPrint(Q0, n, n+m);
  //cout << " Q="; MatrixPrint(Q, n, n+m);

  // (8) if I=0; then term rank of Q=[A,B]=n
  bool QIII = false;
  if(I==0) {
    cout << "term rank = " << n <<"; in the form-III\n";
    QIII = true;
  }
  else
    cout << "term rank = " << n <<"; NOT in the form-III\n";


  return QIII;
}
//////////////////////////////////////////////////////////////////////








////////////////////////////////////////////////////////////////////////
// Bisection search!
// Give an array xx[0..n-1], and given a value x, 
// return a value j such that x is between xx[j-1] and xx[j].
// xx must be monotonic, either increasing or decreasing. 
// j= -1 or n is returned to indicate that x is out of range.  
int BisectionSearch(vector<double>& xx, const double x)
{
  int ju,jm,jl;
  bool ascnd;

  int n=xx.size();
  jl=-1; // initialize lower 
  ju=n;  // and upper limits
  ascnd=(xx[n-1] >= xx[0]); // true if ascending order of table, false otherwise

  while (ju-jl > 1) {
    jm=(ju+jl) >> 1; // compute a midpoint
    if ((x >= xx[jm]) == ascnd) jl=jm;       // replace either the lower limit
    else        	          ju=jm;       // or the upper limit, as appropriate       
  }

  int j;
  if (x == xx[0]) j=0;
  else if (x == xx[n-1]) j=n-1;
  else j=jl+1;

  return j;
}






// Bisection search!
// Give an array xx[1..n], and given a value x, 
// return a value j such that x is between xx[j-1] and xx[j].
// xx must be monotonic, either increasing or decreasing. 
// j= 0 or n+1 is returned to indicate that x is out of range.  
int BisectionSearch(double* xx, int n, const double x)
{
  int ju,jm,jl;
  bool ascnd;

  //int n=xx.size();
  jl=0; // initialize lower 
  ju=n+1;  // and upper limits
  ascnd=(xx[n] >= xx[1]); // true if ascending order of table, false otherwise

  while (ju-jl > 1) {
    jm=(ju+jl) >> 1; // compute a midpoint
    if ((x >= xx[jm]) == ascnd) jl=jm;       // replace either the lower limit
    else        	          ju=jm;       // or the upper limit, as appropriate       
  }

  int j;
  if (x == xx[1]) j=1;
  else if (x == xx[n]) j=n;
  else j=jl+1;

  return j;
}

////////////////////////////////////////////////////////////////////////



/////////////////////////////////////////////////////////////////////////////////////////////////
// find whether the element j is in the list nbl
bool find(Nbl& nbl, int j)
{
  for(Nbl::iterator p = nbl.begin(); p!= nbl.end(); p++) {
    if(j==(*p)) 
      return true;
  }
    
  return false;
}
/////////////////////////////////////////////////////////////////////////////////////////////////



/////////////////////////////////////////////////////////////////////////////////////////////////
bool find(set<int>& nbl, int j)
{
  for(set<int>::iterator p = nbl.begin(); p!= nbl.end(); p++) {
    if(j==(*p)) 
      return true;
  }
    
  return false;
}
/////////////////////////////////////////////////////////////////////////////////////////////////



/////////////////////////////////////////////////////////////////////////////////////////////////
bool find(vector<set<int> >& Vset, set<int>& S)
{
  for(vector<set<int> >::iterator p = Vset.begin(); p!= Vset.end(); p++) {
    if(S==(*p)) 
      return true;
  }
    
  return false;
}
/////////////////////////////////////////////////////////////////////////////////////////////////


/////////////////////////////////////////////////////////////////////////////////////////////////
// find whether the element j is in the list nbl
bool find(vector<vector<int> >& Vv, vector<int>& v)
{
  for(vector<vector<int> >::iterator p = Vv.begin(); p!= Vv.end(); p++) {
    if(v==(*p)) 
      return true;
  }
    
  return false;
}
/////////////////////////////////////////////////////////////////////////////////////////////////



/////////////////////////////////////////////////////////////////////////////////////////////////
// find whether the element j is in the list nbl
bool find(vector<vector<string> >& Vs, vector<string>& v)
{
  for(vector<vector<string> >::iterator p = Vs.begin(); p!= Vs.end(); p++) {
    if(v==(*p)) 
      return true;
  }
    
  return false;
}
/////////////////////////////////////////////////////////////////////////////////////////////////




/////////////////////////////////////////////////////////////////////////////////////////////////
// find whether the state s is in the state pool
bool find(vector< vector<char>  >& Vstate, vector<char>& state)
{
  for(vector< vector<char> >::iterator p = Vstate.begin(); p!= Vstate.end(); p++) {
    if(state==(*p)) 
      return true;
  }
    
  return false;
}
/////////////////////////////////////////////////////////////////////////////////////////////////




/////////////////////////////////////////////////////////////////////////////////////////////////
// find whether the element j is in the vector nbv
bool find(vector<int>& nbv, int j)
{
  int n = nbv.size();
  for(int i=0; i<n; i++)
    if(j==nbv[i]) 
      return true;

  return false;
}
/////////////////////////////////////////////////////////////////////////////////////////////////



/////////////////////////////////////////////////////////////////////////////////////////////////
// find whether the element j is in the vector nbv
bool find(vector<string>& nbs, string s)
{
  int n = nbs.size();
  for(int i=0; i<n; i++)
    if(s==nbs[i]) 
      return true;

  return false;
}
/////////////////////////////////////////////////////////////////////////////////////////////////



/////////////////////////////////////////////////////////////////////////////////////////////////
// find whether the element x is in the deque Q
bool find(deque<int>& Q, int x)
{
  for(deque<int>::iterator q = Q.begin(); q!= Q.end(); q++) {
    if(x==(*q)) 
      return true;
  }
  
  return false;
}
/////////////////////////////////////////////////////////////////////////////////////////////////


/////////////////////////////////////////////////////////////////////////////////////////////////
// find whether the element x is in the vector<Link> L
bool find(list<ULink>& L, ULink e)
{
  for(list<ULink>::iterator q = L.begin(); q!= L.end(); q++) {
    if(e==(*q)) 
      return true;
  }
  
  return false;
}
/////////////////////////////////////////////////////////////////////////////////////////////////



/////////////////////////////////////////////////////////////////////////////////////////////////
// find whether the element x is in the vector<Link> L
bool find(list<Link>& L, Link e)
{
  for(list<Link>::iterator q = L.begin(); q!= L.end(); q++) {
    if(e==(*q)) 
      return true;
  }
  
  return false;
}
/////////////////////////////////////////////////////////////////////////////////////////////////


bool find(list<int>& X, int j) {
  for(list<int>::iterator p = X.begin(); p!= X.end(); p++) {
    if(j==(*p)) 
      return true;
  }
  return false;
}

bool find(vector<string>& X, string j) {
  for(vector<string>::iterator p = X.begin(); p!= X.end(); p++) {
    if(j==(*p)) 
      return true;
  }
  return false;
}

bool find(vector<set<int> >& X, set<int>& j) {
  for(vector<set<int> >::iterator p = X.begin(); p!= X.end(); p++) {
    if(j==(*p)) 
      return true;
  }
  return false;
}

bool find(vector<list<int> >& X, list<int>& j) {
  for(vector<list<int> >::iterator p = X.begin(); p!= X.end(); p++) {
    if(j==(*p)) 
      return true;
  }
  return false;
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void print(deque<int>& Q) 
{
  for(deque<int>::iterator q = Q.begin(); q!= Q.end(); q++) 
    cout << (*q) << "-->";
  cout << endl;
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


/////////////////////////////////////////////////////////////////////////////////////////////////
// find the position of the element j in the list nbl  (assume j is definitely in the list)
int findpos(Nbl& nbl, int j)
{
  int pos=0;
  for(Nbl::iterator p = nbl.begin(); p!= nbl.end(); p++) {
    if(j==(*p)) 
      return pos;
    else
      pos++;
  }
}
/////////////////////////////////////////////////////////////////////////////////////////////////


/////////////////////////////////////////////////////////////////////////////////////////////////
// find the position of the element j in the vector nbv  (assume j is definitely in the vector)
int findpos(vector<int>& nbv, int j)
{
  int n = nbv.size();
  for(int i=0; i<n; i++)
    if(j==nbv[i]) 
      return i;
  
  return -1;
}
/////////////////////////////////////////////////////////////////////////////////////////////////



/////////////////////////////////////////////////////////////////////////////////////////////////
// find the position of the element j in the vector nbv  (assume j is definitely in the vector)
int findpos(vector<string>& nbs, string s)
{
  int n = nbs.size();
  for(int i=0; i<n; i++)
    if(s==nbs[i]) 
      return i;
  
  return -1;
}
/////////////////////////////////////////////////////////////////////////////////////////////////


/////////////////////////////////////////////////////////////////////////////////////////////////
// return the element in the list nbl with a given position (index)
int print(Nbl& nbl)
{
  Nbl::iterator p = nbl.begin();
  cout << *p;
  p++;
  for(; p!= nbl.end(); p++) {
    cout << "-->" << *p;
  }
  cout << endl;
}
/////////////////////////////////////////////////////////////////////////////////////////////////



/////////////////////////////////////////////////////////////////////////////////////////////////
// return the element in the list nbl with a given position (index)
int getnodeinNbl(Nbl& nbl, int index)
{
  int pos=0;
  for(Nbl::iterator p = nbl.begin(); p!= nbl.end(); p++, pos++) {
    if(pos==index) 
      return (*p);
  }
}
/////////////////////////////////////////////////////////////////////////////////////////////////



/////////////////////////////////////////////////////////////////////////////////////////////////
// return the element in the list X with a given position (index)
double getnodeinlist(list<double>& X, int index)
{
  int pos=0;
  for(list<double>::iterator p = X.begin(); p!= X.end(); p++, pos++) {
    if(pos==index) 
      return (*p);
  }
}

int getnodeinlist(list<int>& X, int index)
{
  int pos=0;
  for(list<int>::iterator p = X.begin(); p!= X.end(); p++, pos++) {
    if(pos==index) 
      return (*p);
  }
}
/////////////////////////////////////////////////////////////////////////////////////////////////


///////////////////////////////////////////////////////////////////////////////
template <typename Iter> 
void getNbl(Iter begin, Iter end) //, Nbl& nbl) 
{ 
  int count=1;
  int i = *begin++;
  cout << i << ' ';

  for( ; begin != end ; count++) {
    int j = *begin++;
    cout << j << ' ';
  }
  cout << endl;
  //cout << "count= " << count << endl;
} 
///////////////////////////////////////////////////////////////////////////////




/////////////////////////////////////////////////////////////////////////////////////////////////
const string numbers="0123456789";
// get the number from a string
void GetNumberFromaString(string str, int& x) 
{
  // x = atoi((str.substr(11)).c_str()); // this just works if we know exactly the posisition of the first number
  // in general, we have to do the following:
  int first = str.find_first_of(numbers);
  if(first == string::npos) 
    cout<<"find no numbers"<<endl;

  int last = str.find_last_of(numbers);
  if(last == string::npos) 
    cout<<"find no numbers"<<endl;

  // the following statements have the same function: convert string to int
  // x = atoi((str.substr(first, last-first+1)).c_str());
  istringstream buffer(str.substr(first, last-first+1)); 
  buffer >> x;
  
}
/////////////////////////////////////////////////////////////////////////////////////////////////

#include <sstream>
void Char2String(char& c, string& s) {
  stringstream ss;
  ss << c;
  ss >> s;
}


string Int2String(int i) {
  string s;
  stringstream ss;
  ss << i;
  ss >> s;
  return s;
}



/////////////////////////////////////////////////////////////////////////////////////////////////
/**
 * cpx_polylog_sum -- compute the polylogarithm by direct summation
 *
 * Li_s(z) = sum_{n=1}^infty z^n/ n^s
 * 
 * The magnitude of z must be less than one in order for the 
 * summation to be carried out.
 */
double polylog(double s, double z)
{
  double sum = 0;
  double zn = z;
  int nmax = 1;
  for(int n=1; ; n++) {
    double term = zn/pow(n,s);
    sum += term;
    zn *= z;
    if(term < TOLERANCE) {
      nmax = n;
      break;
    }
  }
  //cout << "n= " << nmax << endl; //test

  return sum;
}
/////////////////////////////////////////////////////////////////////////////////////////////////


/*
/////////////////////////////////////////////////////////////////////////////////////////////////
double polylog(double s, double z, int prec)
{
cpx_t plog;
cpx_init (plog);

cpx_t ess;
cpx_init (ess);
cpx_set_d(ess, s, 0);

cpx_t zee;
cpx_init (zee);
cpx_set_d(zee, z, 0);

cpx_polylog_euler (plog, ess, zee, prec);
//cout << "cpx_polylog_euler: Li_" << s << "[" << z << "] = " << plog[0].re << endl;

return mpf_get_d(plog[0].re);
}
/////////////////////////////////////////////////////////////////////////////////////////////////
*/


 /////////////////////////////////////////////////////////////////////////////////////////////////
 // find the root(s) of f(x)=0 in the interval [X1, X2]
int FindRoots(double fx(const double), double X1, double X2, double& rootmin, double& rootmax, vector<double>& roots)
{
  const int n=100,NBMAX=20;
  //const int n=10000,NBMAX=20;  

  int i,nb=NBMAX;
  DP xacc,root;
  Vec_DP xb1(NBMAX),xb2(NBMAX);
  
  zbrak(fx,X1,X2,n,xb1,xb2,nb);
  cout << endl << "Roots of f:" << endl;
  cout << setw(20) << "x" << setw(16) << "f(x)" << endl;
  //cout << fixed << setprecision(10);

  roots.clear();
  roots.resize(nb,0);

  rootmin = 1e300;
  rootmax = -1e300;
  for (i=0;i<nb;i++) {
    //xacc = (1.0e-10)*(xb1[i]+xb2[i])/2.0;
    xacc = (1.0e-10)*(fabs(xb1[i])+fabs(xb2[i]))/2.0;
    root = rtbis(fx,xb1[i],xb2[i],xacc);
    cout << "root " << setw(3) << (i+1) << setw(15) << root;
    cout << setw(15) << fx(root) << endl;

    if (root < rootmin) {
      rootmin = root;
    }

    if (root > rootmax) {
      rootmax = root;
    }

    roots[i] = root;
  }

  cout << "nb= " << nb << endl; 
  
  return nb;
  //cout << "root_min = " << rootmin << endl;
  //cout << "root_max = " << rootmax << endl;

}
///////////////////////////////////////////////////////////////////////////////////////////////////





///////////////////////////////////////////////////////////////////////////////////////////////////
void zbrak_2D(double fx(const double, const double), 
	      const double a, 
	      const double x1, const double x2, 
	      const int n,
	      vector<double> &xb1, vector<double> &xb2, 
	      int &nroot)
{
  int nb=xb1.size();
  nroot=0;
  
  double dx=(x2-x1)/n;
  double x=x1;
  double fp=fx(x1,a);

  for (int i=0;i<n;i++) {
    x+=dx;
    double fc=fx(x,a);
    //cout << x-dx << ',' << x << ',' << fp << ',' << fc << endl; //debug
    if (fc*fp <= 0.0) {
      xb1[nroot]=x-dx;
      xb2[nroot++]=x;
      if(nroot == nb) return;
    }
    fp=fc;
  }

  //cout << nroot << endl; //debug
}
///////////////////////////////////////////////////////////////////////////////////////////////////



 /////////////////////////////////////////////////////////////////////////////////////////////////
double rtbis_2D(double func(const double, const double), 
		const double a,
		const double x1, const double x2, const double xacc)
{
  const int JMAX=40;
  int j;
  double dx,f,fmid,xmid,rtb;

  f=func(x1,a);
  fmid=func(x2,a);
  //if (f*fmid >= 0.0) nrerror("Root must be bracketed for bisection in rtbis");                                   
  if (f*fmid > 1e-16) nrerror("Root must be bracketed for bisection in rtbis");
  rtb = f < 0.0 ? (dx=x2-x1,x1) : (dx=x1-x2,x2);
  for (j=0;j<JMAX;j++) {
    fmid=func(xmid=rtb+(dx *= 0.5),a);
    if (fmid <= 0.0) rtb=xmid;
    if (fabs(dx) < xacc || fmid == 0.0) return rtb;
  }
  nrerror("Too many bisections in rtbis");
  return 0.0;
}
 /////////////////////////////////////////////////////////////////////////////////////////////////


 /////////////////////////////////////////////////////////////////////////////////////////////////
 // find the root(s) of f(x;a)=0 in the interval [X1, X2]
int FindRoots(double fx(const double, const double), 
	      const double a,
	      double X1, double X2, double& rootmin, double& rootmax, vector<double>& roots)
{
  const int n=100,NBMAX=20;
  int i,nb=NBMAX;
  double xacc,root;
  vector<double> xb1(NBMAX),xb2(NBMAX);
  
  zbrak_2D(fx,a,X1,X2,n,xb1,xb2,nb);
  //cout << endl << "Roots of F:" << endl;
  //cout << setw(20) << "x" << setw(16) << "F(x)" << endl;
  //cout << fixed << setprecision(10);

  roots.clear();
  roots.resize(nb,0);

  rootmin = 1;//1e300;
  rootmax = 0;//-1e300;
  for (i=0;i<nb;i++) {
    //xacc = (1.0e-10)*(xb1[i]+xb2[i])/2.0;
    xacc = (1.0e-11)*(fabs(xb1[i])+fabs(xb2[i]))/2.0;
    root = rtbis_2D(fx,a,xb1[i],xb2[i],xacc);
    //cout << "Root " << setw(3) << (i+1) << setw(15) << root;
    //cout << setw(15) << fx(root,a) << endl;

    if (root < rootmin) {
      rootmin = root;
    }

    if (root > rootmax) {
      rootmax = root;
    }

    roots[i] = root;
  }

  //cout << "nb= " << nb << endl; 
  
  return nb;
  //cout << "root_min = " << rootmin << endl;
  //cout << "root_max = " << rootmax << endl;

}
///////////////////////////////////////////////////////////////////////////////////////////////////




///////////////////////////////////////////////////////////////////////////////////////////////////
void zbrak_3D(double fx(const double, const double, const double), 
	      const double a, const double b, 
	      const double x1, const double x2, 
	      const int n,
	      vector<double> &xb1, vector<double> &xb2, 
	      int &nroot)
{
  int nb=xb1.size();
  nroot=0;

 double dx=(x2-x1)/n;
  double x=x1;
  double fp=fx(x1,a,b);

  for (int i=0;i<n;i++) {
    x+=dx;
    double fc=fx(x,a,b);
    //cout << x-dx << ',' << x << ',' << fp << ',' << fc << endl; //debug
    if (fc*fp <= 0.0) {
      xb1[nroot]=x-dx;
      xb2[nroot++]=x;
      if(nroot == nb) return;
    }
    fp=fc;
  }

}
///////////////////////////////////////////////////////////////////////////////////////////////////



 /////////////////////////////////////////////////////////////////////////////////////////////////
double rtbis_3D(double func(const double, const double, const double), 
		const double a, const double b,
		const double x1, const double x2, const double xacc)
{
  const int JMAX=40;
  int j;
  double dx,f,fmid,xmid,rtb;

  f=func(x1,a,b);
  fmid=func(x2,a,b);
  //if (f*fmid >= 0.0) nrerror("Root must be bracketed for bisection in rtbis");                                   
  if (f*fmid > 1e-16) mynrerror("Root must be bracketed for bisection in rtbis");
  rtb = f < 0.0 ? (dx=x2-x1,x1) : (dx=x1-x2,x2);
  for (j=0;j<JMAX;j++) {
    fmid=func(xmid=rtb+(dx *=0.5), a, b);
    if (fmid <= 0.0) rtb=xmid;
    if (fabs(dx) < xacc || fmid == 0.0) return rtb;
  }
  nrerror("Too many bisections in rtbis");
  return 0.0;
}
 /////////////////////////////////////////////////////////////////////////////////////////////////


 /////////////////////////////////////////////////////////////////////////////////////////////////
 // find the root(s) of f(x;a,b)=0 in the interval [X1, X2]
int FindRoots(double fx(const double, const double, const double), 
	      const double a, const double b,
	      double X1, double X2, double& rootmin, double& rootmax, vector<double>& roots)
{
  const int n=100,NBMAX=20;
  int i,nb=NBMAX;
  double xacc,root;
  vector<double> xb1(NBMAX),xb2(NBMAX);
  
  zbrak_3D(fx,a,b,X1,X2,n,xb1,xb2,nb);
  //cout << endl << "Roots of f:" << endl;
  //cout << setw(20) << "x" << setw(16) << "f(x)" << endl;
  //cout << fixed << setprecision(10);

  roots.clear();
  roots.resize(nb,0);

  rootmin = 1;//1e300;
  rootmax = 0;//-1e300;
  for (i=0;i<nb;i++) {
    //xacc = (1.0e-10)*(xb1[i]+xb2[i])/2.0;
    xacc = (1.0e-11)*(fabs(xb1[i])+fabs(xb2[i]))/2.0;
    root = rtbis_3D(fx,a,b,xb1[i],xb2[i],xacc);
    //cout << "root " << setw(3) << (i+1) << setw(15) << root;
    //cout << setw(15) << fx(root,a,b) << endl;

    if (root < rootmin) {
      rootmin = root;
    }

    if (root > rootmax) {
      rootmax = root;
    }

    roots[i] = root;
  }

  //cout << "nb= " << nb << endl; 
  
  return nb;
  //cout << "root_min = " << rootmin << endl;
  //cout << "root_max = " << rootmax << endl;

}
///////////////////////////////////////////////////////////////////////////////////////////////////




///////////////////////////////////////////////////////////////////////////////////////////////////
void zbrak_long(long double fx(const long double), const long double x1, const long double x2, const int n,
		vector<long double> &xb1, vector<long double> &xb2, int &nroot)
{
  int i;
  long double x,fp,fc,dx;

  int nb=xb1.size();
  nroot=0;
  dx=(x2-x1)/n;
  fp=fx(x=x1);
  for (i=0;i<n;i++) {
    fc=fx(x += dx);
    if (fc*fp <= 0.0) {
      xb1[nroot]=x-dx;
      xb2[nroot++]=x;
      if(nroot == nb) return;
    }
    fp=fc;
  }
}
///////////////////////////////////////////////////////////////////////////////////////////////////



///////////////////////////////////////////////////////////////////////////////////////////////////
long double rtbis_long(long double func(const long double), const long double x1, const long double x2, const long double xacc)
{
  const int JMAX=40;
  int j;
  long double dx,f,fmid,xmid,rtb;

  f=func(x1);
  fmid=func(x2);
  //if (f*fmid >= 0.0) nrerror("Root must be bracketed for bisection in rtbis");
  if (f*fmid > 1e-16) mynrerror("Root must be bracketed for bisection in rtbis");
  rtb = f < 0.0 ? (dx=x2-x1,x1) : (dx=x1-x2,x2);
  for (j=0;j<JMAX;j++) {
    fmid=func(xmid=rtb+(dx *= 0.5));
    if (fmid <= 0.0) rtb=xmid;
    if (fabs(dx) < xacc || fmid == 0.0) return rtb;
  }
  nrerror("Too many bisections in rtbis");
  return 0.0;
}
///////////////////////////////////////////////////////////////////////////////////////////////////




/////////////////////////////////////////////////////////////////////////////////////////////////
// find the root(s) of f(x)=0 in the interval [X1, X2]
int FindRoots_long(long double fx(const long double), long double X1, long double X2, long double& rootmin, long double& rootmax, vector<long double>& roots)
{
  const int n=100,NBMAX=20;
  int i,nb=NBMAX;
  long double xacc,root;
  vector <long double> xb1(NBMAX),xb2(NBMAX);
  
  zbrak_long(fx,X1,X2,n,xb1,xb2,nb);
  cout << endl << "Roots of f:" << endl;
  cout << setw(20) << "x" << setw(16) << "f(x)" << endl;
  //cout << fixed << setprecision(10);

  roots.clear();
  roots.resize(nb,0);

  rootmin = 1e300;
  rootmax = -1e300;
  for (i=0;i<nb;i++) {
    //xacc = (1.0e-10)*(xb1[i]+xb2[i])/2.0;
    xacc = (1.0e-11)*(fabs(xb1[i])+fabs(xb2[i]))/2.0;
    root = rtbis_long(fx,xb1[i],xb2[i],xacc);
    cout << "root " << setw(3) << (i+1) << setw(15) << root;
    cout << setw(15) << fx(root) << endl;

    if (root < rootmin) {
      rootmin = root;
    }

    if (root > rootmax) {
      rootmax = root;
    }

    roots[i] = root;
  }

  cout << "nb= " << nb << endl; 
  
  return nb;
  //cout << "root_min = " << rootmin << endl;
  //cout << "root_max = " << rootmax << endl;

}
///////////////////////////////////////////////////////////////////////////////////////////////////




/////////////////////////////////////////////////////////////////////////////////////////////////
// find the root(s) of f(x)=0 in the interval [X1, X2]
// we first find roots in a slightly bigger interval [X1', X2'] with X1' = X1- (X2-X1)*0.01
//                                                                   X2' = X2+ (X2-X1)*0.01
// then only choose roots within the interval [X1, X2]
int FindRoots_u(double fx(const double), double X1, double X2, vector<double>& roots)
{
  const int n=100,NBMAX=20;
  int i,nb=NBMAX;
  DP xacc,root;
  Vec_DP xb1(NBMAX),xb2(NBMAX);
  
  double X1u = X1 - (X2-X1)*0.01;
  double X2u = X2 + (X2-X1)*0.01;

  double X1u2 = X1 - 1e-10;
  double X2u2 = X2 + 1e+10;

  zbrak(fx,X1u,X2u,n,xb1,xb2,nb);
  cout << endl << "Roots of f:" << endl;
  cout << setw(20) << "x" << setw(16) << "f(x)" << endl;
  //cout << fixed << setprecision(10);

  roots.clear();
  for (i=0;i<nb;i++) {
    //xacc = (1.0e-10)*(xb1[i]+xb2[i])/2.0;
    xacc = (1.0e-10)*(fabs(xb1[i])+fabs(xb2[i]))/2.0;
    root = rtbis(fx,xb1[i],xb2[i],xacc);

    if(fabs(root-X1)<1e-10)
      root = X1;
    if(fabs(root-X2)<1e-10)
      root = X2;
    
    cout << "root " << setw(3) << (i+1) << setw(15) << root;
    cout << setw(15) << fx(root) << endl;

    if(nb>1) {
      if(root<=X2u2 && root>=X1u2)
	roots.push_back(root);
    }
    else 
      roots.push_back(root);
  }

  nb = roots.size();
  cout << "nb= " << nb << endl; 
  return nb;
}
///////////////////////////////////////////////////////////////////////////////////////////////////



void zbrak_rrg(double frx(const int r, const double x), const int r, const double x1, const double x2, const int n,
	       Vec_O_DP &xb1, Vec_O_DP &xb2, int &nroot)
{
  int i;
  DP x,fp,fc,dx;

  int nb=xb1.size();
  nroot=0;
  dx=(x2-x1)/n;
  fp=frx(r, x=x1);
  for (i=0;i<n;i++) {
    fc=frx(r, x += dx);
    if (fc*fp <= 0.0) {
      xb1[nroot]=x-dx;
      xb2[nroot++]=x;
      if(nroot == nb) return;
    }
    fp=fc;
  }
}


double rtbis_rrg(DP func(const int, const DP), const int r, const DP x1, const DP x2, const DP xacc)
{
  const int JMAX=40;
  int j;
  DP dx,f,fmid,xmid,rtb;

  f=func(r, x1);
  fmid=func(r, x2);
  //if (f*fmid >= 0.0) nrerror("Root must be bracketed for bisection in rtbis");
  if (f*fmid > 1e-16) nrerror("Root must be bracketed for bisection in rtbis");
  rtb = f < 0.0 ? (dx=x2-x1,x1) : (dx=x1-x2,x2);
  for (j=0;j<JMAX;j++) {
    fmid=func(r, xmid=rtb+(dx *= 0.5));
    if (fmid <= 0.0) rtb=xmid;
    if (fabs(dx) < xacc || fmid == 0.0) return rtb;
  }
  nrerror("Too many bisections in rtbis");
  return 0.0;
}



/////////////////////////////////////////////////////////////////////////////////////////////////
// designed in particular for random regular graph's population dynamics 
// find the root(s) of f(x)=0 in the interval [X1, X2]
// we first find roots in a slightly bigger interval [X1', X2'] with X1' = X1- (X2-X1)*0.01
//                                                                   X2' = X2+ (X2-X1)*0.01
// then only choose roots within the interval [X1, X2]
int FindRoots_RRG(double frx(const int, const double), const int r, double X1, double X2, vector<double>& roots)
{
  const int n=100,NBMAX=20;
  int i,nb=NBMAX;
  DP xacc,root;
  Vec_DP xb1(NBMAX),xb2(NBMAX);
  
  double X1u = X1 - (X2-X1)*0.01;
  double X2u = X2 + (X2-X1)*0.01;

  double X1u2 = X1 - 1e-10;
  double X2u2 = X2 + 1e+10;

  zbrak_rrg(frx,r,X1u,X2u,n,xb1,xb2,nb);
  cout << endl << "Roots of f:" << endl;
  cout << setw(20) << "x" << setw(16) << "f(x)" << endl;
  //cout << fixed << setprecision(10);

  roots.clear();
  for (i=0;i<nb;i++) {
    //xacc = (1.0e-10)*(xb1[i]+xb2[i])/2.0;
    xacc = (1.0e-10)*(fabs(xb1[i])+fabs(xb2[i]))/2.0;
    root = rtbis_rrg(frx,r,xb1[i],xb2[i],xacc);

    if(fabs(root-X1)<1e-10)
      root = X1;
    if(fabs(root-X2)<1e-10)
      root = X2;
    
    cout << "root " << setw(3) << (i+1) << setw(15) << root;
    cout << setw(15) << frx(r, root) << endl;

    if(nb>1) {
      if(root<=X2u2 && root>=X1u2)
	roots.push_back(root);
    }
    else 
      roots.push_back(root);
  }

  nb = roots.size();
  cout << "nb= " << nb << endl; 

  if(nb==0) {
    cout << "No solutions!\n ";
    exit(0);
  }

  return nb;
}
///////////////////////////////////////////////////////////////////////////////////////////////////





//////////////////////////////////////////////////////////////////////////////////////////////////
// return the binomial coefficient: C(m,k)
double binom(int m, int k)
{
  double c = 1.;
  for(int i=0; i<k; i++) 
    c *= (m-i)/(i+1.);
  return c;
}
//////////////////////////////////////////////////////////////////////////////////////////////////





//////////////////////////////////////////////////////////////////////////////////////////////////
void NR_gser(double &gamser, const double a, const double x)
{
  const int ITMAX=100;
  const DP EPS=numeric_limits<double>::epsilon();
  int n;
  double sum,del,ap;

  //gln=gammln(a);
  if (x <= 0.0) {
    if (x < 0.0) mynrerror("x less than 0 in routine gser");
    gamser=0.0;
    return;
  } else {
    ap=a;
    del=sum=1.0/a;
    for (n=0;n<ITMAX;n++) {
      ++ap;
      del *= x/ap;
      sum += del;
      if (fabs(del) < fabs(sum)*EPS) {
	//gamser=sum*exp(-x+a*log(x)-gln);
	gamser=sum*exp(-x+a*log(x));
	return;
      }
    }
    mynrerror("a too large, ITMAX too small in routine gser");
    return;
  }
}




void NR_gcf(double &gammcf, const double a, const double x)
{
  const int ITMAX=100;
  const double EPS=numeric_limits<double>::epsilon();
  const double FPMIN=numeric_limits<double>::min()/EPS;
  int i;
  double an,b,c,d,del,h;

  //gln=gammln(a);
  b=x+1.0-a;
  c=1.0/FPMIN;
  d=1.0/b;
  h=d;
  for (i=1;i<=ITMAX;i++) {
    an = -i*(i-a);
    b += 2.0;
    d=an*d+b;
    if (fabs(d) < FPMIN) d=FPMIN;
    c=b+an/c;
    if (fabs(c) < FPMIN) c=FPMIN;
    d=1.0/d;
    del=d*c;
    h *= del;
    if (fabs(del-1.0) <= EPS) break;
  }
  if (i > ITMAX) nrerror("a too large, ITMAX too small in gcf");
  //gammcf=exp(-x+a*log(x)-gln)*h;
  gammcf=exp(-x+a*log(x))*h;
}



double NR_gammq(const double a, const double x)
{
  double gamser,gammcf,gln;
  
  double temp;
  if (x < 0.0 || a <= 0.0) {
    if(a<=0.0) cout << "a<=0 \n";
    mynrerror("Invalid arguments in routine gammq");
  }
  if (x < a+1.0) {
    gser(gamser,a,x,gln);
    temp = 1.0-gamser;
  } else {
    gcf(gammcf,a,x,gln);
    temp =gammcf;
  }

  return exp(gln)*temp;
  
}
//////////////////////////////////////////////////////////////////////////////////////////////////








/* gsl_sf_gamma_inc(k-1/a, z*(1-a)); will overflow for large k.
   k-1/a=170; z*(1-a)=3 4.26907e+304
   k-1/a=171; z*(1-a)=3 7.25742e+306
   gsl: gamma.c:1106: ERROR: overflow
   Default GSL error handler invoked.
   Aborted
*/

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Note that P[k] = C * Get_Gammas_StaticModel(double z, double a, int k) 
// with C = pow(z*(1-a), 1/a)/a;
//
double Get_Gammas_StaticModel(double z, double a, int k)
{
  double P = 0;
  double Q = 0;
  double s = k-1/a;
  double x = z*(1-a);
  if(s<=0)
    P = gsl_sf_gamma_inc(s, x) * exp(-gsl_sf_lngamma(k+1));
  else {
    Q = gsl_sf_gamma_inc_Q(s,x);
    P = Q * exp(gsl_sf_lngamma(s)-gsl_sf_lngamma(k+1));
  }
  //cout << "k-1/a=" << k-1/a << "; z*(1-a)=" << z*(1-a) << "; P=" << P << endl;//' ' << NR_gammq(s, x) << endl; //test  
  return P;
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////





///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Adaptive calculation, depending on k, we can choose calculate P[k]
// according to the exact definition or use approximation : Gamma(k-1/a, z(1-a))/Gamma(k+1) ~= k^{-1-1/a}
double Get_Gammas_StaticModel_Adaptive(double z, double a, int k)
{
  if(k>100000) { 
    //cout << k << endl; //test
    return pow((double)k, -1-1/a);
  }
  else 
    return Get_Gammas_StaticModel(z, a, k);
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////






///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// H7:=  <|k/z - 1|> = 2 \sum_{k=0}^z (1-k/z) P(k)
double Get_SF_SM_AbsoluteDeviation(double z, double a)
{
  double sum = 0.0;
  double C = pow(z*(1-a), 1/a)/a;

  for(int k=0; k<=z ; k++)  {
    double Pk = C*Get_Gammas_StaticModel(z, a, k);
    sum += (1.-k/(double)z) * Pk;
  }
  return 2*sum;
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////




///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// H8:=  - \sum_{k=0}^\infty  P(k) logP(k)
double Get_SF_SM_ShannonEntropy(double z, double a)
{
  double sum  = 0.0;
  double term = DBL_MAX;  // Here DBL_MAX is the max double number
  bool   flag = true;
  int kmax;
  double C = pow(z*(1-a), 1/a)/a;

  for(int k=0; flag ; k++)  {
    double Pk = C*Get_Gammas_StaticModel(z, a, k);
    if(flag) {
      term = - Pk * log(Pk);
      //cout << k << ' ' << term << ' ' << sum << endl; //test
      //if(term>DBL_MIN) // too slow
      //if(fabs(term/sum) > TOLERANCE)
      if(fabs(term/sum) > 1e-8)
	sum += term;
      else {
	flag = false;
	kmax = k;
      }
    }
  }

  //cout << "kmax ="  << kmax << endl;

  return sum;
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////






///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// The following function calculate those concave means, such as <log(k+1)> and <k^0.5> for SF network (static model)
// z = 2m = <k>
// a = 1/(gamma-1)

void Get_SF_SM_concavemean(double z, double a, double& logk1, double& ksr)
{
  double Sum_logk1  = 0.0; // k=0 term
  double Sum_ksr    = 0.0; // k=0 term

  double term_logk1 = DBL_MAX;  // Here DBL_MAX is the max double number
  double term_ksr   = DBL_MAX;  // Here DBL_MAX is the max double number
  
  bool flag_logk1 = true;
  bool flag_ksr   = true;
  
  int kmax_logk1;
  int kmax_ksr;

  double C = pow(z*(1-a), 1/a)/a;

  for(int k=1; (flag_logk1 || flag_ksr); k++)  {
    double Pk = C * Get_Gammas_StaticModel(z, a, k);
    //if(k%100==0) cout << k << ' ' << Pk << endl; //test

    if(flag_logk1) {
      term_logk1 = log(k+1.0) * Pk;
      //if(term_logk1>DBL_MIN)
      //if(term_logk1/Sum_logk1 > TOLERANCE)
      if(term_logk1/Sum_logk1 > 1e-8)
	Sum_logk1 += term_logk1;
      else {
	flag_logk1 = false;
	kmax_logk1 = k;
      }
    }

    if(flag_ksr) {
      term_ksr = sqrt((double)k) * Pk;
      //if(term_ksr>DBL_MIN)
      //if(term_ksr/Sum_ksr > TOLERANCE)
      if(term_ksr/Sum_ksr > 1e-8)
	Sum_ksr += term_ksr;
      else {
	flag_ksr = false;
	kmax_ksr = k;
      }
    }
  }

  cout << "kmax_logk1 = "    << kmax_logk1  << endl;
  cout << "kmax_ksr   = "    << kmax_ksr    << endl << endl;

  logk1 = Sum_logk1;
  ksr   = Sum_ksr;
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////






///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// H9:=  relative mean difference = mean difference/mean
double Get_SF_SM_RMD(double z, double a)
{
  double sum = 0; 
  vector<double> Pk;
  bool flag = true;
  int imax;
  double C = pow(z*(1-a), 1/a)/a;

  for(int i=0; flag; i++) {
    //double Pi = C * Get_Gammas_StaticModel(z, a, i);
    double Pi = C * Get_Gammas_StaticModel_Adaptive(z, a, i);
    Pk.push_back(Pi);
    
    double term = 0;
    for(int j=0; j<i; j++) {
      double Pj = Pk[j];
      term += (i-j)*Pi*Pj;
    }
    
    //if(i==0 || fabs(term/sum) > TOLERANCE)
    //if(i<1000 || term > DBL_MIN) {
    //if(i<1000 || fabs(term/sum) > TOLERANCE) {
    //if(i<1000 || term > TOLERANCE) { // slow
    if(i<10000 || fabs(term/sum) > 1e-8) {
      sum += term;
      if(i%10000==0)
	cout << i << ' ' << term <<' ' << sum << " H9= " << 2*sum/z << endl; //test
    }
    else {
      flag = false;
      imax = i;
    }
  }
  
  double MD = 2 * sum;
  //cout << "<k>= " << z << " MD= " << MD << " with imax = " << imax << endl; //test
  return MD/z;
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////






///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int InverseCDF_SF_SM(double z, double a, double x)
{
  double sum = 0;
  int k;
  double Pk;
  double C = pow(z*(1-a), 1/a)/a;

  for(k=0; sum <= x ; k++)  {
    Pk =  C * Get_Gammas_StaticModel(z, a, k);
    sum += Pk;
  }
  if(sum>x){
    sum -= Pk;
    k--;
  }
  
  //cout << "CDF^(-1) (" << x << ") = " << k << " CDF(" << k << ")=" << sum << endl; 
  return k;
}

// H10:=  Quartile coefficient of dispersion = (Q3-Q1)/(Q3+Q1)
double Get_SF_SM_QCD(double z, double a)
{
  double Q1,Q3;
  Q1 = InverseCDF_SF_SM(z, a, 0.25);
  Q3 = InverseCDF_SF_SM(z, a, 0.75);

  return (Q3-Q1)/(Q3+Q1);
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////








///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
double Get_Pk_ChungLu(vector<double>& w, int k)
{ 
  int N = w.size();
  double sum = 0;
  for(int i=0; i<N; i++) 
    sum += GetPoisson(w[i], k);
  
  return sum/N;
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////




///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// H9:=  relative mean difference = mean difference/mean
double Get_CL_RMD(vector<double>& P, double z)
{
  int kmax = P.size()-1;

  double sum = 0; 
  for(int i=0; i<= kmax; i++) {
    for(int j=0; j<i; j++) 
      sum += (i-j)*P[i]*P[j];
  }
    
  return 2*sum/z;
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////




///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// H9:=  relative mean difference = mean difference/mean
// Note that for exponential digraphs, 
// if P(kin) = P(kout) =  (1-e^{-1/kappa}) e^{-k/kappa}
// then P(k) = (1-e^{-1/kappa})^2 e^{-k/kappa} (k+1)
double Get_EXP_D_RMD(double enkappa, double z)
{
  double sum = 0; 
  vector<double> Pk;
  bool   flag = true;
  int imax;
  double C = (1-enkappa)*(1-enkappa);
  
  double P0 = C;
  Pk.push_back(P0);

  for(int i=1; flag; i++) {
    double Pi =  Pk[i-1]/i * enkappa * (i+1);
    Pk.push_back(Pi);
    
    double term = 0;
    for(int j=0; j<i; j++) {
      double Pj = Pk[j];
      term += (i-j)*Pi*Pj;
    }
    
    if(i==0 || fabs(term/sum) > TOLERANCE) {
      sum += term;
      if(i%10000==0)
	cout << i << ' ' << term <<' ' << sum << endl; //test
    }
    else {
      flag = false;
      imax = i;
    }
  }
  
  double MD = 2 * sum;
  //cout << "<k>= " << z << " MD= " << MD << " with imax = " << imax << endl; //test
  return MD/z;

}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////





///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
double GetPoisson(double lambda, int k)
{
  double P = 1.0/exp(lambda);
  for(int i=1;i<=k;i++)
    P *= (lambda/(double)i);

  return P;
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// The following function calculate those concave means, such as <log(k+1)> and <k^0.5> for ER network
void Get_ER_concavemean(double lambda, double& logk1, double& ksr)
{
  double Sum_logk1  = 0.0; // k=0 term
  double Sum_ksr    = 0.0; // k=0 term

  double term_logk1 = DBL_MAX;  // Here DBL_MAX is the max double number
  double term_ksr   = DBL_MAX;  // Here DBL_MAX is the max double number
  
  bool flag_logk1 = true;
  bool flag_ksr   = true;
  
  int kmax_logk1;
  int kmax_ksr;

  for(int k=1; (flag_logk1 || flag_ksr) ; k++)  {
    double Pk = GetPoisson(lambda, k);

    if(flag_logk1) {
      term_logk1 = log(k+1.0) * Pk;
      if(term_logk1>DBL_MIN)
	Sum_logk1 += term_logk1;
      else {
	flag_logk1 = false;
	kmax_logk1 = k;
      }
    }

    if(flag_ksr) {
      term_ksr = sqrt((double)k) * Pk;
      if(term_ksr>DBL_MIN)
	Sum_ksr += term_ksr;
      else {
	flag_ksr = false;
	kmax_ksr = k;
      }
    }
  }

  cout << "kmax_logk1="  << kmax_logk1 << endl;
  cout << "kmax_ksr="    << kmax_ksr   << endl;

  logk1 = Sum_logk1;
  ksr   = Sum_ksr;
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////




///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// H7:=  <|k/z - 1|> = 2 \sum_{k=0}^z (1-k/z) P(k)
double Get_ER_AbsoluteDeviation(double z)
{
  double sum = 0.0;
  for(int k=0; k<=z ; k++)  {
    double Pk = GetPoisson(z, k);
    sum += (1.-k/(double)z) * Pk;
  }
  return 2*sum;
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// H8:=  - \sum_{k=0}^\infty  P(k) logP(k)
double Get_ER_ShannonEntropy(double z)
{
  double sum  = 0.0;
  double term = DBL_MAX;  // Here DBL_MAX is the max double number
  bool   flag = true;
  int kmax;

  for(int k=0; flag ; k++)  {
    double Pk = GetPoisson(z, k);
    if(flag) {
      term = - Pk * log(Pk);
      //if(term>DBL_MIN)
      if(fabs(term/sum) > TOLERANCE)
	sum += term;
      else {
	flag = false;
	kmax = k;
      }
    }
  }

  //cout << "kmax ="  << kmax << endl;

  return sum;
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////




///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// H9:=  relative mean difference = mean difference/mean
double Get_ER_RMD(double z)
{
  double MD = 2*z * exp(-2*z) * (gsl_sf_bessel_I0(2*z) + gsl_sf_bessel_I1(2*z));
  //cout << "<k>= " << z << " MD= " << MD << endl; //test
  return MD/z;
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////






///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int InverseCDF_ER(double z, double x)
{
  double sum = 0;
  int k;
  double Pk;
  for(k=0; sum <= x ; k++)  {
    Pk = GetPoisson(z, k);
    sum += Pk;
  }
  if(sum>x){
    sum -= Pk;
    k--;
  }
  
  //cout << "CDF^(-1) (" << x << ") = " << k << " CDF(" << k << ")=" << sum << endl; 
  return k;
}

// H10:=  Quartile coefficient of dispersion = (Q3-Q1)/(Q3+Q1)
double Get_ER_QCD(double z)
{
  double Q1,Q3;
  Q1 = InverseCDF_ER(z, 0.25);
  Q3 = InverseCDF_ER(z, 0.75);

  return (Q3-Q1)/(Q3+Q1);
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////





///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// H9:=  relative mean difference = mean difference/mean
double Get_SF_EC_RMD(double gamma, double kappa, double z)
{
  double sum = 0; 
  vector<double> Pk;
  Pk.push_back(0); // Pk=0 for k=0

  bool   flag = true;
  int imax;
  double C = 1.0/polylog(gamma, exp(-1.0/kappa));

  for(int i=1; flag; i++) {
    double Pi = C * pow(i,-gamma) * exp(-i/kappa);
    Pk.push_back(Pi);
    
    double term = 0;
    for(int j=1; j<i; j++) {
      double Pj = Pk[j];
      term += (i-j)*Pi*Pj;
    }
    
    //if(i==0 || fabs(term/sum) > TOLERANCE)
    //if(i<1000 || term > DBL_MIN) {
    if(i<1000 || fabs(term/sum) > TOLERANCE) {
      //if(i<1000 || term > TOLERANCE) { // slow
      //if(i<10000 || fabs(term/sum) > 1e-8) {
      sum += term;
      //if(i%1000==0)
      //cout << i << ' ' << term <<' ' << sum << endl; //test
    }
    else {
      flag = false;
      imax = i;
    }
  }
  
  double MD = 2 * sum;
  //cout << "<k>= " << z << " MD= " << MD << " with imax = " << imax << endl; //test
  return MD/z;
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////




///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// H9:=  relative mean difference = mean difference/mean
// Note that the total degree distribution will be not exactly the same functional form (Poisson is a very special case)
// Actually it should be calculated from  
// P(k) = sum_{kin=0}^k Pin(kin) * Pout(k-kin)  
double Get_SF_EC_D_RMD(double gamma, double kappa, double z)
{
  double sum = 0; 
  vector<double> Pkin; // = Pkout 
  Pkin.push_back(0); // Pk=0 for k=0

  vector<double> Pk; // = Ptotal 
  Pk.push_back(0); // Pk=0 for k=0

  // Don't do this. index will mismatch!! And cause calculated H9 wrong!!
  //Pk.push_back(0); // Pk=0 for k=1 in P(k=1) = Pin(0)Pout(1) + Pin(1)Pout(0) = 0

  bool   flag = true;
  int imax;
  double C = 1.0/polylog(gamma, exp(-1.0/kappa));

  for(int i=1; flag; i++) {
    double Pin_i = C * pow(i,-gamma) * exp(-i/kappa);
    Pkin.push_back(Pin_i);
    
    double P_i = 0;
    for(int kin=0; kin<= i; kin++)
      P_i += Pkin[kin] * Pkin[i-kin];
    Pk.push_back(P_i);

    double term = 0;
    for(int j=1; j<i; j++) {
      double P_j = Pk[j];
      term += (i-j)*P_i*P_j;
    }
    
    //if(i==0 || fabs(term/sum) > TOLERANCE)
    //if(i<1000 || term > DBL_MIN) {
    //if(i<100000 || fabs(term/sum) > TOLERANCE) {
    //if(i<1000 || term > TOLERANCE) { // slow
    if(i<1000 || fabs(term/sum) > 1e-8) {
      sum += term;
      //if(i%100==0)
      //cout << i << ' ' << term <<' ' << sum << ' ' << term/sum << endl; //test
    }
    else {
      flag = false;
      imax = i;
    }
  }
  
  double MD = 2 * sum;
  cout << "<k>= " << z << " MD= " << MD << " with imax = " << imax << endl; //test
  return MD/z;
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////




///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// H9:=  relative mean difference = mean difference/mean
double Get_SF_RMD(double gamma, double z)
{
  double sum = 0; 
  vector<double> Pk;
  Pk.push_back(0); // Pk=0 for k=0

  bool   flag = true;
  int imax;
  double C = 1.0/gsl_sf_zeta(gamma);

  for(int i=1; flag; i++) {
    double Pi = C * pow(i,-gamma);
    Pk.push_back(Pi);
    
    double term = 0;
    for(int j=1; j<i; j++) {
      double Pj = Pk[j];
      term += (i-j)*Pi*Pj;
    }
    
    //if(i==0 || fabs(term/sum) > TOLERANCE)
    //if(i<1000 || term > DBL_MIN) {
    if(i<1000 || fabs(term/sum) > TOLERANCE) {
      //if(i<1000 || term > TOLERANCE) { // slow
      //if(i<10000 || fabs(term/sum) > 1e-8) {
      sum += term;
      //if(i%1000==0)
      //cout << i << ' ' << term <<' ' << sum << endl; //test
    }
    else {
      flag = false;
      imax = i;
    }
  }
  
  double MD = 2 * sum;
  //cout << "<k>= " << z << " MD= " << MD << " with imax = " << imax << endl; //test
  return MD/z;
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////







///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// H9:=  relative mean difference = mean difference/mean for a digraph with the same power-law P(kin) and P(kout)
// Note that the total degree distribution will be not exactly the same functional form (Poisson is a very special case)
// Actually it should be calculated from  
// P(k) = sum_{kin=0}^k Pin(kin) * Pout(k-kin)  

double Get_SF_D_RMD(double gamma, double z) // here z is the mean total degree; gamma is the gamma for P(kin) and P(kout)
{
  double sum = 0; 
  vector<double> Pkin; // = Pkout 
  Pkin.push_back(0); // Pk=0 for k=0

  vector<double> Pk; // = Ptotal 
  Pk.push_back(0); // Pk=0 for k=0

  bool   flag = true;
  int imax;
  double C = 1.0/gsl_sf_zeta(gamma);

  for(int i=1; flag; i++) {
    double Pin_i = C * pow(i,-gamma);
    Pkin.push_back(Pin_i);

    double P_i = 0;
    for(int kin=1; kin<= i-1; kin++)
      P_i += Pkin[kin] * Pkin[i-kin];
    Pk.push_back(P_i);

    double term = 0;
    for(int j=1; j<i; j++) {
      double P_j = Pk[j];
      term += (i-j)*P_i*P_j;
    }
    
    //if(i==0 || fabs(term/sum) > TOLERANCE)
    //if(i<1000 || term > DBL_MIN) {
    if(i<1000 || fabs(term/sum) > TOLERANCE) {
      //if(i<1000 || term > TOLERANCE) { // slow
      //if(i<10000 || fabs(term/sum) > 1e-8) {
      sum += term;
      //if(i%1000==0)
      //cout << i << ' ' << term <<' ' << sum << endl; //test
    }
    else {
      flag = false;
      imax = i;
    }
  }
  
  double MD = 2 * sum;
  //cout << "<k>= " << z << " MD= " << MD << " with imax = " << imax << endl; //test
  return MD/z;
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////







///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///   semi-tensor product  (STP)  /////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//  STP of vectors 
//     X: 1xa    Y: bx1   ==> C:1xn with n=a/b or b/a
// (1) X: 1xnp,  Y: px1   ==> C:1xn
// (2) X: 1xp,   Y: npx1  ==> C:nx1
void STP(double* X, int a, double* Y, int b, double* C)
{
  if(a%b!=0 && b%a!=0) 
    mynrerror("Factor dimension condition (X:1xa, Y:bx1) fails");
  
  int p,np,n;
  if(a>b) { 
    np = a; p = b; n = np/p;
    for(int j=0;j<n;j++) {
      C[j] = 0;
      for(int i=0;i<p;i++) {
	C[j] += X[i*n+j] * Y[i];
      }
    }
  }
  else {
    p = a; np = b; n = np/p;
    for(int j=0;j<n;j++) {
      C[j] = 0;
      for(int i=0;i<p;i++) {
	C[j] += X[i] * Y[i*n+j];
      }
    }
  }
}





//  STP of matrices:  A : m x n, B : p x q ===> C : a x b
// (1) n<p and n*t = p , C : m*t x q   , a=m*t, b=q
// (2) n>p and n = p*t , C : m  x t*q  , a=m, b=t*q
void STP(double** A, int m, int n,  double** B, int p, int q, double** C)
{
  if(n%p!=0 && p%n!=0) 
    mynrerror("Factor dimension condition (A:mxn, B:pxq) fails");

  // i-th row of A:  A[i]
  
  // j-th column of B is given by BT[j]
  double** BT= new double* [q];
  for(int i=0;i<q;i++) {
    BT[i] = new double [p];
    for(int j=0; j<p; j++) 
      BT[i][j] = B[j][i];
  }

  int t, a, b;
  if(n<p) { 
    t = p/n; 
    a = m*t; 
    b = q;
    double* Cij = new double [t];
    for(int i=0; i<m; i++) {
      for(int j=0; j<q; j++) {
	STP(A[i], n, BT[j], p, Cij);
	//cout << "Cij=(";  for(int k=0; k<t; k++)    cout << Cij[k] << ',';  cout << ")\n";
	for(int z=0;z<t;z++)
	  C[t*i+z][j] = Cij[z];
      }
    }
    delete [] Cij;
    
  }
  else    { 
    t = n/p; 
    a = m; 
    b = q*t;
    double* Cij = new double [t];
    for(int i=0; i<m; i++) {
      for(int j=0; j<q; j++) {
	STP(A[i], n, BT[j], p, Cij);
	//cout << "Cij=(";  for(int k=0; k<t; k++)    cout << Cij[k] << ',';  cout << ")\n";
	for(int z=0;z<t;z++)
	  C[i][t*j+z] = Cij[z];
      }
    }
    delete [] Cij;
  }

  
  cout << " = \n";
  for(int i=0;i<a;i++) {
    for(int j=0; j<b; j++) {
      cout << C[i][j] << ' ';
    }
    cout << endl;
  }
  

  for(int i=0;i<q;i++) 
    delete [] BT[i];
  delete [] BT;


}


  /*
  cout << "B= \n";
  for(int i=0;i<p;i++) {
    for(int j=0; j<q; j++) {
      cout << B[i][j] << ' ';
    }
    cout << endl;
  }

  cout << "B^T= \n";
  for(int i=0;i<q;i++) {
    for(int j=0; j<p; j++) {
      cout << BT[i][j] << ' ';
    }
    cout << endl;
  }
  */

void print(double** C, int a, int b) 
{
  cout << "C= \n";
  for(int i=0;i<a;i++) {
    for(int j=0; j<b; j++) {
      cout << C[i][j] << ' ';
    }
    cout << endl;
  }
}

void MatrixPrint(double** X, int a, int b) 
{
  cout << "\n";
  for(int i=1;i<=a;i++) {
    for(int j=1; j<=b; j++) {
      cout.width(15);
      cout << X[i][j] << ' ';
    }
    cout << endl;
  }
}


void MatrixPrint(int** X, int a, int b) 
{
  cout << "\n";
  for(int i=1;i<=a;i++) {
    for(int j=1; j<=b; j++) {
      cout << X[i][j] << ' ';
    }
    cout << endl;
  }
}


//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
double Norm(vector<double>& V)
{
  int N = V.size();
  double sum = 0;
  for(int i=0; i<N; i++) 
    sum += V[i]*V[i];

  return sqrt(sum);
} 


double Sum(vector<double>& V)
{
  int N = V.size();
  double sum = 0;
  for(int i=0; i<N; i++) 
    sum += V[i];

  return sum;
} 


double Error(vector<double>& V1, vector<double>& V2) {
  int N = V1.size();
  double sumerror2 = 0;
  for(int i=0; i<N; i++) {
    double error = V1[i]-V2[i];
    sumerror2 += error*error;
  }
  
  return sqrt(sumerror2);
}



int GetOverlap(set<int>& X, set<int>& Y)
{
  set<int> Z(X.begin(), X.end());
  for(set<int>::iterator p = Y.begin(); p!= Y.end(); p++) 
    Z.insert(*p);

  return (X.size()+Y.size()-Z.size());
}
