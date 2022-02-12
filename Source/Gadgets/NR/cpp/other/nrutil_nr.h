#ifndef _NR_UTIL_H_
#define _NR_UTIL_H_

#include <cstdlib>
#include <string>
#include <cmath>
#include <complex>
#include <iostream>
using namespace std;

typedef double DP;

template<class T>
inline const T SQR(const T a) {return a*a;}

template<class T>
inline const T MAX(const T &a, const T &b)
{return b > a ? (b) : (a);}

inline float MAX(const double &a, const float &b)
{return b > a ? (b) : float(a);}

inline float MAX(const float &a, const double &b)
{return b > a ? float(b) : (a);}

template<class T>
inline const T MIN(const T &a, const T &b)
{return b < a ? (b) : (a);}

inline float MIN(const double &a, const float &b)
{return b < a ? (b) : float(a);}

inline float MIN(const float &a, const double &b)
{return b < a ? float(b) : (a);}

template<class T>
inline const T SIGN(const T &a, const T &b)
{return b >= 0 ? (a >= 0 ? a : -a) : (a >= 0 ? -a : a);}

inline float SIGN(const float &a, const double &b)
{return b >= 0 ? (a >= 0 ? a : -a) : (a >= 0 ? -a : a);}

inline float SIGN(const double &a, const float &b)
{return b >= 0 ? (a >= 0 ? a : -a) : (a >= 0 ? -a : a);}

////////////////////////// added by Yang Liu 11/24/04
template<class T>
inline const int SIGN(const T &a)
{return a >= 0 ? 1 : -1;}

inline int SIGN(const float &a)
{return a >= 0 ? 1 : -1;}

inline int SIGN(const double &a)
{return a >= 0 ? 1 : -1;}
////////////////////////// added by Yang Liu 11/24/04


template<class T>
inline void SWAP(T &a, T &b)
{T dum=a; a=b; b=dum;}

namespace NR {
    inline void nrerror(const string error_text)
	// Numerical Recipes standard error handler
	{
	    cerr << "Numerical Recipes run-time error..." << endl;
	    cerr << error_text << endl;
	    cerr << "...now exiting to system..." << endl;
	    exit(1);
	}
}

using namespace NR;

/////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////
// 1D Vector 
// Modified by Yang Liu 
// add a function: resize 10/07/04
// reload operators 10/23/04

template <class T>
class NRVec 
{
 private:
    int nn;	// size of array. upper index is nn-1
    T *v;
 public:
    NRVec();                    // nn(0), v(0)  
    explicit NRVec(int n);	// Zero-based array
    NRVec(const T &a, int n);	// initialize to constant value
    NRVec(const T *a, int n);	// Initialize to array
    NRVec(const NRVec &rhs);	// Copy constructor
    NRVec & operator=(const NRVec &rhs);	//assignment
    NRVec & operator=(const T &a);	//assign a to every element
        
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////
    NRVec & operator+=(const NRVec &rhs); // v += v2 , added by Yang Liu
    NRVec & operator-=(const NRVec &rhs); // v -= v2 , added by Yang Liu
    NRVec & operator+=(const T &a); // v += a*I , added by Yang Liu
    NRVec & operator-=(const T &a); // v -= a*I , added by Yang Liu


    template<class S>
	friend NRVec<S> operator+(const NRVec<S> &rhs); // v = + v2, added by Yang Liu
    template<class S>
	friend NRVec<S> operator-(const NRVec<S> &rhs); // v = - v2, added by Yang Liu
    template<class S>
	friend NRVec<S> operator+(const NRVec<S> &lhs, const NRVec<S> &rhs); // v = v1 + v2, added by Yang Liu
    template<class S>
	friend NRVec<S> operator-(const NRVec<S> &lhs, const NRVec<S> &rhs); // v = v1 - v2, added by Yang Liu
    template<class S>
	friend NRVec<S> operator+(const NRVec<S> &lhs, const S &a); // v = v1 + a * I, added by Yang Liu
    template<class S>
	friend NRVec<S> operator-(const NRVec<S> &lhs, const S &a); // v = v1 - a * I, added by Yang Liu


    template<class S>
	friend NRVec<S> operator*(const S &a, const NRVec<S> &rhs); // v = scalar * v2, added by Yang Liu
    template<class S>
	friend NRVec<S> operator*(const NRVec<S> &lhs, const S &a); // v = v2 * scalar, added by Yang Liu
    template<class S>
	friend NRVec<S> operator/(const NRVec<S> &lhs, const S &a); // v = v2 / scalar, added by Yang Liu
    template<class S>
	friend NRVec<S> operator*(const NRVec<S> &lhs, const NRVec<S> &rhs); // v = v1 * v2, added by Yang Liu
    template<class S>
	friend ostream& operator<<(ostream&, const NRVec<S> &rhs); 

    /////////////////////////////////////////////////////////////////////////////////////////////////////////////

    inline T & operator[](const int i);	//i'th element
    inline const T & operator[](const int i) const;
    inline int size() const;

    /////////////////////////////////////////////////////////////////////////////////////////////////////////////
    inline void resize(int n);           // resize the vector, added by Yang Liu

    template <class S>
      friend S dot (const NRVec<S> &lhs, const NRVec<S> &rhs); // scalar = v1 * v2, added by Yang Liu

    template <class S>
      friend S normalize (NRVec<S> &v); // normalize v, and return the old norm, added by Yang Liu
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////
  

    ~NRVec();
};






template <class T>
NRVec<T>::NRVec() : nn(0), v(0) {}

template <class T>
NRVec<T>::NRVec(int n) : nn(n), v(new T[n]) 
{
     for(int i=0; i<n; i++)
	v[i] = 0;
}

template <class T>
NRVec<T>::NRVec(const T &a, int n) : nn(n), v(new T[n])
{
    for(int i=0; i<n; i++)
	v[i] = a;
}

template <class T>
NRVec<T>::NRVec(const T *a, int n) : nn(n), v(new T[n])
{
    for(int i=0; i<n; i++)
	v[i] = *a++;
}

template <class T>
NRVec<T>::NRVec(const NRVec<T> &rhs) : nn(rhs.nn), v(new T[nn])
{
    for(int i=0; i<nn; i++)
	v[i] = rhs[i];
}

template <class T>
NRVec<T> & NRVec<T>::operator=(const NRVec<T> &rhs)
// postcondition: normal assignment via copying has been performed;
//		if vector and rhs were different sizes, vector
//		has been resized to match the size of rhs
{
    if (this != &rhs)
    {
	if (nn != rhs.nn) {
	    if (v != 0) delete [] (v);
	    nn=rhs.nn;
	    v= new T[nn];
	}
	for (int i=0; i<nn; i++)
	    v[i]=rhs[i];
    }
    return *this;
}

template <class T>
NRVec<T> & NRVec<T>::operator=(const T &a)	//assign a to every element
{
    for (int i=0; i<nn; i++)
	v[i]=a;
    return *this;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////// added by Yang Liu   10/23/04  /////////////////////////////////////////////////////////

///////////////      v += v2     /////////////////////////////////
template <class T>
NRVec<T> & NRVec<T>::operator+=(const NRVec<T> &rhs) 
{
    if (this != &rhs)
    {
	if (nn != rhs.nn) {
	    if (v != 0) delete [] (v);
	    nn=rhs.nn;
	    v= new T[nn];
	}
	for (int i=0; i<nn; i++)
	    v[i] += rhs[i];
    }
    return *this;
    
}

///////////////      v -= v2     /////////////////////////////////
template <class T>
NRVec<T> & NRVec<T>::operator-=(const NRVec<T> &rhs) 
{
    if (this != &rhs)
    {
	if (nn != rhs.nn) {
	    if (v != 0) delete [] (v);
	    nn=rhs.nn;
	    v= new T[nn];
	}
	for (int i=0; i<nn; i++)
	    v[i] -= rhs[i];
    }
    return *this;
}



///////////////      v += a*I     /////////////////////////////////
template <class T>
NRVec<T> & NRVec<T>::operator+=(const T &a) 
{
    for (int i=0; i<nn; i++)
	v[i] += a;

    return *this;

}


///////////////      v -= a*I     /////////////////////////////////
template <class T>
NRVec<T> & NRVec<T>::operator-=(const T &a) 
{
    for (int i=0; i<nn; i++)
	v[i] -= a;

    return *this;

}


///////////////      v = +v2     /////////////////////////////////
template <class T>
inline NRVec<T> operator+(const NRVec<T> &rhs) 
{
    return rhs;
}

///////////////      v = -v2     /////////////////////////////////
template <class T>
inline NRVec<T> operator-(const NRVec<T> &rhs)  
{
    return NRVec<T>(rhs.nn) - rhs;
}


///////////////      v = v1 + v2     /////////////////////////////////
template <class T>
NRVec<T> operator+(const NRVec<T> &lhs, const NRVec<T> &rhs)  
{
    if(lhs.nn != rhs.nn) nrerror("bad vector sizes in operation: v1+v2");
    NRVec<T> sum = lhs;
    sum += rhs;
    return sum;
}

///////////////      v = v1 - v2     /////////////////////////////////
template <class T>
NRVec<T> operator-(const NRVec<T> &lhs, const NRVec<T> &rhs)  
{
    if(lhs.nn != rhs.nn) nrerror("bad vector sizes in operation: v1-v2");
    NRVec<T> sum = lhs;
    sum -= rhs;
    return sum;
}


///////////////      v = v1 + a*I     /////////////////////////////////
template <class T>
NRVec<T> operator+(const NRVec<T> &lhs, const T &a) 
{
    NRVec<T> sum = lhs;
    sum += a;
    return sum;
}


///////////////      v = v1 - a*I     /////////////////////////////////
template <class T>
NRVec<T> operator-(const NRVec<T> &lhs, const T &a) 
{
    NRVec<T> sum = lhs;
    sum -= a;
    return sum;
}


///////////////      v = scalar * v2     /////////////////////////////////
template <class T>
NRVec<T> operator*(const T &scalar, const NRVec<T> &rhs)  
{
    NRVec<T> tmp(rhs.nn);
    for(int i=0; i<rhs.nn; i++)
	tmp[i] = scalar * rhs[i];
    return tmp;
}

///////////////      v = v2 * scalar     /////////////////////////////////
template <class T>
NRVec<T> operator*(const NRVec<T> &lhs, const T &scalar)  
{
    return scalar*lhs;
}

///////////////      v = v2 / scalar     /////////////////////////////////
template <class T>
NRVec<T> operator/(const NRVec<T> &lhs, const T &scalar)  
{
    if(!scalar)
	nrerror("division by zero in vector-scalar division");
    return (1.0/scalar)*lhs;
}

///////////////      v = v1 * v2       /////////////////////////////////
template <class T>
NRVec<T> operator*(const NRVec<T> &lhs, const NRVec<T> &rhs) 
{
    if(lhs.nn != rhs.nn)
	nrerror("bad vector size in operation v = v1 * v2");
    int n = lhs.nn;
    NRVec<T> tmp(n);

    for(int i=0; i<n; i++)
	tmp[i] = lhs[i]*rhs[i];

    return tmp;
}

///////////////      out <<  v        /////////////////////////////////
template <class T>
ostream& operator<<(ostream& s, const NRVec<T> &rhs)  
{
    for(int i=0; i<rhs.size(); i++)
    {
	s << rhs[i] << ' ' ;
    }
    return s;

}

///////////////      scalar = v1 dot v2    /////////////////////////////////
template <class T>
T dot(const NRVec<T> &lhs, const NRVec<T> &rhs)  
{
    if(lhs.nn != rhs.nn)
	nrerror("bad vector size in dot product scalar = v1 . v2");
    T tmp = lhs[0]*rhs[0];

    for(int i=1; i<lhs.nn; i++)
	tmp += lhs[i]*rhs[i];
    return tmp;
}
///////////////////////////////////////////////////////////////////////////////////


///////////////      normalize v, return the old Euclidean norm    //////////////////////////////////
template <class T>
T normalize(NRVec<T> &v)  
{
  T norm = 0;
  for(int i=0; i<v.nn; i++)
    norm += v[i]*v[i];
  
  norm = sqrt(norm);
  for(int i=0; i<v.nn; i++)
    v[i] /= norm;

  return norm;
}
///////////////////////////////////////////////////////////////////////////////////




template <class T>
inline T & NRVec<T>::operator[](const int i)	//subscripting
{
    return v[i];
}

template <class T>
inline const T & NRVec<T>::operator[](const int i) const	//subscripting
{
    return v[i];
}

template <class T>
inline int NRVec<T>::size() const
{
    return nn;
}

//////////////// added by Yang Liu ////////////////
template <class T>
inline void NRVec<T>::resize(const int n) 
{
    nn = n;
    v = new T[n];
}
////////////////////////////////////////////////////


template <class T>
NRVec<T>::~NRVec()
{
    if (v != 0)
	delete[] (v);
}












/////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////
// 2D Matrix 
// Modified by Yang Liu 
// add a function: resize 10/07/04
// reload operators 10/23/04
template <class T>
class NRMat {
 private:
    int nn;
    int mm;
    T **v;
 public:
    NRMat();                            // nn(0),mm(0),v(0)
    NRMat(int n, int m);		// Zero-based array
    NRMat(const T &a, int n, int m);	//Initialize to constant
    NRMat(const T *a, int n, int m);	// Initialize to 1D array
//    NRMat(const T **a, int n, int m);	// Initialize to 2D array
    NRMat(const NRMat &rhs);		// Copy constructor
    NRMat & operator=(const NRMat &rhs);	//assignment
    NRMat & operator=(const T &a);		//assign a to every element

    ///////////////////////////////////////////////////////////////////////////////
    NRMat & operator+=(const NRMat &rhs); // A += A2 , added by Yang Liu
    NRMat & operator-=(const NRMat &rhs); // A -= A2 , added by Yang Liu
    NRMat & operator+=(const T &a); // A += a*I , added by Yang Liu
    NRMat & operator-=(const T &a); // A -= a*I , added by Yang Liu


    template<class S>
	friend NRMat<S> operator+(const NRMat<S> &rhs); // A = + A2, added by Yang Liu
    template<class S>
	friend NRMat<S> operator-(const NRMat<S> &rhs); // A = - A2, added by Yang Liu
    template<class S>
	friend NRMat<S> operator+(const NRMat<S> &lhs, const NRMat<S> &rhs); // A = A1 + A2, added by Yang Liu
    template<class S>
	friend NRMat<S> operator-(const NRMat<S> &lhs, const NRMat<S> &rhs); // A = A1 - A2, added by Yang Liu
    template<class S>
	friend NRMat<S> operator+(const NRMat<S> &lhs, const S &a); // A = A1 + a*I, added by Yang Liu
    template<class S>
	friend NRMat<S> operator-(const NRMat<S> &lhs, const S &a); // A = A1 - a*I, added by Yang Liu


    template<class S>
	friend NRMat<S> operator*(const S &a, const NRMat<S> &rhs); // A = scalar * A2, added by Yang Liu
    template<class S>
	friend NRMat<S> operator*(const NRMat<S> &lhs, const S &a); // A = A2 * scalar, added by Yang Liu
    template<class S>
	friend NRMat<S> operator/(const NRMat<S> &lhs, const S &a); // A = A2 / scalar, added by Yang Liu
    template<class S>
	friend NRMat<S> operator*(const NRMat<S> &lhs, const NRMat<S> &rhs); // A = A1 * A2, added by Yang Liu
    template<class S>
	friend NRVec<S> operator*(const NRMat<S> &lhs, const NRVec<S> &rhs); // A = A1 * v2, added by Yang Liu

    template<class S>
	friend ostream& operator<<(ostream&, const NRMat<S> &rhs); 

  /////////////////////////////////////////////////////////////////////////////////

    inline T* operator[](const int i);	//subscripting: pointer to row i
    inline const T* operator[](const int i) const;
    inline int nrows() const;
    inline int ncols() const;

    /////////////////////////////////////////////////////////////////////////////
    inline void resize(int n, int m);       // resize the matrix, added by Yang Liu
    /////////////////////////////////////////////////////////////////////////////

    ~NRMat();
};

template <class T>
NRMat<T>::NRMat() : nn(0), mm(0), v(0) {}

template <class T>
NRMat<T>::NRMat(int n, int m) : nn(n), mm(m), v(new T*[n])
{
    v[0] = new T[m*n];
    for (int i=1; i< n; i++)
	v[i] = v[i-1] + m;
    
}

template <class T>
NRMat<T>::NRMat(const T &a, int n, int m) : nn(n), mm(m), v(new T*[n])
{
    int i,j;
    v[0] = new T[m*n];
    for (i=1; i< n; i++)
	v[i] = v[i-1] + m;
    for (i=0; i< n; i++)
	for (j=0; j<m; j++)
	    v[i][j] = a;
}

template <class T>
NRMat<T>::NRMat(const T *a, int n, int m) : nn(n), mm(m), v(new T*[n])
{
    int i,j;
    v[0] = new T[m*n];
    for (i=1; i< n; i++)
	v[i] = v[i-1] + m;
    for (i=0; i< n; i++)
	for (j=0; j<m; j++)
	    v[i][j] = *a++;
}


/*// added by Yang Liu 10/23/04
//
// But it doesn't work for the definition of NRMat from a 2D array
// Actually, the definition of NRMat from a 2D array can be performed
// by using:
//  NRMat A(2Darray[0], n ,m)
// because 2Darra[0] stands for the address of the [0][0] element
// in the 2Darray. Similarly, we can use
// NRMat A(*2Darray, n, m); 
//
// The easiest way is to write the 2Darray in a 1D array, then
//  NRMat A(1Darray, n, m); 
//
template <class T>
NRMat<T>::NRMat(const T **a, int n, int m) : nn(n), mm(m), v(new T*[n])
{
    int i,j;
    v[0] = new T[m*n];
    for (i=1; i< n; i++)
	v[i] = v[i-1] + m;

    for (i=0; i< n; i++)
	for (j=0; j<m; j++)
	    v[i][j] = a[i][j];
}
*/

template <class T>
NRMat<T>::NRMat(const NRMat &rhs) : nn(rhs.nn), mm(rhs.mm), v(new T*[nn])
{
    int i,j;
    v[0] = new T[mm*nn];
    for (i=1; i< nn; i++)
	v[i] = v[i-1] + mm;

    for (i=0; i< nn; i++)
	for (j=0; j<mm; j++)
	    v[i][j] = rhs[i][j];
}

template <class T>
NRMat<T> & NRMat<T>::operator=(const NRMat<T> &rhs)
// postcondition: normal assignment via copying has been performed;
//		if matrix and rhs were different sizes, matrix
//		has been resized to match the size of rhs
{
    if (this != &rhs) 
    {
	int i,j;
	if (nn != rhs.nn || mm != rhs.mm) 
	{
	    if (v != 0) 
	    {
		delete[] (v[0]);
		delete[] (v);
	    }
	    nn=rhs.nn;
	    mm=rhs.mm;
	    v = new T*[nn];
	    v[0] = new T[mm*nn];
	    for (i=1; i< nn; i++)
		v[i] = v[i-1] + mm;
	}

	for (i=0; i< nn; i++)
	    for (j=0; j<mm; j++)
		v[i][j] = rhs[i][j];
    }
    return *this;
}

template <class T>
NRMat<T> & NRMat<T>::operator=(const T &a)	//assign a to every element
{
    for (int i=0; i< nn; i++)
	for (int j=0; j<mm; j++)
	    v[i][j] = a;
    return *this;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////// added by Yang Liu   10/23/04  /////////////////////////////////////////////////////////

///////////////      A += A2     /////////////////////////////////
template <class T>
NRMat<T> & NRMat<T>::operator+=(const NRMat<T> &rhs) 
{

    if (this != &rhs) 
    {
	int i,j;
	if (nn != rhs.nn || mm != rhs.mm) 
	{
	    if (v != 0) 
	    {
		delete[] (v[0]);
		delete[] (v);
	    }
	    nn=rhs.nn;
	    mm=rhs.mm;
	    v = new T*[nn];
	    v[0] = new T[mm*nn];
	}
	for (i=1; i< nn; i++)
	    v[i] = v[i-1] + mm;

	for (i=0; i< nn; i++)
	    for (j=0; j<mm; j++)
		v[i][j] += rhs[i][j];
    }
    return *this;
}    

///////////////      A -= A2     /////////////////////////////////
template <class T>
NRMat<T> & NRMat<T>::operator-=(const NRMat<T> &rhs) 
{

    if (this != &rhs) 
    {
	int i,j;
	if (nn != rhs.nn || mm != rhs.mm) 
	{
	    if (v != 0) 
	    {
		delete[] (v[0]);
		delete[] (v);
	    }
	    nn=rhs.nn;
	    mm=rhs.mm;
	    v = new T*[nn];
	    v[0] = new T[mm*nn];
	}
	for (i=1; i< nn; i++)
	    v[i] = v[i-1] + mm;

	for (i=0; i< nn; i++)
	    for (j=0; j<mm; j++)
		v[i][j] -= rhs[i][j];
    }
    return *this;
}    



///////////////      A += a*I     /////////////////////////////////
template <class T>
NRMat<T> & NRMat<T>::operator+=(const T &a) 
{
    for (int i=0; i< nn; i++)
	for (int j=0; j<mm; j++)
	    v[i][j] += a;

    return *this;
}    

///////////////      A -= a*I     /////////////////////////////////
template <class T>
NRMat<T> & NRMat<T>::operator-=(const T &a) 
{
    for (int i=0; i< nn; i++)
	for (int j=0; j<mm; j++)
	    v[i][j] -= a;

    return *this;
}    


///////////////      A = +A2     /////////////////////////////////
template <class T>
inline NRMat<T> operator+(const NRMat<T> &rhs) 
{
    return rhs;
}    


///////////////      A = -A2     /////////////////////////////////
template <class T>
inline NRMat<T> operator-(const NRMat<T> &rhs) 
{
    return NRMat<T>(rhs.nn,rhs.mm) - rhs;
}    


///////////////      A = A1+A2     /////////////////////////////////
template <class T>
NRMat<T> operator+(const NRMat<T> &lhs, const NRMat<T> &rhs)  
{
    if(lhs.nn != rhs.nn || lhs.mm != rhs.mm) 
	nrerror("bad matrix sizes in operation: A = B + C");
    NRMat<T> sum = lhs;
    sum += rhs;
    return sum;
}

///////////////      A = A1-A2     /////////////////////////////////
template <class T>
NRMat<T> operator-(const NRMat<T> &lhs, const NRMat<T> &rhs)  
{
    if(lhs.nn != rhs.nn  || lhs.mm != rhs.mm) 
	nrerror("bad matrix sizes in operation: A = B - C");
    NRMat<T> sum = lhs;
    sum -= rhs;
    return sum;
}

///////////////      A = A1 + a*I     /////////////////////////////////
template <class T>
NRMat<T> operator+(const NRMat<T> &lhs, const T &a)  
{
    NRMat<T> sum = lhs;
    sum += a;
    return sum;
}

///////////////      A = A1 - a*I     /////////////////////////////////
template <class T>
NRMat<T> operator-(const NRMat<T> &lhs, const T &a)  
{
    NRMat<T> sum = lhs;
    sum -= a;
    return sum;
}

///////////////      A = a * A2     /////////////////////////////////
template <class T>
NRMat<T> operator*(const T &a, const NRMat<T> &rhs)  
{
    NRMat<T> tmp(rhs.nn,rhs.mm);
    for(int i=0; i<rhs.nn; i++)
	for(int j=0; j<rhs.mm; j++)
	    tmp[i][j] = a * rhs[i][j];
    return tmp;
}

///////////////      A = A2 * a      /////////////////////////////////
template <class T>
NRMat<T> operator*(const NRMat<T> &lhs, const T &a)  
{
    return a*lhs;
}


///////////////      A = A2 / a      /////////////////////////////////
template <class T>
NRMat<T> operator/(const NRMat<T> &lhs, const T &a) 
{
    if(!a)
	nrerror("division by zero in matrix-scalar division");
    return (1.0/a)*lhs;
}

///////////////      A = A1 * A2      /////////////////////////////////
template <class T>
NRMat<T> operator*(const NRMat<T> &lhs, const NRMat<T> &rhs)  
{
    if(lhs.nn != rhs.nn || lhs.mm != rhs.mm)
	nrerror("bad matrix size in operation: C = A * B");
    int n = lhs.nn;
    int m = lhs.mm;
    NRMat<T> tmp(n,m);

    for(int i=0; i<n; i++)
    {
	for(int j=0; j<m; j++)
	{
	    tmp[i][j] = 0;
	    for(int k=0; k<m; k++)
	    {
		tmp[i][j] += lhs[i][k]*rhs[k][j];
	    }
	}
    }

    return tmp;
}

///////////////      A = A1 * v2      /////////////////////////////////
template <class T>
NRVec<T> operator*(const NRMat<T> &lhs, const NRVec<T> &rhs)  // A = A1 * v2
{
    if(lhs.mm != rhs.size()) // we have to use rhs.size(), since rhs.nn is private member
	nrerror("bad vector size in operation: C = A * v2");
    int n = lhs.nn;
    int m = lhs.mm;
    NRVec<T> tmp(n);

    for(int i=0; i<n; i++)
	for(int j=0; j<m; j++)
	    tmp[i] += lhs[i][j]*rhs[j];

    return tmp;
}

///////////////      out << A       /////////////////////////////////
template <class T>
ostream& operator<<(ostream& s, const NRMat<T> &rhs)  
{
    for(int i=0; i<rhs.nrows(); i++)
    {
	for(int j=0; j<rhs.ncols(); j++)
	{
	    cout.width(8);
	    s << rhs[i][j] << ' ';
	}
	s << endl;
    }
    return s;

}

/////////////////////////////////////////////////////////////////////////////////////////////////////////




template <class T>
inline T* NRMat<T>::operator[](const int i)	//subscripting: pointer to row i
{
    return v[i];
}

template <class T>
inline const T* NRMat<T>::operator[](const int i) const
{
    return v[i];
}

template <class T>
inline int NRMat<T>::nrows() const
{
    return nn;
}

template <class T>
inline int NRMat<T>::ncols() const
{
    return mm;
}


//////////////// added by Yang Liu ////////////////
// Note that the allocation is a little different from
// the standard method. Because the descontruction function 
// is designed in a different way. I mean the way of free the 
// allocated memory space. 
//
// This kind of difference can also been seen in the 
// class NRMat3D. I designed the resize function in the same taste,
// similar to its own construction function.

template <class T>
inline void NRMat<T>::resize(const int n, const int m) 
{
    nn = n;
    mm = m;
    v = new T* [n];
    v[0] = new T [m*n];
    for (int i=1; i< n; i++)
	v[i] = v[i-1] + m;

}
////////////////////////////////////////////////////



template <class T>
NRMat<T>::~NRMat()
{
    if (v != 0) {
	delete[] (v[0]);
	delete[] (v);
    }
}


















/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
// 3D Matrix
// Modified by Yang Liu 10/23/04
// add a function: resize
// add copy constructor adn reload operators = <<  10/24/04
// I didn't reload the operators, such as * += -= / , for the 3D matrix
// because I think it is not efficient at all in 3D. 
template <class T>
class NRMat3d {
 private:
    int nn;
    int mm;
    int kk;
    T ***v;
 public:
    NRMat3d();
    NRMat3d(int n, int m, int k);

    ////  added by Yang Liu 10/24/04   //////////////////////////////////////////////////////////////////////////
    NRMat3d(const T &a, int n, int m, int k);	//Initialize to constant
    NRMat3d(const T *a, int n, int m, int k);	// Initialize to 1D array
    NRMat3d(const NRMat3d &rhs);		// Copy constructor
    NRMat3d & operator=(const NRMat3d &rhs);	//assignment
    NRMat3d & operator=(const T &a);		//assign a to every element
    template<class S>
	friend ostream& operator<<(ostream&, const NRMat3d<S> &rhs); 
    /////////////////////////////////////////////////////////////////////////////


    inline T** operator[](const int i);	//subscripting: pointer to row i
    inline const T* const * operator[](const int i) const;
    inline int dim1() const;
    inline int dim2() const;
    inline int dim3() const;

    /////////////////////////////////////////////////////////////////////////////
    inline void resize(int n, int m, int k);       // resize the 3Dmatrix, added by Yang Liu
    /////////////////////////////////////////////////////////////////////////////

    ~NRMat3d();
};

template <class T>
NRMat3d<T>::NRMat3d(): nn(0), mm(0), kk(0), v(0) {}

template <class T>
NRMat3d<T>::NRMat3d(int n, int m, int k) : nn(n), mm(m), kk(k), v(new T**[n])
{
    int i,j;
    v[0] = new T*[n*m];
    v[0][0] = new T[n*m*k];
    for(j=1; j<m; j++)
	v[0][j] = v[0][j-1] + k;
    for(i=1; i<n; i++) 
    {
	v[i] = v[i-1] + m;
	v[i][0] = v[i-1][0] + m*k;
	for(j=1; j<m; j++)
	    v[i][j] = v[i][j-1] + k;
    }
}


template <class T>
NRMat3d<T>::NRMat3d(const T &a, int n, int m, int k) : nn(n), mm(m), kk(k), v(new T**[n])
{
    int i,j,s;
    v[0] = new T*[n*m];
    v[0][0] = new T[n*m*k];
    for(j=1; j<m; j++)
	v[0][j] = v[0][j-1] + k;
    for(i=1; i<n; i++) 
    {
	v[i] = v[i-1] + m;
	v[i][0] = v[i-1][0] + m*k;
	for(j=1; j<m; j++)
	    v[i][j] = v[i][j-1] + k;
    }
   
    for (i=0; i< n; i++)
	for (j=0; j<m; j++)
	    for (s=0; s<k; s++)
		v[i][j][s] = a;
}

template <class T>
NRMat3d<T>::NRMat3d(const T *a, int n, int m, int k) : nn(n), mm(m), kk(k), v(new T**[n])
{
    int i,j,s;
    v[0] = new T*[n*m];
    v[0][0] = new T[n*m*k];
    for(j=1; j<m; j++)
	v[0][j] = v[0][j-1] + k;
    for(i=1; i<n; i++) 
    {
	v[i] = v[i-1] + m;
	v[i][0] = v[i-1][0] + m*k;
	for(j=1; j<m; j++)
	    v[i][j] = v[i][j-1] + k;
    }
   
    for (i=0; i< n; i++)
	for (j=0; j<m; j++)
	    for (s=0; s<k; s++)
		v[i][j][s] = *a++;
}


/////////////////////////////////////////////////////////////
template <class T>
NRMat3d<T>::NRMat3d(const NRMat3d &rhs) : nn(rhs.nn), mm(rhs.mm), kk(rhs.kk), v(new T**[nn])
{
    int i,j,k;
    v[0] = new T*[nn*mm];
    v[0][0] = new T[nn*mm*kk];
    for(j=1; j<mm; j++)
	v[0][j] = v[0][j-1] + kk;
    for(i=1; i<nn; i++) 
    {
	v[i] = v[i-1] + mm;
	v[i][0] = v[i-1][0] + mm*kk;
	for(j=1; j<mm; j++)
	    v[i][j] = v[i][j-1] + kk;
    }

    for (i=0; i< nn; i++)
	for (j=0; j<mm; j++)
	    for (k=0; k<kk; k++)
		v[i][j][k] = rhs[i][j][k];

}

template <class T>
NRMat3d<T> & NRMat3d<T>::operator=(const NRMat3d<T> &rhs)
// postcondition: normal assignment via copying has been performed;
//		if matrix and rhs were different sizes, matrix
//		has been resized to match the size of rhs
{
    if (this != &rhs) 
    {
	int i,j,k;
	if (nn != rhs.nn || mm != rhs.mm || kk != rhs.kk) 
	{
	    if (v != 0) 
	    {
		delete[] (v[0][0]);
		delete[] (v[0]);
		delete[] (v);
	    }
	    nn=rhs.nn;
	    mm=rhs.mm;
	    kk=rhs.kk;

	    v = new T**[nn];
   	    v[0] = new T*[nn*mm];
	    v[0][0] = new T[nn*mm*kk];
	    for(j=1; j<mm; j++)
		v[0][j] = v[0][j-1] + kk;
	    for(i=1; i<nn; i++) 
	    {
		v[i] = v[i-1] + mm;
		v[i][0] = v[i-1][0] + mm*kk;
		for(j=1; j<mm; j++)
		    v[i][j] = v[i][j-1] + kk;
	    }
	}

	for (i=0; i< nn; i++)
	    for (j=0; j<mm; j++)
		for (k=0; k<kk; k++)
		    v[i][j][k] = rhs[i][j][k];
	
    }
    return *this;
}

template <class T>
NRMat3d<T> & NRMat3d<T>::operator=(const T &a)	//assign a to every element
{
    for (int i=0; i< nn; i++)
	for (int j=0; j<mm; j++)
	    for (int k=0; k<kk; k++)
		v[i][j][k] = a;
    return *this;
}



///////////////      out << A       /////////////////////////////////
template <class T>
ostream& operator<<(ostream& s, const NRMat3d<T> &rhs)  
{
    for(int i=0; i<rhs.dim1(); i++)
    {
	for(int j=0; j<rhs.dim2(); j++)
	    {
		for(int k=0; k<rhs.dim3(); k++)
		    s << rhs[i][j][k] << ' ';
		s << endl;
	    }
	s << endl;
    }
    return s;

}

/////////////////////////////////////////////////////////////////////////////////////////////////////////









//////////////////////////////////////////////////////
template <class T>
inline T** NRMat3d<T>::operator[](const int i) //subscripting: pointer to row i
{
    return v[i];
}

template <class T>
inline const T* const * NRMat3d<T>::operator[](const int i) const
{
    return v[i];
}

template <class T>
inline int NRMat3d<T>::dim1() const
{
    return nn;
}

template <class T>
inline int NRMat3d<T>::dim2() const
{
    return mm;
}

template <class T>
inline int NRMat3d<T>::dim3() const
{
    return kk;
}


//////////////// added by Yang Liu ////////////////
template <class T>
inline void NRMat3d<T>::resize(const int n, const int m, const int k) 
{
    nn = n;
    mm = m;
    kk = k;
    v = new T** [n];
	
    // Anormal allocation method
    int i,j;
    v[0] = new T*[n*m];
    v[0][0] = new T[n*m*k];
    for(j=1; j<m; j++)
	v[0][j] = v[0][j-1] + k;
    for(i=1; i<n; i++) {
	v[i] = v[i-1] + m;
	v[i][0] = v[i-1][0] + m*k;
	for(j=1; j<m; j++)
	    v[i][j] = v[i][j-1] + k;
    }

    /*for (int i=0; i<nn;i++)
      {
      v[i] = new T* [mm];
      for (int j=0; j<mm;j++)
      v[i][j] = new T [kk];
      }*/ //standard allocation method
	
}


	
////////////////////////////////////////////////////




template <class T>
NRMat3d<T>::~NRMat3d()
{
    if (v != 0) {
	delete[] (v[0][0]);
	delete[] (v[0]);
	delete[] (v);
    }
}





//The next 3 classes are used in artihmetic coding, Huffman coding, and
//wavelet transforms respectively. This is as good a place as any to put them!

class arithcode {
 private:
    NRVec<unsigned long> *ilob_p,*iupb_p,*ncumfq_p;
 public:
    NRVec<unsigned long> &ilob,&iupb,&ncumfq;
    unsigned long jdif,nc,minint,nch,ncum,nrad;
    arithcode(unsigned long n1, unsigned long n2, unsigned long n3)
	: ilob_p(new NRVec<unsigned long>(n1)),
	iupb_p(new NRVec<unsigned long>(n2)),
	ncumfq_p(new NRVec<unsigned long>(n3)),
	ilob(*ilob_p),iupb(*iupb_p),ncumfq(*ncumfq_p) {}
    ~arithcode() {
	if (ilob_p != 0) delete ilob_p;
	if (iupb_p != 0) delete iupb_p;
	if (ncumfq_p != 0) delete ncumfq_p;
    }
};

class huffcode {
 private:
    NRVec<unsigned long> *icod_p,*ncod_p,*left_p,*right_p;
 public:
    NRVec<unsigned long> &icod,&ncod,&left,&right;
    int nch,nodemax;
    huffcode(unsigned long n1, unsigned long n2, unsigned long n3,
	     unsigned long n4) :
	icod_p(new NRVec<unsigned long>(n1)),
	ncod_p(new NRVec<unsigned long>(n2)),
	left_p(new NRVec<unsigned long>(n3)),
	right_p(new NRVec<unsigned long>(n4)),
	icod(*icod_p),ncod(*ncod_p),left(*left_p),right(*right_p) {}
    ~huffcode() {
	if (icod_p != 0) delete icod_p;
	if (ncod_p != 0) delete ncod_p;
	if (left_p != 0) delete left_p;
	if (right_p != 0) delete right_p;
    }
};

/*
class wavefilt {
 private:
    NRVec<DP> *cc_p,*cr_p;
 public:
    int ncof,ioff,joff;
    NRVec<DP> &cc,&cr;
    wavefilt() : cc(*cc_p),cr(*cr_p) {}
    wavefilt(const DP *a, const int n) :  //initialize to array
	cc_p(new NRVec<DP>(n)),cr_p(new NRVec<DP>(n)),
	ncof(n),ioff(-(n >> 1)),joff(-(n >> 1)),cc(*cc_p),cr(*cr_p) {
	int i;
	for (i=0; i<n; i++)
	    cc[i] = *a++;
	DP sig = -1.0;
	for (i=0; i<n; i++) {
	    cr[n-1-i]=sig*cc[i];
	    sig = -sig;
	}
    }
    ~wavefilt() {
	if (cc_p != 0) delete cc_p;
	if (cr_p != 0) delete cr_p;
    }
};
*/

//Overloaded complex operations to handle mixed float and double
//This takes care of e.g. 1.0/z, z complex<float>

inline const complex<float> operator+(const double &a,
				      const complex<float> &b) { return float(a)+b; }

inline const complex<float> operator+(const complex<float> &a,
				      const double &b) { return a+float(b); }

inline const complex<float> operator-(const double &a,
				      const complex<float> &b) { return float(a)-b; }

inline const complex<float> operator-(const complex<float> &a,
				      const double &b) { return a-float(b); }

inline const complex<float> operator*(const double &a,
				      const complex<float> &b) { return float(a)*b; }

inline const complex<float> operator*(const complex<float> &a,
				      const double &b) { return a*float(b); }

inline const complex<float> operator/(const double &a,
				      const complex<float> &b) { return float(a)/b; }

inline const complex<float> operator/(const complex<float> &a,
				      const double &b) { return a/float(b); }

//some compilers choke on pow(float,double) in single precision. also atan2

inline float pow (float x, double y) {return pow(double(x),y);}
inline float pow (double x, float y) {return pow(x,double(y));}
inline float atan2 (float x, double y) {return atan2(double(x),y);}
inline float atan2 (double x, float y) {return atan2(x,double(y));}
#endif /* _NR_UTIL_H_ */
