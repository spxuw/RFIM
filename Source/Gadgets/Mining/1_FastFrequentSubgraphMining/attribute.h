#ifndef ATTRIBUTES_H
#define ATTRIBUTES_H

#include "common.h"
#include "argraph.h"
#include "argedit.h"

class Point
{ 
public:
	int x;
    Point(int x){ 
		this->x=x;
    }
};

class PointDestroyer: public AttrDestroyer
{ 
public:
	virtual void destroy(void *p){ delete (Point *) p;	}
};

class PointComparator: public AttrComparator
{ 

public:
  PointComparator(){};
  virtual bool compatible(void *pa, void *pb)
    { 
     
		bool f;
#ifdef INDUCEDS
		f =  ( ((Point*)pa)->x == ((Point*)pb)->x );
#else
		f =  (((Point*)pa)->x == EMPTY_EDGE) || ( ((Point*)pa)->x == ((Point*)pb)->x );	
#endif
		return f;
    }
};

/*class EdgeComparator: public AttrComparator
{ 

public:
  EdgeComparator(){};
  virtual bool compatible(void *pa, void *pb)
    { 
     
		int x= ((Point*)pa)->x;
		bool f =  (x == EMPTY_EDGE) || ( x == ((Point*)pb)->x );
		//if (f ) cout <<" matching " << ((Point*)pa)->x  <<((Point*)pb)->x ;
		return f;
    }
};*/

#endif

