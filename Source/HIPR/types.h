#ifndef TYPES_H
#define TYPES_H
#include <vector>
#include <iostream>
using namespace std;


typedef vector<double>  Vec;         // for h[i]
typedef long long      hType;        // for RESOLUTION, J, CC, heff[i], Wi, // Yang Liu 06/19/06  
typedef unsigned long long cType;    // for capacity // Yang Liu 06/19/06

class myarc               // added by Yang Liu 06/21/06
{
 public:
    long tail;
    long head;
    cType           cap;    
    myarc(long t, long h, cType C)
	{tail=t; head=h; cap=C;}
};


////////////////////////////////////////////////////////////////////////////////
#ifdef EXCESS_TYPE_LONG
typedef unsigned long excessType;
#else
typedef unsigned long long excessType; /* change to double if not supported */
#endif


/* arc */
typedef  struct arcSt
{
   cType           resCap;          /* residual capacity */
   struct nodeSt   *head;           /* arc head */
   struct arcSt    *rev;            /* reverse arc */
} arc;


/* node */
typedef  struct nodeSt
{
   arc             *first;           /* first outgoing arc */
   arc             *current;         /* current outgoing arc */
   excessType      excess;           /* excess at the node 
				        change to double if needed */
   long            d;                /* distance label */
   struct nodeSt   *bNext;           /* next node in bucket */
   struct nodeSt   *bPrev;           /* previous node in bucket */
} node;


/* bucket */
typedef  struct bucketSt
{
  node             *firstActive;      /* first node with positive excess */
  node             *firstInactive;    /* first node with zero excess */
} bucket;

#endif /* !TYPE_H */
