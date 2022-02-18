#ifndef TYPES_H
#define TYPES_H
#include <vector>
#include <iostream>
using namespace std;


typedef vector<double>  Vec;         
typedef long long      hType;        
typedef unsigned long long cType;    

class myarc               
{
 public:
    long tail;
    long head;
    cType           cap;    
    myarc(long t, long h, cType C)
	{tail=t; head=h; cap=C;}
};



#ifdef EXCESS_TYPE_LONG
typedef unsigned long excessType;
#else
typedef unsigned long long excessType;
#endif



typedef  struct arcSt
{
   cType           resCap;         
   struct nodeSt   *head;          
   struct arcSt    *rev;           
} arc;



typedef  struct nodeSt
{
   arc             *first;          
   arc             *current;        
   excessType      excess;          
   long            d;               
   struct nodeSt   *bNext;          
   struct nodeSt   *bPrev;          
} node;



typedef  struct bucketSt
{
  node             *firstActive;     
  node             *firstInactive;   
} bucket;

#endif
