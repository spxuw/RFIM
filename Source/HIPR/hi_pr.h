#ifndef HIPR_H
#define HIPR_H
#include <sys/time.h>
#include <sys/resource.h>
#include <sys/stat.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include <limits.h>
#include "State.h"   
#include "types.h"                      

extern  const hType RESOLUTION;
#define J     RESOLUTION          
#define CC    (4*J)               

#define GLOB_UPDT_FREQ 0.5
#define ALPHA 6
#define BETA 12
#define WHITE 0
#define GREY 1
#define BLACK 2

#define forAllNodes(i) for ( i = nodes; i != sentinelNode; i++ )
#define forAllArcs(i,a) for (a = i->first, stopA = (i+1)->first; a != stopA; a++)

#define nNode( i ) ( (i) - nodes ) 
#define nArc( a )  ( ( a == NULL )? -1 : (a) - arcs )
#define min( a, b ) ( ( (a) < (b) ) ? a : b )

#define aAdd(l,i)\
{\
  i->bNext = l->firstActive;\
  l->firstActive = i;\
  i_dist = i->d;\
  if (i_dist < aMin)\
    aMin = i_dist;\
  if (i_dist > aMax)\
    aMax = i_dist;\
  if (dMax < aMax)\
    dMax = aMax;\
}


#define aRemove(l,i)\
{\
  l->firstActive = i->bNext;\
}

#define iAdd(l,i)\
{\
  i_next = l->firstInactive;\
  i->bNext = i_next;\
  i->bPrev = sentinelNode;\
  i_next->bPrev = i;\
  l->firstInactive = i;\
}

#define iDelete(l,i)\
{\
  i_next = i->bNext;\
  if (l->firstInactive == i) {\
    l->firstInactive = i_next;\
    i_next->bPrev = sentinelNode;\
  }\
  else {\
    i_prev = i->bPrev;\
    i_prev->bNext = i_next;\
    i_next->bPrev = i_prev;\
  }\
}

int allocDS();
void init();
void checkMax();
void globalUpdate();
void stageTwo();
int gap(bucket* emptyB);
long relabel(node* i);
void discharge(node* i);
void wave();
void stageOne();


void initilize_pointer();
float timer();

void DeleTempfiles(int D, int L, double R, int seed);
void DeleTempfiles(int N, double c, double R, int seed);
void DeleTempfiles(string file);

void Calculate_GS_M(State &W, bool uES, char FZdir); 


#endif



