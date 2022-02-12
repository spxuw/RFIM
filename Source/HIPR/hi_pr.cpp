//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// This file is added by Yang Liu 06/02/06
// Goal: get the ground state of the RFIM from the min-cut of the corresponding network 
// using the HIPR algorithm 
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#include "hi_pr.h"                      // subroutines used in hi_pr.c, Yang Liu 05/24/06 
#include "parser.cpp"                   // parser 

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/* global variables */
long   n;                    /* number of nodes */
long   m;                    /* number of arcs */
long   nm;                   /* n + ALPHA * m */
long   dMax;                 /* maximum label */
long   aMax;                 /* maximum actie node label */
long   aMin;                 /* minimum active node label */
double flow;                 /* flow value */
long pushCnt   = 0;          /* number of pushes */
long relabelCnt= 0;          /* number of relabels */
long updateCnt = 0;          /* number of updates */
long gapCnt    = 0;          /* number of gaps */
long gNodeCnt  = 0;          /* number of nodes after gap */  
long workSinceUpdate=0;      /* the number of arc scans since last update */
float globUpdtFreq;          /* global update frequency */
long i_dist;

node   *nodes;        /* array of nodes */
arc    *arcs;        /* array of arcs */
cType  *caps;        /* array of capacities */
bucket *buckets;        /* array of buckets */

node   *source;        /* source node pointer */
node   *sink;        /* sink node pointer */
node   *sentinelNode;   /* end of the node list marker */
arc    *stopA;        /* used in forAllArcs */
node   *i_next;
node   *i_prev;
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void initilize_pointer() // added by Yang Liu 06/15/06
{
    nodes        = NULL;        /* array of nodes */
    arcs         = NULL;        /* array of arcs */
    caps         = NULL;        /* array of capacities */
    buckets      = NULL;        /* array of buckets */

    source       = NULL;        /* source node pointer */
    sink         = NULL;        /* sink node pointer */
    sentinelNode = NULL;   /* end of the node list marker */
    stopA        = NULL;        /* used in forAllArcs */
    i_next       = NULL;
    i_prev       = NULL;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/* allocate datastructures, initialize related variables */
int allocDS()
{

    nm = ALPHA * n + m;
    /*
      queue = (node**) calloc ( n, sizeof (node*) );
      if ( queue == NULL ) return ( 1 );
      qLast = queue + n - 1;
      qInit();
    */
    buckets = (bucket*) calloc (n+2, sizeof(bucket) );
    if ( buckets == NULL ) return (1);

    sentinelNode = nodes + n;
    sentinelNode->first = arcs + 2*m;

    return  (0);

} /* end of allocate */
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////





//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void init()
{
    int overflowDetected;

    arc    *a = NULL;
    node   *i = NULL;        /* current node */
    bucket *l = NULL;


#ifdef EXCESS_TYPE_LONG
    double testExcess;
#endif

#ifndef OLD_INIT
    //unsigned long delta;
    unsigned long long delta; // Yang Liu 06/19/06
#endif

    // initialize excesses

    forAllNodes(i) 
	{
	    i->excess = 0;
	    i->current = i->first;
	    forAllArcs(i, a)
		a->resCap = caps[a-arcs];
	}


    for (l = buckets; l <= buckets + n-1; l++) {
	l->firstActive   = sentinelNode;
	l->firstInactive  = sentinelNode;
    }
    
    overflowDetected = 0;



#ifdef EXCESS_TYPE_LONG
    testExcess = 0;
    forAllArcs(source,a) {
	if (a->head != source) {
	    testExcess += a->resCap;
	}
    }
    if (testExcess > LONG_MAX) {
	printf("c WARNING: excess overflow. See README for details.\nc\n");
	overflowDetected = 1;
    }
#endif


#ifdef OLD_INIT
    source -> excess = LONG_MAX;
#else
    if (overflowDetected) {
	source -> excess = LONG_MAX;
    }
    else {
	source->excess = 0;
	forAllArcs(source,a) {
	    if (a->head != source) {
		pushCnt ++;
		delta = a -> resCap;
		a -> resCap -= delta;
		(a -> rev) -> resCap += delta;
		a->head->excess += delta;
	    }
	}
    }

    /*  setup labels and buckets */
    l = buckets + 1;
    
    aMax = 0;
    aMin = n;
    
    forAllNodes(i) {
	if (i == sink) {
	    i->d = 0;
	    iAdd(buckets,i);
	    continue;
	}
	if ((i == source) && (!overflowDetected)) {
	    i->d = n;
	}
	else
	    i->d = 1;
	if (i->excess > 0) {
	    /* put into active list */
	    aAdd(l,i);
	}
	else { /* i -> excess == 0 */
	    /* put into inactive list */
	    if (i->d < n)
		iAdd(l,i);
	}
    }
    dMax = 1;
#endif

    //  dMax = n-1;
    //  flow = 0.0;

} /* end of init */
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////





//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void checkMax()
{
    bucket *l = NULL;

    for (l = buckets + dMax + 1; l < buckets + n; l++) {
	assert(l->firstActive == sentinelNode);
	assert(l->firstInactive == sentinelNode);
    }
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////






//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/* global update via backward breadth first search from the sink */
void globalUpdate()
{

    node   *i = NULL;
    node   *j = NULL;        /* node pointers */
    arc    *a = NULL;        /* current arc pointers  */
    bucket *l = NULL;
    bucket *jL= NULL;        /* bucket */

    long curDist, jD;
    long state;

    updateCnt ++;

    /* initialization */

    forAllNodes(i)
	i -> d = n;
    sink -> d = 0;

    for (l = buckets; l <= buckets + dMax; l++) {
	l -> firstActive   = sentinelNode;
	l -> firstInactive  = sentinelNode;
    }

    dMax = aMax = 0;
    aMin = n;

    /* breadth first search */

    // add sink to bucket zero

    iAdd(buckets, sink);
    for (curDist = 0; 1; curDist++) {

	state = 0;
	l = buckets + curDist;
	jD = curDist + 1;
	jL = l + 1;
	/*
	  jL -> firstActive   = sentinelNode;
	  jL -> firstInactive  = sentinelNode;
	*/

	if ((l->firstActive == sentinelNode) && 
	    (l->firstInactive == sentinelNode))
	    break;

	while (1) {

	    switch (state) {
		case 0: 
		    i = l->firstInactive;
		    state = 1;
		    break;
		case 1:
		    i = i->bNext;
		    break;
		case 2:
		    i = l->firstActive;
		    state = 3;
		    break;
		case 3:
		    i = i->bNext;
		    break;
		default: 
		    assert(0);
		    break;
	    }
      
	    if (i == sentinelNode) {
		if (state == 1) {
		    state = 2;
		    continue;
		}
		else {
		    assert(state == 3);
		    break;
		}
	    }

	    /* scanning arcs incident to node i */
	    forAllArcs(i,a) {
		if (a->rev->resCap > 0 ) {
		    j = a->head;
		    if (j->d == n) {
			j->d = jD;
			j->current = j->first;
			if (jD > dMax) dMax = jD;
	    
			if (j->excess > 0) {
			    /* put into active list */
			    aAdd(jL,j);
			}
			else {
			    /* put into inactive list */
			    iAdd(jL,j);
			}
		    }
		}
	    } /* node i is scanned */ 
	}
    }

} /* end of global update */
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////






//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/* second stage -- preflow to flow */
void stageTwo()
/*
  do dsf in the reverse flow graph from nodes with excess
  cancel cycles if found
  return excess flow in topological order
*/

/*
  i->d is used for dfs labels 
  i->bNext is used for topological order list
  buckets[i-nodes]->firstActive is used for DSF tree
*/

{
    node *i, *j, *tos, *bos, *restart, *r;
    arc *a;
    cType delta;

    /* deal with self-loops */
    forAllNodes(i) {
	forAllArcs(i,a)
	    if ( a -> head == i ) {
		a -> resCap = caps[a - arcs];
	    }
    }

    /* initialize */
    tos = bos = NULL;
    forAllNodes(i) {
	i -> d = WHITE;
	//    buckets[i-nodes].firstActive = NULL;
	buckets[i-nodes].firstActive = sentinelNode;
	i -> current = i -> first;
    }

    /* eliminate flow cycles, topologicaly order vertices */
    forAllNodes(i)
	if (( i -> d == WHITE ) && ( i -> excess > 0 ) &&
	    ( i != source ) && ( i != sink )) {
	    r = i;
	    r -> d = GREY;
	    do {
		for ( ; i->current != (i+1)->first; i->current++) {
		    a = i -> current;
		    if (( caps[a - arcs] == 0 ) && ( a -> resCap > 0 )) { 
			j = a -> head;
			if ( j -> d == WHITE ) {
			    /* start scanning j */
			    j -> d = GREY;
			    buckets[j-nodes].firstActive = i;
			    i = j;
			    break;
			}
			else
			    if ( j -> d == GREY ) {
				/* find minimum flow on the cycle */
				delta = a -> resCap;
				while ( 1 ) {
				    delta = min ( delta, j -> current -> resCap );
				    if ( j == i )
					break;
				    else
					j = j -> current -> head;
				}

				/* remove delta flow units */
				j = i;
				while ( 1 ) {
				    a = j -> current;
				    a -> resCap -= delta;
				    a -> rev -> resCap += delta;
				    j = a -> head;
				    if ( j == i )
					break;
				}
	  
				/* backup DFS to the first saturated arc */
				restart = i;
				for ( j = i -> current -> head; j != i; j = a -> head ) {
				    a = j -> current;
				    if (( j -> d == WHITE ) || ( a -> resCap == 0 )) {
					j -> current -> head -> d = WHITE;
					if ( j -> d != WHITE )
					    restart = j;
				    }
				}
	  
				if ( restart != i ) {
				    i = restart;
				    i->current++;
				    break;
				}
			    }
		    }
		}

		if (i->current == (i+1)->first) {
		    /* scan of i complete */
		    i -> d = BLACK;
		    if ( i != source ) {
			if ( bos == NULL ) {
			    bos = i;
			    tos = i;
			}
			else {
			    i -> bNext = tos;
			    tos = i;
			}
		    }

		    if ( i != r ) {
			i = buckets[i-nodes].firstActive;
			i->current++;
		    }
		    else
			break;
		}
	    } while ( 1 );
	}


    /* return excesses */
    /* note that sink is not on the stack */
    if ( bos != NULL ) {
	for ( i = tos; i != bos; i = i -> bNext ) {
	    a = i -> first;
	    while ( i -> excess > 0 ) {
		if (( caps[a - arcs] == 0 ) && ( a -> resCap > 0 )) {
		    if (a->resCap < i->excess)
			delta = a->resCap;
		    else
			delta = i->excess;
		    a -> resCap -= delta;
		    a -> rev -> resCap += delta;
		    i -> excess -= delta;
		    a -> head -> excess += delta;
		}
		a++;
	    }
	}
	/* now do the bottom */
	i = bos;
	a = i -> first;
	while ( i -> excess > 0 ) {
	    if (( caps[a - arcs] == 0 ) && ( a -> resCap > 0 )) {
		if (a->resCap < i->excess)
		    delta = a->resCap;
		else
		    delta = i->excess;
		a -> resCap -= delta;
		a -> rev -> resCap += delta;
		i -> excess -= delta;
		a -> head -> excess += delta;
	    }
	    a++;
	}
    }
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////





//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/* gap relabeling */
int gap(bucket* emptyB)
{

    bucket *l;
    node  *i; 
    long  r;           /* index of the bucket before l  */
    int   cc;          /* cc = 1 if no nodes with positive excess before
			  the gap */

    gapCnt ++;
    r = ( emptyB - buckets ) - 1;

    /* set labels of nodes beyond the gap to "infinity" */
    for ( l = emptyB + 1; l <= buckets + dMax; l ++ ) {
	/* this does nothing for high level selection 
	   for (i = l -> firstActive; i != sentinelNode; i = i -> bNext) {
	   i -> d = n;
	   gNodeCnt++;
	   }
	   l -> firstActive = sentinelNode;
	*/

	for ( i = l -> firstInactive; i != sentinelNode; i = i -> bNext ) {
	    i -> d = n;
	    gNodeCnt ++;
	}

	l -> firstInactive = sentinelNode;
    }

    cc = ( aMin > r ) ? 1 : 0;

    dMax = r;
    aMax = r;

    return ( cc );

}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////







//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*--- relabelling node i */
long relabel (node* i)  //node *i;   /* node to relabel */
{

    node  *j;
    long  minD;     /* minimum d of a node reachable from i */
    arc   *minA;    /* an arc which leads to the node with minimal d */
    arc   *a;

    assert(i->excess > 0);

    relabelCnt++;
    workSinceUpdate += BETA;

    i->d = minD = n;
    minA = NULL;

    /* find the minimum */
    forAllArcs(i,a) {
	workSinceUpdate++;
	if (a -> resCap > 0) {
	    j = a -> head;
	    if (j->d < minD) {
		minD = j->d;
		minA = a;
	    }
	}
    }

    minD++;
      
    if (minD < n) {

	i->d = minD;
	i->current = minA;

	if (dMax < minD) dMax = minD;

    } /* end of minD < n */
      
    return ( minD );

} /* end of relabel */
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////





//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/* discharge: push flow out of i until i becomes inactive */
void discharge (node* i) //node  *i;
{

    node  *j;                 /* sucsessor of i */
    long  jD;                 /* d of the next bucket */
    bucket *lj;               /* j's bucket */
    bucket *l;                /* i's bucket */
    arc   *a;                 /* current arc (i,j) */
    cType  delta;
    arc *stopA;

    assert(i->excess > 0);
    assert(i != sink);
    do {

	jD = i->d - 1;
	l = buckets + i->d;

	/* scanning arcs outgoing from  i  */
	for (a = i->current, stopA = (i+1)->first; a != stopA; a++) {
	    if (a -> resCap > 0) {
		j = a -> head;

		if (j->d == jD) {
		    pushCnt ++;
		    if (a->resCap < i->excess)
			delta = a->resCap;
		    else
			delta = i->excess;
		    a->resCap -= delta;
		    a->rev->resCap += delta;

		    if (j != sink) {

			lj = buckets + jD;

			if (j->excess == 0) {
			    /* remove j from inactive list */
			    iDelete(lj,j);
			    /* add j to active list */
			    aAdd(lj,j);
			}
		    }

		    j -> excess += delta;
		    i -> excess -= delta;
	  
		    if (i->excess == 0) break;

		} /* j belongs to the next bucket */
	    } /* a  is not saturated */
	} /* end of scanning arcs from  i */

	if (a == stopA) {
	    /* i must be relabeled */
	    relabel (i);

	    if (i->d == n) break;
	    if ((l -> firstActive == sentinelNode) && 
		(l -> firstInactive == sentinelNode)
		)
		gap (l);

	    if (i->d == n) break;
	}
	else {
	    /* i no longer active */
	    i->current = a;
	    /* put i on inactive list */
	    iAdd(l,i);
	    break;
	}
    } while (1);
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////






//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// go from higher to lower buckets, push flow
void wave() 
{
    node   *i;
    bucket  *l;

    for (l = buckets + aMax; l > buckets; l--) {
	for (i = l->firstActive; i != sentinelNode; i = l->firstActive) {
	    aRemove(l,i);

	    assert(i->excess > 0);
	    discharge (i);

	}
    }
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////








//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/* first stage  -- maximum preflow*/
void stageOne()
{
    node   *i;
    bucket  *l;             /* current bucket */


#if defined(INIT_UPDATE) || defined(OLD_INIT) || defined(WAVE_INIT)
    globalUpdate ();
#endif

    workSinceUpdate = 0;

#ifdef WAVE_INIT
    wave();
#endif  

    /* main loop */
    while ( aMax >= aMin ) {
	l = buckets + aMax;
	i = l->firstActive;

	if (i == sentinelNode)
	    aMax--;
	else {
	    aRemove(l,i);

	    assert(i->excess > 0);
	    discharge (i);

	    if (aMax < aMin)
		break;

	    /* is it time for global update? */
	    if (workSinceUpdate * globUpdtFreq > nm) {
		globalUpdate ();
		workSinceUpdate = 0;
	    }

	}
    
    } /* end of the main loop */
    
    flow = sink -> excess;

} 
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////






//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// added by Yang Liu 06/26/06
void DeleTempfiles(int D, int L, double R, int seed)
{
    struct stat s;

    char jobtmp[256]; 
    sprintf(jobtmp,"./data/GS-D%d-L%d-R%.3lf-seed%d",D,L,R,seed);
    char cmd[256];
    sprintf(cmd,"rm -rf %s",jobtmp);
    if (-1!=stat(jobtmp,&s))
    {
	int  err = system(cmd);
	if (err == -1) 
	    cout << "Failed to delete dir" << jobtmp << endl;
    }
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// added by Yang Liu 04/17/2014
void DeleTempfiles(int N, double c, double R, int seed)
{
    struct stat s;

    char jobtmp[256]; 
    sprintf(jobtmp,"./data/GS-N%d-c%e-R%.3lf-seed%d",N,c,R,seed);
    char cmd[256];
    sprintf(cmd,"rm -rf %s",jobtmp);
    if (-1!=stat(jobtmp,&s))
    {
	int  err = system(cmd);
	if (err == -1) 
	    cout << "Failed to delete dir" << jobtmp << endl;
    }
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// added by Yang Liu 04/18/2014
void DeleTempfiles(string file)
{
    struct stat s;

    char jobtmp[256]; 
    sprintf(jobtmp,"./data/GS-%s",file.c_str());
    char cmd[256];
    sprintf(cmd,"rm -rf %s",jobtmp);
    if (-1!=stat(jobtmp,&s))
    {
	int  err = system(cmd);
	if (err == -1) 
	    cout << "Failed to delete dir" << jobtmp << endl;
    }
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
float timer ()
{
  struct rusage r;

  getrusage(0, &r);
  return (float)(r.ru_utime.tv_sec+r.ru_utime.tv_usec/(float)1000000);
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

