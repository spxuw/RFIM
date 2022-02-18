#include "hi_pr.h"                      
#include "parser.cpp"                   


long   n;                    
long   m;                    
long   nm;                   
long   dMax;                 
long   aMax;                 
long   aMin;                 
double flow;                 
long pushCnt   = 0;          
long relabelCnt= 0;          
long updateCnt = 0;          
long gapCnt    = 0;          
long gNodeCnt  = 0;           
long workSinceUpdate=0;      
float globUpdtFreq;          
long i_dist;

node   *nodes;        
arc    *arcs;        
cType  *caps;       
bucket *buckets;        

node   *source;        
node   *sink;       
node   *sentinelNode;   
arc    *stopA;        
node   *i_next;
node   *i_prev;





void initilize_pointer() 
{
    nodes        = NULL;        
    arcs         = NULL;        
    caps         = NULL;        
    buckets      = NULL;        

    source       = NULL;        
    sink         = NULL;        
    sentinelNode = NULL;   
    stopA        = NULL;        
    i_next       = NULL;
    i_prev       = NULL;
}





int allocDS()
{

    nm = ALPHA * n + m;

    buckets = (bucket*) calloc (n+2, sizeof(bucket) );
    if ( buckets == NULL ) return (1);

    sentinelNode = nodes + n;
    sentinelNode->first = arcs + 2*m;

    return  (0);

} 







void init()
{
    int overflowDetected;

    arc    *a = NULL;
    node   *i = NULL;        
    bucket *l = NULL;


#ifdef EXCESS_TYPE_LONG
    double testExcess;
#endif

#ifndef OLD_INIT
    
    unsigned long long delta; 
#endif

    

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
	    aAdd(l,i);
	}
	else { 
	    if (i->d < n)
		iAdd(l,i);
	}
    }
    dMax = 1;
#endif

    
    

} 





void checkMax()
{
    bucket *l = NULL;

    for (l = buckets + dMax + 1; l < buckets + n; l++) {
	assert(l->firstActive == sentinelNode);
	assert(l->firstInactive == sentinelNode);
    }
}








void globalUpdate()
{

    node   *i = NULL;
    node   *j = NULL;       
    arc    *a = NULL;        
    bucket *l = NULL;
    bucket *jL= NULL;       

    long curDist, jD;
    long state;

    updateCnt ++;


    forAllNodes(i)
	i -> d = n;
    sink -> d = 0;

    for (l = buckets; l <= buckets + dMax; l++) {
	l -> firstActive   = sentinelNode;
	l -> firstInactive  = sentinelNode;
    }

    dMax = aMax = 0;
    aMin = n;


    

    iAdd(buckets, sink);
    for (curDist = 0; 1; curDist++) {

	state = 0;
	l = buckets + curDist;
	jD = curDist + 1;
	jL = l + 1;


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

	    forAllArcs(i,a) {
		if (a->rev->resCap > 0 ) {
		    j = a->head;
		    if (j->d == n) {
			j->d = jD;
			j->current = j->first;
			if (jD > dMax) dMax = jD;
	    
			if (j->excess > 0) {
			    aAdd(jL,j);
			}
			else {
			    iAdd(jL,j);
			}
		    }
		}
	    } 
	}
    }

} 





void stageTwo()


{
    node *i, *j, *tos, *bos, *restart, *r;
    arc *a;
    cType delta;

    forAllNodes(i) {
	forAllArcs(i,a)
	    if ( a -> head == i ) {
		a -> resCap = caps[a - arcs];
	    }
    }

    tos = bos = NULL;
    forAllNodes(i) {
	i -> d = WHITE;
	
	buckets[i-nodes].firstActive = sentinelNode;
	i -> current = i -> first;
    }

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
			    j -> d = GREY;
			    buckets[j-nodes].firstActive = i;
			    i = j;
			    break;
			}
			else
			    if ( j -> d == GREY ) {
				delta = a -> resCap;
				while ( 1 ) {
				    delta = min ( delta, j -> current -> resCap );
				    if ( j == i )
					break;
				    else
					j = j -> current -> head;
				}

				j = i;
				while ( 1 ) {
				    a = j -> current;
				    a -> resCap -= delta;
				    a -> rev -> resCap += delta;
				    j = a -> head;
				    if ( j == i )
					break;
				}
	  
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







int gap(bucket* emptyB)
{

    bucket *l;
    node  *i; 
    long  r;           
    int   cc;          

    gapCnt ++;
    r = ( emptyB - buckets ) - 1;

    for ( l = emptyB + 1; l <= buckets + dMax; l ++ ) {


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









long relabel (node* i)  
{

    node  *j;
    long  minD;     
    arc   *minA;    
    arc   *a;

    assert(i->excess > 0);

    relabelCnt++;
    workSinceUpdate += BETA;

    i->d = minD = n;
    minA = NULL;

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

    } 
      
    return ( minD );

} 




void discharge (node* i) 
{

    node  *j;                
    long  jD;                
    bucket *lj;              
    bucket *l;               
    arc   *a;                
    cType  delta;
    arc *stopA;

    assert(i->excess > 0);
    assert(i != sink);
    do {

	jD = i->d - 1;
	l = buckets + i->d;

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
			    iDelete(lj,j);
			    aAdd(lj,j);
			}
		    }

		    j -> excess += delta;
		    i -> excess -= delta;
	  
		    if (i->excess == 0) break;

		} 
	    } 
	} 

	if (a == stopA) {
	    relabel (i);

	    if (i->d == n) break;
	    if ((l -> firstActive == sentinelNode) && 
		(l -> firstInactive == sentinelNode)
		)
		gap (l);

	    if (i->d == n) break;
	}
	else {
	    i->current = a;
	    iAdd(l,i);
	    break;
	}
    } while (1);
}









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










void stageOne()
{
    node   *i;
    bucket  *l;             


#if defined(INIT_UPDATE) || defined(OLD_INIT) || defined(WAVE_INIT)
    globalUpdate ();
#endif

    workSinceUpdate = 0;

#ifdef WAVE_INIT
    wave();
#endif  

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

	    if (workSinceUpdate * globUpdtFreq > nm) {
		globalUpdate ();
		workSinceUpdate = 0;
	    }

	}
    
    } 
    
    flow = sink -> excess;

} 









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




float timer ()
{
  struct rusage r;

  getrusage(0, &r);
  return (float)(r.ru_utime.tv_sec+r.ru_utime.tv_usec/(float)1000000);
}


