#include "hi_pr.cpp"



void Calculate_GS_M(State &W, bool uES, char FZdir)
{

    initilize_pointer(); 

#if (defined(PRINT_FLOW) || defined(CHECK_SOLUTION))
    node *i = NULL;
    arc  *a = NULL;
#endif

#ifdef PRINT_FLOW
    long ni, na;
#endif

#ifdef PRINT_CUT
    node *j = NULL;
#endif

    int  cc;

#ifdef CHECK_SOLUTION
    excessType sum;
    bucket *l = NULL;
#endif

    globUpdtFreq = GLOB_UPDT_FREQ;

    
    
    parse(W, uES, FZdir, n, m, nodes, arcs, caps, source, sink);
    
    

    cc = allocDS();
    if ( cc ) {fprintf ( stderr, "Allocation error\n");	exit ( 1 );}

    init(); 

    
    stageOne ( );

    

 
#ifndef CUT_ONLY
    stageTwo ( );
#endif


#ifdef CHECK_SOLUTION

    
    
    forAllNodes(i) {
	forAllArcs(i,a) {
	    if (caps[a - arcs] > 0) 
		if ((a->resCap + a->rev->resCap != caps[a - arcs]) 
		    || (a->resCap < 0)
		    || (a->rev->resCap < 0)) {
		    printf("ERROR: bad arc flow\n");
		    exit(2);
		}
	}
    }

    
    forAllNodes(i)
	if ((i != source) && (i != sink)) {
#ifdef CUT_ONLY
	    if (i->excess < 0) {
		printf("ERROR: nonzero node excess\n");
		exit(2);
	    }
#else
	    if (i->excess != 0) {
		printf("ERROR: nonzero node excess\n");
		exit(2);
	    }
#endif

	    sum = 0;
	    forAllArcs(i,a) {
		if (caps[a - arcs] > 0) 
		    sum -= caps[a - arcs] - a->resCap;
		else
		    sum += a->resCap;
	    }

	    if (i->excess != sum) {
		printf("ERROR: conservation constraint violated\n");
		exit(2);
	    }
	}

    
    aMax = dMax = 0;
    for (l = buckets; l < buckets + n; l++) {
	l->firstActive = sentinelNode;
	l->firstInactive = sentinelNode;
    }
    globalUpdate();
    if (source->d < n) {
	printf("ERROR: the solution is not optimal\n");
	exit(2);
    }

    printf("\n Solution checks (feasible and optimal)\n");
#endif


#ifdef PRINT_STAT
    printf ("c pushes:      %10ld\n", pushCnt);
    printf ("c relabels:    %10ld\n", relabelCnt);
    printf ("c updates:     %10ld\n", updateCnt);
    printf ("c gaps:        %10ld\n", gapCnt);
    printf ("c gap nodes:   %10ld\n", gNodeCnt);
    printf ("c\n");
#endif


#ifdef PRINT_FLOW
    printf ("c flow values\n");
    forAllNodes(i) {
	ni = nNode(i);
	forAllArcs(i,a) {
	    na = nArc(a);
	    if ( caps[na] > 0 )
		printf ( "f %7ld %7ld %12ld\n",
			 ni, nNode( a -> head ), caps[na] - ( a -> resCap )
		    );
	}
    }
    printf("c\n");
#endif
  




#ifdef PRINT_CUT  
    globalUpdate();
    

    
    
    
    
    

    if(uES) 
    {
	if(FZdir==1) 
	{            
	    int Nup = W.Get_Nup();                 
	    forAllNodes(j){
		int nj = nNode(j);
		if (j->d >= n && nj != 0 && nj != n-1)  
		{
		    int loc = W.LocVec[nj];        
		    
		    W.FlipUP(loc);                 
		    Nup++;                        
		}
	    }
	    W.CalMbyNup(Nup);
	}
	else  
	{     
	    int Ndn = W.Get_Ndn();                
	    forAllNodes(j) {
		int nj = nNode(j);
		if (j->d < n && nj != 0 && nj != n-1)  
		{
		    int loc = W.LocVec[nj];       
		    
		    W.FlipDN(loc);                
		    Ndn++;                        
		}
	    }
  
	    W.CalMbyNdn(Ndn);
	}
    }

    else 
    {

      
      int Nup = 0;
      forAllNodes(j){
	int nj = nNode(j);
	if(nj != 0 && nj != n-1) {
	  if (j->d >= n) { 
	    W.SetUP(nj-1); 
	    Nup++;                        
	  }
	  else { 
	    W.SetDN(nj-1);
	  }
	}
      }
      W.CalMbyNup(Nup);
      

    }

#endif  





    free(nodes);   nodes   = NULL;
    free(arcs);    arcs    = NULL;
    free(caps);    caps    = NULL;
    free(buckets); buckets = NULL;
}


