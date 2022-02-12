//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// The following function is added by Yang Liu 06/12/06. Goal: get the ground state of the RFIM from the min-cut of the 
// corresponding network using the HIPR algorithm. But here we use the earlier solution as an input to remove some spins 
// (nodes) from the ground state calculation.  
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#include "hi_pr.cpp"

//uES: use earlier solution
//FZdir: frozen spin direction, for increasing field, UP spins are frozen; for decreasing field, DOWN spins are frozen.
void Calculate_GS_M(State &W, bool uES, char FZdir)
{

    initilize_pointer(); // to avoid abnormal-pointer

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

    ////////////////////////////////////////////////////////////////
    //cout<< "\n start parsing: "; // debug
    parse(W, uES, FZdir, n, m, nodes, arcs, caps, source, sink);
    //cout<< "\n parse successfully!"; // debug
    ////////////////////////////////////////////////////////////////

    cc = allocDS();
    if ( cc ) {fprintf ( stderr, "Allocation error\n");	exit ( 1 );}

    init(); 

    //cout<< "\n Stage One : "; // debug
    stageOne ( );

    //printf ("c flow:       %12.01f\n", flow);

 
#ifndef CUT_ONLY
    stageTwo ( );
#endif


#ifdef CHECK_SOLUTION

    // check if you have a flow (pseudoflow) 
    // check arc flows 
    forAllNodes(i) {
	forAllArcs(i,a) {
	    if (caps[a - arcs] > 0) // original arc 
		if ((a->resCap + a->rev->resCap != caps[a - arcs]) 
		    || (a->resCap < 0)
		    || (a->rev->resCap < 0)) {
		    printf("ERROR: bad arc flow\n");
		    exit(2);
		}
	}
    }

    // check conservation 
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
		if (caps[a - arcs] > 0) // original arc 
		    sum -= caps[a - arcs] - a->resCap;
		else
		    sum += a->resCap;
	    }

	    if (i->excess != sum) {
		printf("ERROR: conservation constraint violated\n");
		exit(2);
	    }
	}

    // check if mincut is saturated 
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
  


//////////////////////////////////////////////////////////////////////////////////////////////////////
//Get the Ground state spin configuration from the min-cut
#ifdef PRINT_CUT  
    globalUpdate();
    //printf ("\n Check nodes on the sink side or not (i.e. spin DOWN or UP )\n"); 

    // This is done by checking whether j->d >= n. if j->d >= n, then this node 
    // (index given by nNode(j)) is on the source side, in other words, according 
    // to the definition given in Hartmann's book, x[nNode(j)] = 1, and the 
    // corresponding spin value is given by  S[nNode(j)] = 2*x[nNode(j)]-1 = +1, 
    // i.e. spin UP!  Comments by YL 05/24/06. 

    if(uES) // using earlier solution /////////////////////////////////////////////////////////////////
    {
	if(FZdir==1) // this means there are many frozen UP spins
	{            // we just need to flip UP a small number of spins
	    int Nup = W.Get_Nup();                 // number of earlier UP spins
	    forAllNodes(j){
		int nj = nNode(j);
		if (j->d >= n && nj != 0 && nj != n-1)  // ATTENTION: j->d >= n (source side) means SPINUP
		{
		    int loc = W.LocVec[nj];        // the node index is given by nj-nMin, the 
		    // corresponding spin loc is then given by LocVec[nj-nMin]
		    W.FlipUP(loc);                 // flip this spin to be UP.
		    Nup++;                        
		}
	    }
	    W.CalMbyNup(Nup);
	}//end of if frozen direction=UP
	else  // this means there are many frozen DOWN spins
	{     // we just need to flip DOWN a small number of spins
	    int Ndn = W.Get_Ndn();                // number of earlier DOWN spins
	    forAllNodes(j) {
		int nj = nNode(j);
		if (j->d < n && nj != 0 && nj != n-1)  // ATTENTION: j->d < n (sink side) means SPINDN
		{
		    int loc = W.LocVec[nj];       // the node index is given by nj-nMin, the 
		    // corresponding spin loc is then given by LocVec[nj-nMin]
		    W.FlipDN(loc);                // flip this spin to be DOWN.
		    Ndn++;                        
		}
	    }
  
	    W.CalMbyNdn(Ndn);
	}//end of if frozen direction=DOWN
    }// end of if use earlier solution

    else // don't use earlier solution  //////////////////////////////////////////////////////////////////////
    {

      /*------------------------ See my debug note on 04/18/2014 ------------------------------------------------------------------------------
	double H = W.Get_Hext();

	if(H <= 0.0) // this means the GS will have M<0, i.e. Nup<Ndn
	{
	 
	    int Nup=0;  // number of initial UP spins
	    forAllNodes(j){
		int nj = nNode(j);
		if (j->d >= n && nj != 0 && nj != n-1) { // ATTENTION: j->d >= n (source side) means SPINUP
		    W.FlipUP(nj-1); // node index : nj ---> spin loc= nj-1
		    // flip the spin with location: nj-1 to be UP.
		    // Note that in this case, the initial spins must all be DOWN
		    // If initially all spins UP, then finnally all spins stay UP
		    // So pay attention to the initial value of spins. Y.L. 06/12/06
		    Nup++;                        
		}
	    }
	    W.CalMbyNup(Nup);

	    //cout << "M=" << W.Getm() << endl;//test

	}// end of if H<=0
	else // this means the GS will have M>0, i.e. Ndn<Nup
	{
	    int Ndn=0;  // number of initial DOWN spins
	    forAllNodes(j){
		int nj = nNode(j);
		if (j->d < n && nj != 0 && nj != n-1) {  // ATTENTION: j->d < n (sink side) means SPINDN
		    W.FlipDN(nj-1); // node index : nj ---> spin loc= nj-1
		    // flip the spin with location: nj-1 to be DOWN.
		    // Note that in this case, the initial spins must all be UP
		    Ndn++;                        
		}
	    }
	    W.CalMbyNdn(Ndn);

	    //cout << "M=" << W.Getm() << endl;//test
	    
	}// end of if H>0
	------------------------------------------------------------------------------------------------------*/

      
      int Nup = 0;
      forAllNodes(j){
	int nj = nNode(j);
	if(nj != 0 && nj != n-1) {
	  if (j->d >= n) { // ATTENTION: j->d >= n (source side) means SPINUP
	    W.SetUP(nj-1); // node index : nj ---> spin loc= nj-1, set the spin with location: nj-1 to be UP.
	    Nup++;                        
	  }
	  else { // ATTENTION: j->d < n (sink side) means SPINDN
	    W.SetDN(nj-1);
	  }
	}
      }
      W.CalMbyNup(Nup);
      

    }// end of if don't use earlier solution

#endif  
///////////////////////////////////////////////////////////////////////////////////////////////////////



// free internal memory , added by Yang Liu
    free(nodes);   nodes   = NULL;
    free(arcs);    arcs    = NULL;
    free(caps);    caps    = NULL;
    free(buckets); buckets = NULL;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

