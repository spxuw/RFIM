/************************************************************************************************************************/
/*  parse (...) :
    1. Reads maximal flow problem in extended DIMACS format directly from the program (without reading from a input file).
    2. Prepares internal data representation.
       
    types: 'arc' and 'node' must be predefined
    type arc  must contain fields 'head', 'rev', 'resCap'
    type node must contain the field 'first': */
/************************************************************************************************************************/
#include "parser.h"

// all parameters are output, except the first one
int parse( State  &W,             /* address of the State W */
	   bool   &uES,               // whether using earlier solution
	   char   &FZdir,             // if uES=true, then what is the direction of the frozen spins 
	   long   &n,              /* address of the number of nodes */
	   long   &m,              /* address of the number of arcs */
	   node*  &nodes,         /* address of the array of nodes */ 
	   arc*   &arcs,           /* address of the array of arcs */
	   cType* &caps,          /* address of the array of capasities */  
	   node*  &source,        /* address of the pointer to the source */
	   node*  &sink)          /* address of the pointer to the sink */
{

    long    //n,                      /* internal number of nodes */
        node_max=0,                 /* maximal No. of nodes */
	*arc_first=NULL,            /* internal array for holding - node degree - position of the first outgoing arc */
	*arc_tail=NULL,             /* internal array: tails of the arcs */
        s=0,                   /* No. of the source */
        t=0,                     /* No. of the sink */
	head, tail, i;              /* temporary variables carrying No. of nodes */

    long    //m,                      /* internal number of arcs */
        last, arc_num, arc_new_num; /* temporary variables carrying No. of arcs */

    node //*nodes=NULL,               /* pointer to the node structure */
	*head_p=NULL, 
	*ndp=NULL;

    arc //*arcs=NULL,                 /* pointer to the arc structure */
        *arc_current=NULL,
        *arc_new=NULL,
        *arc_tmp=NULL;

    cType //*caps=NULL,    /* array of capasities */
      cap;                        /* capasity of the current arc */

    long  no_lines=0;               /* No. of current input line */
    long  pos_current=0;            /* 2*No._alines */
    int err_no;                     /* No. of detected error */




/////////////////////////////////////////////////////////////////////////////////////////////////
// The main part:  mapping the RFIM into a network, Yang Liu 05/26/06 
    // get the basic information of network from the random field realization 

    //cout << "start Mapping: with uES= " << (int)uES << endl; //debug 
   
    if(uES)
    {
	W.Mapping_Frozenspins(FZdir, n, m, s, t);
	//W.Mapping_Frozenspins(FZdir); //test
    }
    else
    {
	W.Mapping(n, m, s, t);
	//W.Mapping();  //test
    }

    //cout << "Mapping done!\n"; //debug 

    // I realize that using `new` to allocate memory, the program will have trouble when dealing with  
    // some special case, like D=3, L=17, 18,....,31. I am very glad to know this! :) Yang Liu 06/08/06
    /*// allocating memory for  'nodes', 'arcs' ,'acap' and internal arrays 
    nodes    = new node [n+2];          if(nodes==NULL)     MEMORY_ALLOCATION_FAILS;
    arcs     = new arc [2*m+1];         if(arcs==NULL)      MEMORY_ALLOCATION_FAILS;
    acap     = new unsigned long [2*m]; if(acap==NULL)      MEMORY_ALLOCATION_FAILS;
    arc_tail = new long [2*m];          if(arc_tail==NULL)  MEMORY_ALLOCATION_FAILS;
    arc_first= new long [n+2];          if(arc_first==NULL) MEMORY_ALLOCATION_FAILS;
    */ 

    /* allocating memory for  'nodes', 'arcs'  and internal arrays */
    arc_tail = (long*) calloc ( 2*m,   sizeof(long) ); 
    arc_first= (long*) calloc ( n+2, sizeof(long) );

    nodes    = (node*) calloc ( n+2, sizeof(node) );
    arcs     = (arc*)  calloc ( 2*m+1, sizeof(arc) );
    caps     = (cType*) calloc ( 2*m, sizeof(cType) );
     
    if ( caps  == NULL || nodes == NULL ||  arcs  == NULL || 
	 arc_first == NULL  || arc_tail == NULL ) 
    { 
	cout << "memory allocation fails!";
	exit(1);
    }
	
    // setting pointer to the first arc 
    arc_current = arcs;   
    node_max = t;

    list<myarc>::iterator lmi; 
    for(lmi=W.ArcList.begin(); lmi!=W.ArcList.end(); lmi++)
    {
	tail = (*lmi).tail;
	head = (*lmi).head;
	cap  = (*lmi).cap;
	STORE_ARC;
    }
/////////////////////////////////////////////////////////////////////////////////////////////////////////////




/////////////////////////////////////////////////////////////////////////////////////////////////////////////
// ordering arcs - linear time algorithm 

    nodes->first = arcs;  // first arc from the first node 

/* before below loop arc_first[i+1] is the number of arcs outgoing from i;
   after this loop arc_first[i] is the position of the first 
   outgoing from node i arcs after they would be ordered;
   this value is transformed to pointer and written to node.first[i]
*/
    for(i=1; i<=node_max + 1; i++) 
    {
	arc_first[i]    += arc_first[i-1];
	(nodes+i)->first = arcs + arc_first[i];
    }

    for(i=0; i<node_max; i++) /* scanning all the nodes  except the last*/
    {
	last = ((nodes+i+1)->first) - arcs;
	/* arcs outgoing from i must be cited from position arc_first[i] to the position 
           equal to initial value of arc_first[i+1]-1  */ 

	for(arc_num=arc_first[i]; arc_num<last; arc_num++)
	{ 
	    tail = arc_tail[arc_num];

            /* the arc No. arc_num is not in place because arc cited here
	       must go out from i; we'll put it to its place and continue this process
     	       until an arc in this position would go out from i */
	    while(tail != i)
	    { 
		arc_new_num  = arc_first[tail];
		arc_current  = arcs + arc_num;
		arc_new      = arcs + arc_new_num;
	    
		/* arc_current must be cited in the position arc_new    
		   swapping these arcs:                                 */
		head_p             = arc_new->head;
		arc_new->head      = arc_current->head;
		arc_current->head  = head_p;

		cap                 = arc_new->resCap;
		arc_new->resCap     = arc_current->resCap;
		arc_current->resCap = cap;

		if(arc_new != arc_current->rev)
		{
		    arc_tmp          = arc_new->rev;
		    arc_new->rev     = arc_current->rev;
		    arc_current->rev = arc_tmp;

		    (arc_current->rev)->rev = arc_current;
		    (arc_new->rev)->rev     = arc_new;
		}

		arc_tail[arc_num]     = arc_tail[arc_new_num];
		arc_tail[arc_new_num] = tail;

		/* we increase arc_first[tail]  */
		arc_first[tail]++;

		tail = arc_tail[arc_num];
	    }
	}
	/* all arcs outgoing from  i  are in place */
    }       

/* -----------------------  arcs are ordered  ------------------------- */
/////////////////////////////////////////////////////////////////////////////////////////////////////////////


/////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*----------- constructing lists ---------------*/

    for ( ndp = nodes; ndp <= nodes + node_max;  ndp ++ )
	ndp->first = (arc*) NULL;

    for ( arc_current = arcs + (2*m-1); arc_current >= arcs; arc_current -- )
    {
	arc_num = arc_current - arcs;
	tail = arc_tail [arc_num];
	ndp = nodes + tail;
	/* avg
	   arc_current -> next = ndp -> first;
	*/
	ndp->first = arc_current;
    }
/////////////////////////////////////////////////////////////////////////////////////////////////////////////


/////////////////////////////////////////////////////////////////////////////////////////////////////////////
/* ----------- assigning output values ------------*/
    source = nodes + s;
    sink   = nodes + t;
  
    for ( arc_current = arcs, arc_num = 0; arc_num < 2*m;  arc_current ++, arc_num ++)
	caps [ arc_num ] = arc_current -> resCap; 

    if ( s < 0 || s > node_max ) /* bad value of the source */
    { err_no = EN20; goto error; }
  
    if ( (source) -> first == (arc*) NULL || (sink) -> first == (arc*) NULL ) 
    { err_no = EN20; goto error; } 	/* no arc goes out of the source */
/////////////////////////////////////////////////////////////////////////////////////////////////////////////


/* free internal memory */
    free(arc_first);    arc_first = NULL;
    free(arc_tail);     arc_tail  = NULL;


/* Thanks God! all is done */
    return (0);


/* error found reading input */
 error:  
    printf ( "\nline %ld of input - %s\n", no_lines, err_message[err_no] );
    exit (1);

}
/* --------------------   end of parser  -------------------*/
