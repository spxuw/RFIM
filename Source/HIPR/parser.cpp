#include "parser.h"

int parse( State  &W,            
	   bool   &uES,               
	   char   &FZdir,           
	   long   &n,              
	   long   &m,              
	   node*  &nodes,         
	   arc*   &arcs,          
	   cType* &caps,         
	   node*  &source,        
	   node*  &sink)          
{

    long    
        node_max=0,                 
	*arc_first=NULL,            
	*arc_tail=NULL,             
        s=0,                   
        t=0,                     
	head, tail, i;              

    long    
        last, arc_num, arc_new_num; 

    node 
	*head_p=NULL, 
	*ndp=NULL;

    arc 
        *arc_current=NULL,
        *arc_new=NULL,
        *arc_tmp=NULL;

    cType 
      cap;                        

    long  no_lines=0;              
    long  pos_current=0;           
    int err_no;                    

   
    if(uES)
    {
	W.Mapping_Frozenspins(FZdir, n, m, s, t);
	
    }
    else
    {
	W.Mapping(n, m, s, t);
	
    }

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


    nodes->first = arcs;  

    for(i=1; i<=node_max + 1; i++) 
    {
	arc_first[i]    += arc_first[i-1];
	(nodes+i)->first = arcs + arc_first[i];
    }

    for(i=0; i<node_max; i++) 
    {
	last = ((nodes+i+1)->first) - arcs;


	for(arc_num=arc_first[i]; arc_num<last; arc_num++)
	{ 
	    tail = arc_tail[arc_num];

	    while(tail != i)
	    { 
		arc_new_num  = arc_first[tail];
		arc_current  = arcs + arc_num;
		arc_new      = arcs + arc_new_num;
                            
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

		arc_first[tail]++;

		tail = arc_tail[arc_num];
	    }
	}
    }       



    for ( ndp = nodes; ndp <= nodes + node_max;  ndp ++ )
	ndp->first = (arc*) NULL;

    for ( arc_current = arcs + (2*m-1); arc_current >= arcs; arc_current -- )
    {
	arc_num = arc_current - arcs;
	tail = arc_tail [arc_num];
	ndp = nodes + tail;

	ndp->first = arc_current;
    }

    source = nodes + s;
    sink   = nodes + t;
  
    for ( arc_current = arcs, arc_num = 0; arc_num < 2*m;  arc_current ++, arc_num ++)
	caps [ arc_num ] = arc_current -> resCap; 

    if ( s < 0 || s > node_max ) 
    { err_no = EN20; goto error; }
  
    if ( (source) -> first == (arc*) NULL || (sink) -> first == (arc*) NULL ) 
    { err_no = EN20; goto error; } 	


    free(arc_first);    arc_first = NULL;
    free(arc_tail);     arc_tail  = NULL;


    return (0);


 error:  
    printf ( "\nline %ld of input - %s\n", no_lines, err_message[err_no] );
    exit (1);

}
