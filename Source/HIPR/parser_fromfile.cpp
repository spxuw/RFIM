#define MAXLINE       100	
#define ARC_FIELDS      3	
#define NODE_FIELDS     2	
#define P_FIELDS        3       
#define PROBLEM_TYPE "max"      

int parse( long* n_ad, long* m_ad, node** nodes_ad, arc** arcs_ad, 
	   unsigned long** cap_ad,  node** source_ad, node** sink_ad, long* node_min_ad )
{

    long    n,                      
        node_min=0,                 
        node_max=0,                 
	*arc_first=NULL,            
	*arc_tail=NULL,             
        source=0,                   
        sink=0,                     
        head, tail, i;

    long    m,                     
        last, arc_num, arc_new_num;

    node    *nodes=NULL,           
        *head_p,
        *ndp;

    arc     *arcs=NULL,             
        *arc_current=NULL,
        *arc_new,
        *arc_tmp;

    unsigned long    *acap=NULL,             
	cap;                    

    long    no_lines=0,             
        no_plines=0,            
        no_nslines=0,           
        no_nklines=0,           
        no_alines=0,            
        pos_current=0;         

    char    in_line[MAXLINE],      
        pr_type[3],             
        nd;                     

    int     k,                  
        err_no;                



#define EN1   0
#define EN2   1
#define EN3   2
#define EN4   3
#define EN6   4
#define EN10  5
#define EN7   6
#define EN8   7
#define EN9   8
#define EN11  9
#define EN12 10
#define EN13 11
#define EN14 12
#define EN16 13
#define EN15 14
#define EN17 15
#define EN18 16
#define EN21 17
#define EN19 18
#define EN20 19
#define EN22 20

    static char *err_message[] = 
	{ 
	    /* 0*/    "more than one problem line.",
	    /* 1*/    "wrong number of parameters in the problem line.",
	    /* 2*/    "it is not a Max Flow problem line.",
	    /* 3*/    "bad value of a parameter in the problem line.",
	    /* 4*/    "can't obtain enough memory to solve this problem.",
	    /* 5*/    "more than one line with the problem name.",
	    /* 6*/    "can't read problem name.",
	    /* 7*/    "problem description must be before node description.",
	    /* 8*/    "this parser doesn't support multiply sources and sinks.",
	    /* 9*/    "wrong number of parameters in the node line.",
	    /*10*/    "wrong value of parameters in the node line.",
	    /*11*/    " ",
	    /*12*/    "source and sink descriptions must be before arc descriptions.",
	    /*13*/    "too many arcs in the input.",
	    /*14*/    "wrong number of parameters in the arc line.",
	    /*15*/    "wrong value of parameters in the arc line.",
	    /*16*/    "unknown line type in the input.",
	    /*17*/    "reading error.",
	    /*18*/    "not enough arcs in the input.",
	    /*19*/    "source or sink doesn't have incident arcs.",
	    /*20*/    "can't read anything from the input file."
	};


    while (fgets(in_line, MAXLINE, stdin) != NULL )
    {
	no_lines ++;

	switch (in_line[0])
	{
	    case 'c':                  
	    case '\n':                 
	    case '\0':                
                break;


	    case 'p':                 

                if ( no_plines > 0 )   // more than one problem line 
		{ err_no = EN1 ; goto error; }

                no_plines = 1;
   
                if (sscanf ( in_line, "%*c %3s %ld %ld", pr_type, &n, &m ) != P_FIELDS  )
		{ err_no = EN2; goto error; } //wrong number of parameters in the problem line

                if ( strcmp ( pr_type, PROBLEM_TYPE ) )  //wrong problem type
		{ err_no = EN3; goto error; }

                if ( n <= 0  || m <= 0 )    //wrong value of no of arcs or nodes
		{ err_no = EN4; goto error; }

                nodes    = (node*) calloc ( n+2, sizeof(node) );
		arcs     = (arc*)  calloc ( 2*m+1, sizeof(arc) );
	        arc_tail = (long*) calloc ( 2*m,   sizeof(long) ); 
		arc_first= (long*) calloc ( n+2, sizeof(long) );
                acap     = (unsigned long*) calloc ( 2*m, sizeof(long) );

                if ( nodes == NULL || arcs == NULL || arc_first == NULL || arc_tail == NULL ) 
		{ err_no = EN6; goto error; }   // memory is not allocated 
		     
		arc_current = arcs; 		// setting pointer to the first arc 

                break;



	    case 'n':		         

		if ( no_plines == 0 )    // there was not problem line above 
		{ err_no = EN8; goto error; }

		k = sscanf ( in_line,"%*c %ld %c", &i, &nd );  // reading source or sink 
 
		if ( k < NODE_FIELDS )  // node line is incorrect 
		{ err_no = EN11; goto error; }

		if ( i < 0 || i > n )   // wrong value of node 
		{ err_no = EN12; goto error; }

		switch ( nd )
		{
		    case 's': 
			
			if ( no_nslines != 0)   // more than one source line 
			{ err_no = EN9; goto error; }

			no_nslines = 1;
			source = i;
			break;

		    case 't':  

			if ( no_nklines != 0)   // more than one sink line 
			{ err_no = EN9; goto error; }

			no_nklines = 1;
			sink = i;
			break;

		    default:   
			err_no = EN12; goto error; 
			break;
		}

                node_max = 0;
                node_min = n;
		break;



	    case 'a':                    
		if ( no_nslines == 0 || no_nklines == 0 ) // there was not source and sink description above 
		{ err_no = EN14; goto error; }

		if ( no_alines >= m )	// too many arcs on input
		{ err_no = EN16; goto error; }
		
		if (sscanf ( in_line,"%*c %ld %ld %ld", &tail, &head, &cap ) != ARC_FIELDS ) 
		{ err_no = EN15; goto error; }   // arc description is not correct 

		if ( tail < 0  ||  tail > n  || head < 0  ||  head > n  )
		{ err_no = EN17; goto error; }   // wrong value of nodes 

		arc_first[tail + 1] ++; 
		arc_first[head + 1] ++;

		arc_tail[pos_current]        = tail;
		arc_tail[pos_current+1]      = head;
		arc_current       -> head    = nodes + head;
		arc_current       -> resCap    = cap;
		arc_current       -> rev  = arc_current + 1;
		( arc_current + 1 ) -> head    = nodes + tail;
		( arc_current + 1 ) -> resCap    = 0;
		( arc_current + 1 ) -> rev  = arc_current;

                if ( head < node_min ) node_min = head;
                if ( tail < node_min ) node_min = tail;
                if ( head > node_max ) node_max = head;
                if ( tail > node_max ) node_max = tail;

		no_alines   ++;
		arc_current += 2;
		pos_current += 2;
		break;


	    default:  
		err_no = EN18; goto error;
		break;

	} 
    }     
/////////////////////////////////////////////////////////////////////////////////////////////////////////////



    if ( feof (stdin) == 0 ) 
    { err_no=EN21; goto error; } 

    if ( no_lines == 0 )     
    { err_no = EN22; goto error; } 

    if ( no_alines < m )     
    { err_no = EN19; goto error; } 




    ( nodes + node_min ) -> first = arcs;


    for ( i = node_min + 1; i <= node_max + 1; i ++ ) 
    {
	arc_first[i]          += arc_first[i-1];
	( nodes + i ) -> first = arcs + arc_first[i];
    }


    for ( i = node_min; i < node_max; i ++ ) 
    {
	last = ( ( nodes + i + 1 ) -> first ) - arcs;


	for ( arc_num = arc_first[i]; arc_num < last; arc_num ++ )
	{ 
	    tail = arc_tail[arc_num];

	    while ( tail != i )
	    { 
		arc_new_num  = arc_first[tail];
		arc_current  = arcs + arc_num;
		arc_new      = arcs + arc_new_num;
                               */
		head_p               = arc_new -> head;
		arc_new -> head      = arc_current -> head;
		arc_current -> head  = head_p;

		cap                 = arc_new -> resCap;
		arc_new -> resCap     = arc_current -> resCap;
		arc_current -> resCap = cap;

		if ( arc_new != arc_current -> rev )
		{
		    arc_tmp                = arc_new -> rev;
		    arc_new  -> rev     = arc_current -> rev;
		    arc_current -> rev  = arc_tmp;

		    ( arc_current -> rev ) -> rev = arc_current;
		    ( arc_new     -> rev ) -> rev = arc_new;
		}

		arc_tail[arc_num] = arc_tail[arc_new_num];
		arc_tail[arc_new_num] = tail;

		arc_first[tail] ++ ;

		tail = arc_tail[arc_num];
	    }
	}
    }       


    for ( ndp = nodes + node_min; ndp <= nodes + node_max;  ndp ++ )
	ndp -> first = (arc*) NULL;

    for ( arc_current = arcs + (2*m-1); arc_current >= arcs; arc_current -- )
    {
	arc_num = arc_current - arcs;
	tail = arc_tail [arc_num];
	ndp = nodes + tail;

	ndp -> first = arc_current;
    }

    *m_ad = m;
    *n_ad = node_max - node_min + 1;
    *source_ad = nodes + source;
    *sink_ad   = nodes + sink;
    *node_min_ad = node_min;
    *nodes_ad = nodes + node_min;
    *arcs_ad = arcs;
    *cap_ad = acap;

    for ( arc_current = arcs, arc_num = 0; arc_num < 2*m;  arc_current ++, arc_num ++)
	acap [ arc_num ] = arc_current -> resCap; 

    if ( source < node_min || source > node_max ) 
    { err_no = EN20; goto error; }
  
    if ( (*source_ad) -> first == (arc*) NULL || (*sink_ad  ) -> first == (arc*) NULL ) 
    { err_no = EN20; goto error; } 
/////////////////////////////////////////////////////////////////////////////////////////////////////////////


    free ( arc_first ); 
    free ( arc_tail );

    return (0);


 error:  
    printf ( "\nline %ld of input - %s\n", no_lines, err_message[err_no] );
    exit (1);

}
