#ifndef PARSER_H
#define PARSER_H
#define MAXLINE       100	/* max line length in the input file */
#define ARC_FIELDS      3	/* No. of fields in arc line  */
#define NODE_FIELDS     2	/* No. of fields in node line  */
#define P_FIELDS        3       /* No. of fields in problem line */
#define PROBLEM_TYPE "max"      /* name of problem type */
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

#define MEMORY_ALLOCATION_FAILS {cout<<"Memory Allocation Fails"<<endl; exit(1);}

#define STORE_ARC \
{\
arc_first[tail + 1] ++;                             \
arc_first[head + 1] ++;                             \
arc_tail[pos_current]         = tail;               \
arc_tail[pos_current+1]       = head;               \
arc_current       -> head     = nodes + head;       \
arc_current       -> resCap   = cap;                \
arc_current       -> rev      = arc_current + 1;    \
( arc_current + 1 ) -> head   = nodes + tail;       \
( arc_current + 1 ) -> resCap = 0;                  \
( arc_current + 1 ) -> rev    = arc_current;        \
arc_current += 2;                                   \
pos_current += 2;                                   \
}

/* -------------- error numbers & error messages ---------------- */
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
/* --------------------------------------------------------------- */
#endif /* !PARSER_H */
