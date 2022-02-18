#ifndef PARSER_H
#define PARSER_H
#define MAXLINE       100
#define ARC_FIELDS      3
#define NODE_FIELDS     2
#define P_FIELDS        3      
#define PROBLEM_TYPE "max"     
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


    static char *err_message[] = 
	{ 
	       "more than one problem line.",
	       "wrong number of parameters in the problem line.",
	       "it is not a Max Flow problem line.",
	       "bad value of a parameter in the problem line.",
	       "can't obtain enough memory to solve this problem.",
	       "more than one line with the problem name.",
	       "can't read problem name.",
	       "problem description must be before node description.",
	       "this parser doesn't support multiply sources and sinks.",
	       "wrong number of parameters in the node line.",
	       "wrong value of parameters in the node line.",
	       " ",
	       "source and sink descriptions must be before arc descriptions.",
	       "too many arcs in the input.",
	       "wrong number of parameters in the arc line.",
	       "wrong value of parameters in the arc line.",
	       "unknown line type in the input.",
	       "reading error.",
	       "not enough arcs in the input.",
	       "source or sink doesn't have incident arcs.",
	       "can't read anything from the input file."
	};

#endif
