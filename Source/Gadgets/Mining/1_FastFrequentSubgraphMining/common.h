#ifndef COMMON_H
#define COMMON_H

//OS type

//induced or not induced
//#define INDUCEDS 1

//maximal or not
//#define MAXIMALPATTERN 1

//validation related macros
//#define VALIDATION 1

//running in safemode
//#define SAFEMODE	1

//requiring density
#define DENSITY 1

//keep the search results
#define SAVERESULTS   2

//search the frequent subgraphs based on spanning tree, 1 for closed tree
//2 for normal tree
//#define TREEBASEDSCAN 1

//using VF lib
//#define VFLIB 1

//use the minimal label
//#define MINIMAL_LABEL 0 

//output the parital results
#define BUFFERSIZE 100000
#define PATTERN_SIZE_LOW 4
#define PATTERN_SIZE_UP  10

//basic set up
#define EMPTY_EDGE '0'
#define NODE_TOKEN "node"
#define EDGE_TOKEN "edge"
#define LOGFILE    "sym.log"

//debug flags
//debug 1 for carnoical form related debug information
//debug 2 for candidate proposing related debug
//debug 3 for flow related debug
//debug 4 for sub carnonical form related debug
//debug 5 for initialization related debug
//debug 6 for large dataset debugging
//debug 7 for detailed candidate proposing related debug
//debug 8 for outer candidate proposing related debug
//debug 9 ?
//debug 10 for filter debug
//debug 11 for basic level maximal debug
//debug 12 for level 2 maximal pattern debug

//#define DEBUG1 1
//#define DEBUG2 1
//#define DEBUG3 1
//#define DEBUG4 1
//#define DEBUG5 1
//#define DEBUG6 1
//#define DEBUG7 1
//#define DEBUG8 1
//#define DEBUG9 1
//#define DEBUG10 1
//#define DEBUG11 1
//#define DEBUG12 1

//timer related macros
//#define TIMERFLAG 1 
#define TIMER1   0
#define TIMER2   1
#define TIMER3   2
#define TIMER4	 3
#define TIMER5   4
#define TIMER6   5
#define TIMER7   6

//registers related macros
//#define REGISTERFLAG 1
#define REGISTER1   0
#define REGISTER2   1
#define REGISTER3   2
#define REGISTER4   3
#define REGISTER5   4
#define REGISTER6   5
#define REGISTER7   6
#define REGISTER8   7
#define REGISTER9   8
#define REGISTER10  9
#define REGISTER11  10
#define REGISTER12  11
#define REGISTER13  12

//important types
#include <iostream>
#include <fstream>
#include <set>
#include <map>
#include <vector>
//#include <list>

#ifdef _WIN32
	//#include "stl_hash.h"
	//using std::hash_map;
	//#include <hash_map>
	#define ASSOCONTAINOR map
#else
	#include <ext/hash_map>
	#include <ext/hash_set>
    using __gnu_cxx::hash_map;
	using __gnu_cxx::hash_set;
	#define ASSOCONTAINOR hash_map
#endif
using namespace std;

//common types
#ifndef _WIN32
	typedef long long DLONG;
#else
	typedef _int64 DLONG;
#endif

#define  HALFDLONG 32
#define  QUADLONG  16 
#define  GNULLINDEX 0x0000ffff
#define  INULLINDEX 0x7fffffff
	
typedef vector<DLONG>  COCCS;
typedef int   BTYPE;
typedef unsigned char BYTE;

//typedef unsigned char  GNSIZE;
typedef short GNSIZE;
typedef char  GLTYPE;

//macros use Ullman's backtract
#define   UEDGE     '1'
#define   UNEDGE    '0'

typedef vector<GNSIZE>      GRAPHS1;
typedef map<DLONG, int>     TRA; 
typedef vector<int>         IIDTYPE;
typedef vector<IIDTYPE*>    FREELE;
typedef ASSOCONTAINOR<int, int>  CGOCC;
typedef ASSOCONTAINOR<int, IIDTYPE*>  ELEFRENQC;

//related macros
#define MAKEKEY(x, y, z) ( ( ( (DLONG)(x) ) << HALFDLONG ) |  ((DLONG)(y) << QUADLONG ) | ((DLONG)(z)) )
#define FIRSTKEY(x)   ( (int) ( (x)  >> HALFDLONG ) )
#define SECONDKEY(x)  ( (int) ( ((x)  << HALFDLONG  >>  HALFDLONG) >> QUADLONG ) )
#define THIRDKEY(x)   ( ((int) ((x)  << HALFDLONG  >>  HALFDLONG) ) & 0x0000ffff)

typedef int COSTTYPE; 
#define COSTLENGTH  16
#define LMASK   0x000000ff
#define MAKECOST(x, y, z)  ( ( ((COSTTYPE) (x))  << COSTLENGTH ) | ( (y) << (COSTLENGTH>>1) ) | (z)  )
#define COSTINDEX(x)  ( (x)  >> COSTLENGTH ) 
#define COSTEL(x)     ( ((x) >> (COSTLENGTH>>1) ) & LMASK )
#define COSTNL(x)     ( (x) & LMASK )

#include "myTimer.h"
#include "register.h"

inline void error(const char *p, const char *p2 = " ")
{
	cerr << p << ' ' << p2 << endl;
	exit(1);
}

inline void log(const char *p, const char *p2 = " ")
{

#ifndef _WIN32
	ofstream outfile(LOGFILE, ios::app);
#else
	ofstream outfile(LOGFILE, ios_base::app);
#endif

	if( !outfile )
		error("log file is not ready for write: ", LOGFILE);

	outfile << p << ' ' << p2 << endl;
	outfile.close();

}

#endif
