#ifndef COST1
#define COST1

#include <set>
#include <map>
#include <vector>

class cost1
{
private:
	int v;
public:

	cost1 () { v = 0; };

	cost1 (int index, char elabel, char vlabel){
		v = index << 16;
		v |= ((int) elabel) << 8;
		v |= (int) vlabel;
	}

	cost1 (int u) { v = u ; }

	void setCost(int index, char elabel, char vlabel){
		v = index << 16;
		v |= ((int) elabel) << 8;
		v |= (int) vlabel;
	}


	inline void setCost(int u){ v = u; };

	inline int intValue(){ return v; }

	inline int index()   { return v >> 16 ; }
	inline char elabel() { return (char) ( (v >> 8) & 0x000000FF );  }

	inline char nlabel() { return (char) ( v & 0x000000FF ); }

};

#endif

