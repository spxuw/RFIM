#ifndef CHILDREN_H
#define CHILDREN_H

#include <map>
#include <vector>
#include "pattern.h"
#include "cost1.h"

class children
{
private:
	map<int, pattern *> pmap;
	vector< pattern *> ps;
	pattern * par;

public:
	//constructor
	children(pattern *p ){ 
		par = p;
	}

	~children(){ pmap.clear(); ps.clear(); }

	void put(occur * oc, cost1 * c1, int s, bool f);

	inline int count( cost1 * c1) { return pmap.count(c1.intValue()) ; }

	void scan(int f);
};

#endif