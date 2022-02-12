/*********************************************************************
 * cost.h
 *
 *********************************************************************/

#ifndef COST_H
#define COST_H

#include <string>
#include "common.h"

class Cost
{
private:
	vector<GLTYPE> costStr;

public:
	Cost(){ }
	~Cost(){ }

	string getCost() { string s; for(int i=0; i< costStr.size(); i++) s+= costStr[i]; return s; }
	void    append(GLTYPE c) { costStr.push_back(c); }
	void    insert(int i, GLTYPE c){ costStr.insert((costStr.begin()+i), c); }
	int     size() {return costStr.size(); } 
	GLTYPE  node(int i) {return costStr[i]; } 
};

#endif

