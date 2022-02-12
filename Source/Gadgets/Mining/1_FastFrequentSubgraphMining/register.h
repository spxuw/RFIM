#ifndef REGISTER_H
#define REGISTER_H

#include <map>
#include <string>
#include "common.h"

class Register
{
private:
	int totalRegisters;
	vector<DLONG> registers;
	vector<string>   annotations; 

	//int get_index(string s);
public:

	Register(){ totalRegisters = 0; }
	~Register(){};
	void allocate_register(int i, string s);
	void peg(int i);
	//void peg(int s, int n){ for (int i=0; i< n; i++) peg(s); }
	long int read(int s);
	void print();

};
#endif
