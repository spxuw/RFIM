#pragma warning (disable:4786 4018 4267)
#include <vector>
#include <string>
#include <map>
using namespace std;

#include "register.h"
/*
int Register :: get_index(string s)
{	
	if( index.count(s) ) 
		return index[s];
	else 
		return -1;

}
*/

void Register :: allocate_register(int i, string s)
{
	if( i == totalRegisters  ) {
		annotations.push_back(s);
		registers.push_back(0);
		totalRegisters++;
	}
	else
		log("the register has been allocated");
}

void Register :: peg(int i)
{
	//cout << "peg" << i << " " << totalRegisters << endl;
	if( i >= 0  && i < totalRegisters) {
		registers[i]++;
	}
	else
		log("no such register for pegging");
}

long Register :: read(int i)
{
	if( i >= 0  && i < totalRegisters) {
		return registers[i];
	}
	else
		log("no such register for reading");

	return -1;
}

void Register :: print()
{
	for( int i=0; i< totalRegisters; i++){
		cout << annotations[i] ;
		cout << ": " << read(i) << endl;
	}
}


