#pragma warning (disable:4786 4018 4267)
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
using namespace std;
#include "myTimer.h"

/*
int MyTimer :: get_index(string  s)
{
	int i;
	//map<string, int> ::iterator ip = timers.begin();
	//for( ip; ip != timers.end(); ip++)
	//	cout << "looking for " << s << "  " << ip->first << endl;

	if( timers.count(s))
		i = timers[s];
	else
		i = -1;
	return i;
}
*/

void MyTimer :: clear_timer(int i)
{

	if( i>= 0 && i < totalTimers){	
		elapsed[i] = 0.0;
		start[i] = -1.0;
	}
	else
		log("no timer for clear");
	
}
void MyTimer :: allocate_timer(int i, string s)
{
	if( i == totalTimers ){ //not allocated yet
		elapsed.push_back(0.0);
		start.push_back(-1.0);
		totalTimers++;
		timers.push_back(s);
		//clear_timer(s);
	}
	else
		log("timer already there");
}

void MyTimer :: start_timer(int i)
{
	if( i>= 0 && i < totalTimers)		
		start[i] = get_currentTime();
	else
		log("no such timer to start: ");
}

double MyTimer :: get_currentTime()
{

#ifdef _WIN32
	double e = GetTickCount();
	double r = (e - zero)/1000;
#else
	struct timeval e ;
	int i = gettimeofday(&e, &tz);
	double r = (e.tv_sec - zero.tv_sec) + 1.0e-6*(e.tv_usec-zero.tv_usec);
#endif

	return r;
}

void MyTimer :: stop_timer(int i)
{
	double e = get_currentTime();

	if( i>= 0 && i < totalTimers){	
		if( start[i] >= 0 ) 
			elapsed[i] += e - start[i];

		start[i] = -1.0;
	}
	else
		log("no such timer to stop: ");	
}

double  MyTimer::value_timer(int i)
{
	if( i>= 0 && i < totalTimers)	
		return elapsed[i];
	else
		log("no such timer for evaluation: ");
	return 0;
}

void MyTimer::print()
{
	vector<string> ::iterator ip = timers.begin();
	for( ip; ip != timers.end(); ip++){
		cout << dec << *ip ; 
		//cout << "index is "   << ip-timers.begin() << endl;
		cout << dec << ": " << elapsed[ip-timers.begin()] << endl;
	}
}
