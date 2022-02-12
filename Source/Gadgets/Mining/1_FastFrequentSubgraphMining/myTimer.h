#ifndef  MYTIMER
#define  MYTIMER
#include <map>
#include <string>
#include "common.h"

#ifndef _WIN32 
	#include <sys/time.h>
	#include <unistd.h> 
#else
	#include <windows.h>
#endif

using namespace std;

class MyTimer
{
private:
	vector<double> elapsed;
	vector<double> start;
	//map<string, int> timers;
	vector<string>  timers;
	int totalTimers;

#ifndef _WIN32 
	struct timeval zero;
	struct timezone tz; 
#else 
	double zero;
#endif

	double get_currentTime();
	//int    get_index(string  s);

public:
	MyTimer(){

#ifndef _WIN32 
	 int i = gettimeofday(&zero, &tz);
#else 
	 zero = GetTickCount();
#endif	
	 totalTimers = 0;

	}

	void clear_timer(int s);
	void allocate_timer(int  s, string s1);
	void start_timer(int  s);
	void stop_timer(int   s);
	double  value_timer(int  s);
	void print();
};

#endif
