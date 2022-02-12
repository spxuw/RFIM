#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cmath>
#include <vector>
#include <list>
#include <queue> 
#include <deque>
#include <map>
#include <stack>
#include <set>
#include <limits.h>
#include <stack> 
#include <time.h> 

using namespace std;

// My version
double Px_HypergeometricDistribution(int a, int b, int c, int d);
double pvalue_HypergeometricDistribution(int a, int b, int c, int d);


// version of  Mark Von Tress: http://lists.mcgill.ca/scripts/wa.exe?A2=ind9808b&L=stat-l&P=5654
void fishertest (double a, double b, double c, double d, double &leftpval, double &rightpval, double &twopval);
