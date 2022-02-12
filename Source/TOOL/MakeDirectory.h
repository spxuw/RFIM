#ifndef MAKEDIRECTORIES_H
#define MAKEDIRECTORIES_H

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <ctime>
#include <limits.h>
#include <iterator>
#include <sstream>

#include <list>
#include <vector> 
#include <queue> 
#include <stack> 
#include <map> 
#include <algorithm>

#include <string>

using namespace std;


void MakeDirectory();
void MakeDirectory_Default();

void MakeDirectory(int D, int L, double R, int seed);
void MakeDirectory(int N, double c, double R, int seed);
void MakeDirectory(string file);
void MakeDirectories(string file);


#endif // MAKEDIRECTORIES_H

