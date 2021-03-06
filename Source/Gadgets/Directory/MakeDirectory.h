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


// Makes directory "data" and "data/average", both under Unix and Windows
void MakeDirectories();
void MakeDirectories(double m);
void MakeDirectories(string file);
void MakeDirectory(string path, int type);
void MakeDirectory(string path, double p);

void MakeDirectory(string file);
void MakeDirectories(char* file);
void MakeDirectories(char* file, int seed);

//void MakeDirectories(int L);

#endif // MAKEDIRECTORIES_H

