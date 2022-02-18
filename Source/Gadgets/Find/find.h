#ifndef _MYNR_H_
#define _MYNR_H_

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <string>
#include <list>
#include <set>
#include <stack>
#include <deque>
#include <limits>
#include <vector>



using namespace std;

typedef list<int> Nbl; // neighbor list

bool find(vector<set<int> >& Vset, set<int>& s);
bool find(vector< vector<char>  >& Vstate, vector<char>& state);
bool find(vector< vector<int> >& Vv, vector<int>& v);
bool find(vector< vector<string> >& Vs, vector<string>& v);

bool find(vector<set<int> >& X, set<int>& j);
bool find(vector<list<int> >& X, list<int>& j);
bool find(set<int>& nbl, int j);

bool find(Nbl& nbl, int j);
int findpos(Nbl& nbl, int j);

bool find(vector<int>& nbv, int j);
int findpos(vector<int>& nbv, int j);

bool find(vector<string>& nbs, string s);
int findpos(vector<string>& nbs, string s);

bool find(deque<int>& S, int x);
void print(deque<int>& Q);



int getnodeinNbl(Nbl& nbl, int j);
double getnodeinlist(list<double>& X, int index);
int getnodeinlist(list<int>& X, int index);

#endif
