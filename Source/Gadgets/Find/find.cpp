#include <iostream>
#include <string>
#include "find.h"

using namespace std;


bool find(Nbl& nbl, int j)
{
  for(Nbl::iterator p = nbl.begin(); p!= nbl.end(); p++) {
    if(j==(*p)) 
      return true;
  }
    
  return false;
}

bool find(set<int>& nbl, int j)
{
  for(set<int>::iterator p = nbl.begin(); p!= nbl.end(); p++) {
    if(j==(*p)) 
      return true;
  }
    
  return false;
}

bool find(vector<vector<int> >& Vv, vector<int>& v)
{
  for(vector<vector<int> >::iterator p = Vv.begin(); p!= Vv.end(); p++) {
    if(v==(*p)) 
      return true;
  }
    
  return false;
}

bool find(vector<vector<string> >& Vs, vector<string>& v)
{
  for(vector<vector<string> >::iterator p = Vs.begin(); p!= Vs.end(); p++) {
    if(v==(*p)) 
      return true;
  }
    
  return false;
}

bool find(vector< vector<char>  >& Vstate, vector<char>& state)
{
  for(vector< vector<char> >::iterator p = Vstate.begin(); p!= Vstate.end(); p++) {
    if(state==(*p)) 
      return true;
  }
    
  return false;
}

bool find(vector<int>& nbv, int j)
{
  int n = nbv.size();
  for(int i=0; i<n; i++)
    if(j==nbv[i]) 
      return true;

  return false;
}

bool find(vector<string>& nbs, string s)
{
  int n = nbs.size();
  for(int i=0; i<n; i++)
    if(s==nbs[i]) 
      return true;

  return false;
}

bool find(deque<int>& Q, int x)
{
  for(deque<int>::iterator q = Q.begin(); q!= Q.end(); q++) {
    if(x==(*q)) 
      return true;
  }
  
  return false;
}





bool find(vector<set<int> >& X, set<int>& j) {
  for(vector<set<int> >::iterator p = X.begin(); p!= X.end(); p++) {
    if(j==(*p)) 
      return true;
  }
  return false;
}

bool find(vector<list<int> >& X, list<int>& j) {
  for(vector<list<int> >::iterator p = X.begin(); p!= X.end(); p++) {
    if(j==(*p)) 
      return true;
  }
  return false;
}


int findpos(Nbl& nbl, int j)
{
  int pos=0;
  for(Nbl::iterator p = nbl.begin(); p!= nbl.end(); p++) {
    if(j==(*p)) 
      return pos;
    else
      pos++;
  }

  return -1;
}

int findpos(vector<int>& nbv, int j)
{
  int n = nbv.size();
  for(int i=0; i<n; i++)
    if(j==nbv[i]) 
      return i;
  
  return -1;
}

int findpos(vector<string>& nbs, string s)
{
  int n = nbs.size();
  for(int i=0; i<n; i++)
    if(s==nbs[i]) 
      return i;
  
  return -1;
}

int print(Nbl& nbl)
{
  Nbl::iterator p = nbl.begin();
  cout << *p;
  p++;
  for(; p!= nbl.end(); p++) {
    cout << "-->" << *p;
  }
  cout << endl;

  return 1;
}

int getnodeinNbl(Nbl& nbl, int index)
{
  int pos=0;
  for(Nbl::iterator p = nbl.begin(); p!= nbl.end(); p++, pos++) {
    if(pos==index) 
      return (*p);
  }
  
  return -1;
}

double getnodeinlist(list<double>& X, int index)
{
  int pos=0;
  for(list<double>::iterator p = X.begin(); p!= X.end(); p++, pos++) {
    if(pos==index) 
      return (*p);
  }

  return -1;
}

int getnodeinlist(list<int>& X, int index)
{
  int pos=0;
  for(list<int>::iterator p = X.begin(); p!= X.end(); p++, pos++) {
    if(pos==index) 
      return (*p);
  }

  return -1;
}
