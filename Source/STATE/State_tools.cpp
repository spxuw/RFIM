#include "State.h"

using namespace std;

//////////////////////////////////////////////////////////////////////////////////////////////////////
int State::GetLoc(const int* coords) const
{
    int loc=0;
    loc+=coords[D-1];
    for(int i=D-2;i>=0;i--)
    {
	loc+=coords[i]*stride[i];
    }
    return loc;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////



//////////////////////////////////////////////////////////////////////////////////////////////////////
void State::GetCoords(int loc, int* newCoords) const
{
    long remainder;
    remainder = loc%stride[0];
    newCoords[0] = loc/stride[0];

    for(int i=1;i<D;i++) {
	newCoords[i] = remainder/stride[i];
	remainder %= stride[i];
    }
}
//////////////////////////////////////////////////////////////////////////////////////////////////////



//////////////////////////////////////////////////////////////////////////////////////////////////////
void State::GetNeighbors(int loc) const  
{
    int* coords = new int[D];
    GetCoords(loc, coords);        

    for(int i=0;i<D;i++) {
	for(int j=0;j<D;j++) {
	    if(j==i) {
		if(coords[j]==size[j]-1)//if(coords[j]==L-1)
		{
		    neighborLocs[2*j]=loc-(size[j]-1)*stride[j];
		    neighborLocs[2*j+1]=loc-stride[j];
		}
		else if(coords[j]==0)
		{
		    neighborLocs[2*j]=loc+stride[j];
		    neighborLocs[2*j+1]=loc+(size[j]-1)*stride[j];
		}
		else
		{
		    neighborLocs[2*j]=loc+stride[j];
		    neighborLocs[2*j+1]=loc-stride[j];
		}
	    }
	}
    }

    delete [] coords;

}
//////////////////////////////////////////////////////////////////////////////////////////////////////





//////////////////////////////////////////////////////////////////////////////////////////////////////
void Statistics(const vector<double>& data, double& ave, double& var)
{
  int k;
  int n = data.size();

  ////////// average value
  double s = 0.0;
  for (k=0;k<n;k++)    {
    s +=  data[k];
  }
  ave = s/n;
 
  ////////// variance 
  double ep = 0.0;
  var = 0.0;
  for (k=0;k<n;k++) {
    ep += (s=data[k]-ave);
    var += s*s;
  }
  var=(var-ep*ep/n)/(n-1); // Corrected two-pass algorithm
}
///////////////////////////////////////////////////////////////////////////////


#include <sstream>
void Char2String(char& c, string& s) {
  stringstream ss;
  ss << c;
  ss >> s;
}


string Int2String(int i) {
  string s;
  stringstream ss;
  ss << i;
  ss >> s;
  return s;
}
