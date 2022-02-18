#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <vector>
#include <set>
#include <algorithm>

using namespace std;

void GetDensityMap(const vector<double>& X, const vector<double>& Y, 
		   int Nbinx, int Nbiny, 
		   double xmin, double xmax, double ymin, double ymax, 
		   char* fname);






int main(int argc, char** argv)
{
  if(argc!=8) {
    cout << "\nNumber of arguments is incorrect.\n";
    cout << "\nCommand format : Get_DensityMap file Nbinx Nbiny xmin xmax ymin ymax \n ";
    exit(0);
  }
   
  string file   = argv[1];
  int    Nbinx  = atoi(argv[2]);
  int    Nbiny  = atoi(argv[3]);
  double xmin   = atof(argv[4]);
  double xmax   = atof(argv[5]);
  double ymin   = atof(argv[6]);
  double ymax   = atof(argv[7]);


  char fname [256];
  sprintf(fname, "%s", file.c_str());

  ifstream fin(fname, ios::in);         
  if (!fin) {cout << "cannot open file: " << fname << " for read\n"; exit (0);}

  vector<double> X;   
  vector<double> Y;  
  double x;
  double y;
  while (fin >> x >> y) {
    X.push_back(x);
    Y.push_back(y);
  }
  fin.close();
    
  
  char filename [256];
  sprintf(filename, "%s.densitymap", fname);
  GetDensityMap(X,Y,Nbinx,Nbiny,xmin,xmax,ymin,ymax,filename);

  char cmd [256];
  sprintf(cmd,"./showps.densitymap %s &", filename);
  system(cmd);
  
 
  exit(1);

}




class Bin2D
{
 public:
    double xbelow;    
    double xabove;    

    double ybelow;    
    double yabove;    

    double count;    
    double normalized_count;       
};




void GetDensityMap(const vector<double>& X, const vector<double>& Y, 
		   int Nbins_x, int Nbins_y, 
		   double Xmin, double Xmax, double Ymin, double Ymax, 
		   char* fname)
{
  int N = X.size();
    
  int Nbins = Nbins_x * Nbins_y; 
  cout << Nbins_x << "," << Nbins_y << endl; 

  Bin2D** Binarray = new Bin2D* [Nbins_y];
  for(int b_y=0; b_y<Nbins_y; b_y++) {
    Binarray[b_y] = new Bin2D [Nbins_x];
  }

  double dx = (Xmax-Xmin)/(double)Nbins_x;
  double dy = (Ymax-Ymin)/(double)Nbins_y;

  
  for(int b_y=0;b_y<Nbins_y;b_y++) {
    for(int b_x=0; b_x<Nbins_x; b_x++){   
    
      if(b_x==0) Binarray[b_y][b_x].xbelow = Xmin;
      else       Binarray[b_y][b_x].xbelow = Binarray[b_y][b_x-1].xabove;

      if(b_y==0) Binarray[b_y][b_x].ybelow = Ymin;
      else       Binarray[b_y][b_x].ybelow = Binarray[b_y-1][b_x].yabove;

      Binarray[b_y][b_x].xabove = Binarray[b_y][b_x].xbelow + dx;
      Binarray[b_y][b_x].yabove = Binarray[b_y][b_x].ybelow + dy;

      Binarray[b_y][b_x].count = 0;
    }
  }


  int b_y = 0;
  int b_x = 0;
  for(int i=0;i<N;i++)    {	
    if(X[i]>=Xmin && Y[i]>=Ymin && X[i]<=Xmax && Y[i]<=Ymax) {

      b_y = (int)floor((Y[i]-Ymin)/dy);
      b_x = (int)floor((X[i]-Xmin)/dx);

      Binarray[b_y][b_x].count++; 
    }
  }


  
  char filename[256];
  sprintf(filename,"%s", fname);
  ofstream fout(filename, ios_base::out);
    

  for(int b_y=0;b_y<Nbins_y;b_y++) {
    for(int b_x=0; b_x<Nbins_x; b_x++){   
      Binarray[b_y][b_x].normalized_count = Binarray[b_y][b_x].count/(double)N;
      fout << Binarray[b_y][b_x].xbelow + 0.5*dx << ' ';        
      fout << Binarray[b_y][b_x].ybelow + 0.5*dy << ' ';        
      fout << Binarray[b_y][b_x].normalized_count << endl;
    }
    fout << endl;
  }
  fout.close();

  for(int i=0;i<Nbins_y;i++) 
    delete [] Binarray[i];
  
  delete [] Binarray;

}


