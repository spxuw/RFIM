#include <iostream>
#include <fstream>
#include <vector>
#include "Histogram.h"

using namespace std;

///////////////////////////////////////////////////////////////////////////////////////
// This universal routine will read a single column from a data file to calculate its histogram.
// 08/20/12
///////////////////////////////////////////////////////////////////////////////////////
int main_histogram_single_column(int argc, char** argv)
//int main(int argc, char** argv)
{
  if(argc!=4 && argc!=5) {
    cout << "\nNumber of arguments is incorrect.\n";
    cout << "\nCommand format : histogram file col type [ratio or Nbins]\n ";
    cout << " type==1 : linear-bin for discrete dist.  binwidth=1).\n ";
    cout << " type==2 : linear-bin for continuous dist. Nbins should be given.\n ";
    cout << " type==3 : log-bin for discrete dist. ratio should be given.\n ";
    cout << " type==4 : log-bin for continuous dist. ratio should be given.\n ";
    cout << " type==5 : discrete dist (X=number) without binning.\n ";
    cout << " type==6 : discrete dist (X=string) without binning.\n\n ";
    exit(0);
  }
   
  string file   = argv[1];
  int    col    = atoi(argv[2]);
  int    type   = atoi(argv[3]); // 0: linear; 1: log; 2: discrete;
  int    Nbins  = 100;
  double ratio  = 1.5;

  if(type==2) {
    if(argc!=5) {
      cout << "Nbin must be given.\n";
      exit(0);
    }
    Nbins = atoi(argv[4]);
    cout << "Nbins = " << Nbins << endl;
  }

  if(type==3 || type==4) {
    if(argc!=5) {
      cout << "ratio must be given.\n";
      exit(0);
    }
    ratio = atof(argv[4]);
    cout << "ratio = " << ratio << endl;
  }

  char fname [256];
  sprintf(fname, "%s", file.c_str());
  cout << fname << endl; //test

  char cmd [256];
  char ftemp [256];
  sprintf(ftemp, "%s.col%d", fname, col);
  //sprintf(cmd, "awk '$%d !~ /nan/ {print $%d}' %s > %s", col, col, fname, ftemp);
  sprintf(cmd, "awk '{if($%d !~ /nan/ && $%d !~ /inf/) print $%d}' %s > %s", col, col, col, fname, ftemp);
  system(cmd);


  ifstream fin(ftemp, ios::in);         
  if (!fin) {cout << "cannot open file: " << ftemp << " for read\n"; exit (0);}

  char filename [256];
  sprintf(filename, "%s.histogram", ftemp);


  if(type<6) {
    vector<double> X;    
    double x;
    while (fin >> x) {
      X.push_back(x);
    }
    fin.close();
    vector<int> Xint(X.begin(), X.end());
      
    switch (type) {
    case 1 :  // linear-bin, discrete 
      GetHistogram(Xint, filename);
      break;
    case 2 :  // linear-bin, continuous
      GetHistogram(X, Nbins, filename);
      break;
    case 3 :  // log-bin, discrete 
      GetHistogram(Xint, ratio, filename);
      break;
    case 4 :  // log-bin, continuous
      GetHistogram(X, ratio, filename);
      break;
    case 5 :  // discrete, without binning
      GetDiscreteDistribution(X, filename);
      break;
    default:
      cout << "type must be 1,2,3,4,5\n";
      exit(0);
    }

    char cmd[256];
    sprintf(cmd,"showps.histogram %s %d %e &", filename, type, ratio);
    system(cmd);
  }
  else if(type==6){
    vector<string> X;    
    char x[512];

    while(fin.getline(x,512)) {
      //cout << x << endl; //test
      X.push_back(x);
    }
    fin.close();

    GetDiscreteDistribution(X, filename);

    char cmd[256];
    sprintf(cmd,"gnuplot %s.discrete.plotscript; open %s.discrete.eps &", filename, filename);
    system(cmd);
  }



  sprintf(cmd,"rm %s", ftemp);
  system(cmd);

  exit(1);

}



///////////////////////////////////////////////////////////////////////////////////////
// This universal routine will read three columns from a data file to calculate 
// its 2D histogram, i.e., 2D binned average.
// 08/21/12
///////////////////////////////////////////////////////////////////////////////////////
int main_2D_histogram_logX_linearY_percentiles(int argc, char** argv)
//int main(int argc, char** argv)
{
  if(argc!=7) {
    cout << "\nNumber of arguments is incorrect.\n";
    cout << "\nCommand format : binnedaverage file colx coly colz ratiox Nbins_y\n ";
    cout << " log-bin for x; linear-bin for y with given # of bins.\n ";
    exit(0);
  }
   
  string file   = argv[1];
  int    colx   = atoi(argv[2]);
  int    coly   = atoi(argv[3]);
  int    colz   = atoi(argv[4]);

  double ratiox  = atof(argv[5]);
  int    Nbins_y = atoi(argv[6]);


  char fname [256];
  sprintf(fname, "%s", file.c_str());

  char cmd [256];
  char ftemp [256];
  sprintf(ftemp, "%s.col%d-col%d-col%d", fname, colx, coly, colz);
  sprintf(cmd, "awk '{print $%d, $%d, $%d}' %s > %s", colx, coly, colz, fname, ftemp);
  system(cmd);

  ifstream fin(ftemp, ios::in);         
  if (!fin) {cout << "cannot open file: " << ftemp << " for read\n"; exit (0);}

  char filename [256];
  sprintf(filename, "%s.2Dbinnedaverage", ftemp);


  vector<double> X;   
  vector<double> Y;  
  vector<double> Z;  

  double x;
  double y;
  double z;

  while (fin >> x >> y >> z) {
    X.push_back(x);
    Y.push_back(y);
    Z.push_back(z);
  }
  fin.close();
    
  //Get2DHistogram(X, Y, ratiox, ratioy, filename, Z);
  //Get2DHistogram_u(X, Y, ratiox, yinterval, filename, Z);
  Get2DHistogram_u2(X, Y, ratiox, Nbins_y, filename, Z);

 
  exit(1);

}
///////////////////////////////////////////////////////////////////////////////////////





///////////////////////////////////////////////////////////////////////////////////////
// This universal routine will read three columns from a data file to calculate 
// its 2D histogram, i.e., 2D binned average.
// 08/21/12
///////////////////////////////////////////////////////////////////////////////////////
int main_2D_histogram_logX_linearY(int argc, char** argv)
//int main(int argc, char** argv)
{
  if(argc!=7) {
    cout << "\nNumber of arguments is incorrect.\n";
    cout << "\nCommand format : binnedaverage file colx coly colz ratiox yinterval\n ";
    cout << " log-bin for x; linear-bin for y with given yinterval.\n ";
    exit(0);
  }
   
  string file   = argv[1];
  int    colx   = atoi(argv[2]);
  int    coly   = atoi(argv[3]);
  int    colz   = atoi(argv[4]);

  double ratiox    = atof(argv[5]);
  double yinterval = atof(argv[6]);


  char fname [256];
  sprintf(fname, "%s", file.c_str());

  char cmd [256];
  char ftemp [256];
  sprintf(ftemp, "%s.col%d-col%d-col%d", fname, colx, coly, colz);
  sprintf(cmd, "awk '{print $%d, $%d, $%d}' %s > %s", colx, coly, colz, fname, ftemp);
  system(cmd);

  ifstream fin(ftemp, ios::in);         
  if (!fin) {cout << "cannot open file: " << ftemp << " for read\n"; exit (0);}

  char filename [256];
  sprintf(filename, "%s.2Dbinnedaverage", ftemp);


  vector<double> X;   
  vector<double> Y;  
  vector<double> Z;  

  double x;
  double y;
  double z;

  while (fin >> x >> y >> z) {
    X.push_back(x);
    Y.push_back(y);
    Z.push_back(z);
  }
  fin.close();
    
  //Get2DHistogram(X, Y, ratiox, ratioy, filename, Z);
  Get2DHistogram_u(X, Y, ratiox, yinterval, filename, Z);

  sprintf(cmd,"./showps.2Dhistogram_binnedaverage_logX_linearY %s %e %e &", filename, ratiox, yinterval);
  system(cmd);
  
 
  exit(1);

}
///////////////////////////////////////////////////////////////////////////////////////





///////////////////////////////////////////////////////////////////////////////////////
// This universal routine will read three columns from a data file to calculate 
// its 2D histogram, i.e. 2D binned average.
// 08/21/12
///////////////////////////////////////////////////////////////////////////////////////
int main_2D_histogram(int argc, char** argv)
//int main(int argc, char** argv)
{
  if(argc!=7) {
    cout << "\nNumber of arguments is incorrect.\n";
    cout << "\nCommand format : binnedaverage file colx coly colz ratiox ratioy\n ";
    cout << " type==4 : log-bin for continuous dist. ratiox, ratioy should be given.\n ";
    exit(0);
  }
   
  string file   = argv[1];
  int    colx   = atoi(argv[2]);
  int    coly   = atoi(argv[3]);
  int    colz   = atoi(argv[4]);

  double ratiox = atof(argv[5]);
  double ratioy = atof(argv[6]);


  char fname [256];
  sprintf(fname, "%s", file.c_str());

  char cmd [256];
  char ftemp [256];
  sprintf(ftemp, "%s.col%d-col%d-col%d", fname, colx, coly, colz);
  sprintf(cmd, "awk '{print $%d, $%d, $%d}' %s > %s", colx, coly, colz, fname, ftemp);
  system(cmd);


  ifstream fin(ftemp, ios::in);         
  if (!fin) {cout << "cannot open file: " << ftemp << " for read\n"; exit (0);}

  char filename [256];
  sprintf(filename, "%s.2Dbinnedaverage", ftemp);


  vector<double> X;   
  vector<double> Y;  
  vector<double> Z;  

  double x;
  double y;
  double z;

  while (fin >> x >> y >> z) {
    X.push_back(x);
    Y.push_back(y);
    Z.push_back(z);
  }
  fin.close();
    
  Get2DHistogram(X, Y, ratiox, ratioy, filename, Z);

  sprintf(cmd,"./showps.2Dhistogram_binnedaverage %s %e %e &", filename, ratiox, ratioy);
  system(cmd);
  
 
  exit(1);

}
///////////////////////////////////////////////////////////////////////////////////////




/*
///////////////////////////////////////////////////////////////////////////////////////
// This universal routine will read two columns from a data file to calculate 
// its binned average.
// 08/21/12
///////////////////////////////////////////////////////////////////////////////////////
int main1_binned_average(int argc, char** argv)
{
  if(argc!=5 && argc!=6) {
    cout << "\nNumber of arguments is incorrect.\n";
    cout << "\nCommand format : binnedaverage file colx coly  type [ratio or Nbins]\n ";
    cout << " type==1 : linear-bin for discrete dist.  binwidth=1).\n ";
    cout << " type==2 : linear-bin for continuous dist. Nbins should be given.\n ";
    cout << " type==3 : log-bin for discrete dist. ratio should be given.\n ";
    cout << " type==4 : log-bin for continuous dist. ratio should be given.\n ";
    cout << " type==5 : discrete dist (X=number) without binning.\n ";
    exit(0);
  }
   
  string file   = argv[1];
  int    colx   = atoi(argv[2]);
  int    coly   = atoi(argv[3]);
  int    type   = atoi(argv[4]); // 0: linear; 1: log; 2: discrete;
  int    Nbins  = 100;
  double ratio  = 1.5;

  if(type==2) {
    if(argc!=6) {
      cout << "Nbin must be given.\n";
      exit(0);
    }
    Nbins = atoi(argv[5]);
    cout << "Nbins = " << Nbins << endl;
  }

  if(type==3 || type==4) {
    if(argc!=6) {
      cout << "ratio must be given.\n";
      exit(0);
    }
    ratio = atof(argv[5]);
    cout << "ratio = " << ratio << endl;
  }

  char fname [256];
  sprintf(fname, "%s", file.c_str());

  char cmd [256];
  char ftemp [256];
  sprintf(ftemp, "%s.col%d-col%d", fname, colx, coly);
  sprintf(cmd, "awk '{print $%d, $%d}' %s > %s", colx, coly, fname, ftemp);
  system(cmd);


  ifstream fin(ftemp, ios::in);         
  if (!fin) {cout << "cannot open file: " << ftemp << " for read\n"; exit (0);}

  char filename [256];
  sprintf(filename, "%s.binnedaverage", ftemp);



  vector<double> X;   
  vector<double> Y;   
  double x;
  double y;
  while (fin >> x >> y) {
    X.push_back(x);
    Y.push_back(y);
  }
  fin.close();
  vector<int> Xint(X.begin(), X.end());
      
  switch (type) {
  case 1 :  // linear-bin, discrete 
    GetHistogram(Xint, filename, Y);
    break;
  case 2 :  // linear-bin, continuous
    GetHistogram(X, Nbins, filename, Y);
    break;
  case 3 :  // log-bin, discrete 
    GetHistogram(Xint, ratio, filename, Y);
    break;
  case 4 :  // log-bin, continuous
    GetHistogram(X, ratio, filename, Y);
    break;
  case 5 :  // discrete, without binning
    GetDiscreteDistribution(X, filename, Y);
    break;
  default:
    cout << "type must be 1,2,3,4,5\n";
    exit(0);
  }

  sprintf(cmd,"./showps.histogram_binnedaverage %s %d %e &", filename, type, ratio);
  system(cmd);
  
 
  exit(1);

}
///////////////////////////////////////////////////////////////////////////////////////






///////////////////////////////////////////////////////////////////////////////////////
// This universal routine will read two columns from a data file to 
// bin the x value, and count in each x-bin, how many y values are positive
// This is use to calculate the winner fraction in the etoro project
// 11/29/12
///////////////////////////////////////////////////////////////////////////////////////
int main_binned_fraction(int argc, char** argv)
//int main(int argc, char** argv)
{
  if(argc!=5 && argc!=6) {
    cout << "\nNumber of arguments is incorrect.\n";
    cout << "\nCommand format : binnedaverage file colx coly  type [ratio or Nbins]\n ";
    cout << " type==2 : linear-bin for continuous dist. Nbins should be given.\n ";
    cout << " type==4 : log-bin for continuous dist. ratio should be given.\n ";
    exit(0);
  }
   
  string file   = argv[1];
  int    colx   = atoi(argv[2]);
  int    coly   = atoi(argv[3]);
  int    type   = atoi(argv[4]); 
  int    Nbins  = 100;
  double ratio  = 1.5;

  if(type==2) {
    if(argc!=6) {
      cout << "Nbin must be given.\n";
      exit(0);
    }
    Nbins = atoi(argv[5]);
    cout << "Nbins = " << Nbins << endl;
  }
  
  if(type==4) {
    if(argc!=6) {
      cout << "ratio must be given.\n";
      exit(0);
    }
    ratio = atof(argv[5]);
    cout << "ratio = " << ratio << endl;
  }

  char fname [256];
  sprintf(fname, "%s", file.c_str());

  char cmd [256];
  char ftemp [256];
  sprintf(ftemp, "%s.col%d-col%d", fname, colx, coly);
  sprintf(cmd, "awk '{print $%d, $%d}' %s > %s", colx, coly, fname, ftemp);
  system(cmd);


  ifstream fin(ftemp, ios::in);         
  if (!fin) {cout << "cannot open file: " << ftemp << " for read\n"; exit (0);}

  char filename [256];
  sprintf(filename, "%s.binnedaverage", ftemp);

  vector<double> X;   
  vector<double> Y;   
  double x;
  double y;
  while (fin >> x >> y) {
    X.push_back(x);
    Y.push_back(y);
  }
  fin.close();
  vector<int> Xint(X.begin(), X.end());
      
  switch (type) {
  case 2 :  // linear-bin, continuous
    GetHistogram(X, Nbins, filename, Y);
    break;
  case 4 :  // log-bin, continuous
    GetHistogram(X, ratio, filename, Y);
    break;
  default:
    cout << "type must be 2 or 4\n";
    exit(0);
  }

  sprintf(cmd,"./showps.histogram_binnedaverage %s %d %e &", filename, type, ratio);
  system(cmd);
  
 
  exit(1);

}
///////////////////////////////////////////////////////////////////////////////////////
*/


///////////////////////////////////////////////////////////////////////////////////////
// This universal routine will read two columns from a data file to 
// bin the x value in log bin. Note that x could be negative.
// Then in each x-bin, calculate the average of y values. 
// This is use in the etoro project
// 12/24/13
///////////////////////////////////////////////////////////////////////////////////////
//int main1_binned_average_negative_log(int argc, char** argv)
int main(int argc, char** argv)
{
  if(argc!=5 && argc!=6) {
    cout << "\nNumber of arguments is incorrect.\n";
    cout << "\nCommand format : binnedaverage file colx coly  type [ratio or Nbins]\n ";
    cout << " type==1 : linear-bin for discrete dist.  binwidth=1).\n ";
    cout << " type==2 : linear-bin for continuous dist. Nbins should be given.\n ";
    cout << " type==3 : log-bin for discrete dist. ratio should be given.\n ";
    cout << " type==4 : log-bin for continuous dist. ratio should be given.\n ";
    cout << " type==5 : discrete dist (X=number) without binning.\n ";
    exit(0);
  }
   
  string file   = argv[1];
  int    colx   = atoi(argv[2]);
  int    coly   = atoi(argv[3]);
  int    type   = atoi(argv[4]); // 0: linear; 1: log; 2: discrete;
  int    Nbins  = 100;
  double ratio  = 1.5;

  if(type==2) {
    if(argc!=6) {
      cout << "Nbin must be given.\n";
      exit(0);
    }
    Nbins = atoi(argv[5]);
    cout << "Nbins = " << Nbins << endl;
  }

  if(type==3 || type==4) {
    if(argc!=6) {
      cout << "ratio must be given.\n";
      exit(0);
    }
    ratio = atof(argv[5]);
    cout << "ratio = " << ratio << endl;
  }

  char fname [256];
  sprintf(fname, "%s", file.c_str());

  char cmd [256];
  char ftemp [256];
  sprintf(ftemp, "%s.col%d-col%d", fname, colx, coly);
  sprintf(cmd, "awk '{print $%d, $%d}' %s > %s", colx, coly, fname, ftemp);
  system(cmd);


  ifstream fin(ftemp, ios::in);         
  if (!fin) {cout << "cannot open file: " << ftemp << " for read\n"; exit (0);}

  char filename [256];
  sprintf(filename, "%s.binnedaverage", ftemp);



  vector<double> X;   
  vector<double> Y;   
  double x;
  double y;
  while (fin >> x >> y) {
    X.push_back(x);
    Y.push_back(y);
  }
  fin.close();
  vector<int> Xint(X.begin(), X.end());
      
  switch (type) {
  case 1 :  // linear-bin, discrete 
    GetHistogram(Xint, filename, Y);
    break;
  case 2 :  // linear-bin, continuous
    GetHistogram(X, Nbins, filename, Y);
    break;
  case 3 :  // log-bin, discrete 
    GetHistogram(Xint, ratio, filename, Y);
    break;
  case 4 :  // log-bin, continuous
    GetHistogram(X, ratio, filename, Y);
    break;
  case 5 :  // discrete, without binning
    GetDiscreteDistribution(X, filename, Y);
    break;
  default:
    cout << "type must be 1,2,3,4,5\n";
    exit(0);
  }

  sprintf(cmd,"./showps.histogram_binnedaverage %s %d %e &", filename, type, ratio);
  system(cmd);
  
 
  exit(1);

}
///////////////////////////////////////////////////////////////////////////////////////



