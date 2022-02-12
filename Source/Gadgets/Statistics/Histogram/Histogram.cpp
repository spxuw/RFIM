#include "Histogram.h"

// P(k) = C(N-1,k) p^k (1-p)^(N-1-k)
double Binomial(int N, double p, int k)
{
  double Pk = 1;
  for (int i=1; i<=k ; i++)
    Pk *= (N-i)*p/(double)i;

  for (int i=1; i<=N-1-k; i++)
    Pk *= (1-p);

  return Pk;
}


///////////////////////////////////////////////////////////////////////////
// case 3: log bin for discrete distribution
void GetHistogram(const vector<int>& X, double ratio, char* fname)
{
  int N = X.size();
    
  double Xmax, Xmin, Xave, Xvar, Xerror;
  Statistics(X, Xmax, Xmin, Xave, Xvar, Xerror);

  if(Xmin==0) {
    cout << "Note that Xmin = 0.\n";
    Xmin=1;
  }

  int Nbins = (int)(floor(log(Xmax/Xmin)/log(ratio))) + 1; 
  Bin* Binarray = new Bin [Nbins];
  for(int i=0; i<Nbins; i++)     // initialize these bins
    {
      if(i==0) 
	Binarray[i].below = Xmin;
      else 
	Binarray[i].below = Binarray[i-1].above;
      Binarray[i].above = Binarray[i].below*ratio;
      Binarray[i].count = 0;
    }// end of initializing these bins


  //cout << "test" << endl; 
  int b = 0;
  for(int i=0;i<N;i++)  {	
    if(X[i]>0) { // Note that we don't consider isolated nodes (X[i]=0)
      b = (int)( floor( log(X[i]/Xmin)/log(ratio) ) ); // get the index of the areabin
      Binarray[b].count++; // for this bin, update it with new data point 
    }
  }


  int Dmax = 0;
  int Dsum = 0;
  for(int b=0;b<Nbins;b++)   {
    Dsum += Binarray[b].count;
    if(Binarray[b].count>Dmax)
      Dmax = Binarray[b].count;
  }
  //cout << "\n Dsum = " << Dsum;
  //cout << "\n Dmax = " << Dmax << endl;

  // save the distribution to file
  char filename[256];
  sprintf(filename,"%s", fname);
  ofstream fout(filename, ios_base::out);
    
  for(int b=0; b<Nbins; b++)    {
    Binarray[b].normalized_count = Binarray[b].count/(double)N;
    
    fout << b << ' ';
    fout.width(8);   fout << Binarray[b].below << ' ';
    fout.width(8);   fout << Binarray[b].above << ' ';
    fout.width(8);   fout << Binarray[b].count << ' ';
    fout.width(8);   fout << Binarray[b].normalized_count << ' ';
    fout << endl;
    
  }// end of writing data for each bin
  fout.close();

  delete [] Binarray;


}
///////////////////////////////////////////////////////////////////////////





///////////////////////////////////////////////////////////////////////////
// case 3: log bin for discrete distribution
void GetHistogram_Accumulatively(const vector<int>& X, double ratio, char* label, char* fname)
{
  int N = X.size();
    
  double Xmax, Xmin, Xave, Xvar, Xerror;
  Statistics(X, Xmax, Xmin, Xave, Xvar, Xerror);

  if(Xmin==0) {
    cout << "Note that Xmin = 0.\n";
    Xmin=1;
  }

  int Nbins = (int)(floor(log(Xmax/Xmin)/log(ratio))) + 1; 
  Bin* Binarray = new Bin [Nbins];
  for(int i=0; i<Nbins; i++)     // initialize these bins
    {
      if(i==0) 
	Binarray[i].below = Xmin;
      else 
	Binarray[i].below = Binarray[i-1].above;
      Binarray[i].above = Binarray[i].below*ratio;
      Binarray[i].count = 0;
    }// end of initializing these bins


  //cout << "test" << endl; 
  int b = 0;
  for(int i=0;i<N;i++)  {	
    if(X[i]>0) { // Note that we don't consider isolated nodes (X[i]=0)
      b = (int)( floor( log(X[i]/Xmin)/log(ratio) ) ); // get the index of the areabin
      Binarray[b].count++; // for this bin, update it with new data point 
    }
  }


  int Dmax = 0;
  int Dsum = 0;
  for(int b=0;b<Nbins;b++)   {
    Dsum += Binarray[b].count;
    if(Binarray[b].count>Dmax)
      Dmax = Binarray[b].count;
  }
  //cout << "\n Dsum = " << Dsum;
  //cout << "\n Dmax = " << Dmax << endl;

  // save the distribution to file
  char filename[256];
  sprintf(filename,"%s", fname);
  ofstream fout(filename, ios_base::app);
  fout << "\n#########  " << label << "  ##########" << endl; 
   
  for(int b=0; b<Nbins; b++)    {
    Binarray[b].normalized_count = Binarray[b].count/(double)N;
    
   
    fout << b << ' ';
    fout.width(8);   fout << Binarray[b].below << ' ';
    fout.width(8);   fout << Binarray[b].above << ' ';
    fout.width(8);   fout << Binarray[b].count << ' ';
    fout.width(8);   fout << Binarray[b].normalized_count << ' ';
    fout << endl;
    
  }// end of writing data for each bin
  fout.close();

  delete [] Binarray;


}
///////////////////////////////////////////////////////////////////////////




///////////////////////////////////////////////////////////////////////////
// log bin for discrete distribution, for GKK model 
void GetHistogram(const vector<int>& X, double ratio, char* fname, double& Kmean, double& K2mean)
{
  int N = X.size();
    
  double Xmax, Xmin, Xave, Xvar, Xerror;
  Statistics(X, Xmax, Xmin, Xave, Xvar, Xerror);

  Kmean  = Xave;
  K2mean = Xvar + Xave*Xave;
    
  //cout << "\n Kmean = " << Kmean << endl; //test

  if(Xmin==0) {
    cout << "Note that Xmin = 0.\n";
    Xmin=1;
  }

  int Nbins = (int)(floor(log(Xmax/Xmin)/log(ratio))) + 1; 
  Bin* Binarray = new Bin [Nbins];
  for(int i=0; i<Nbins; i++) {    // initialize these bins
    if(i==0) 
      Binarray[i].below = Xmin;
    else 
      Binarray[i].below = Binarray[i-1].above;
    Binarray[i].above = Binarray[i].below*ratio;
    Binarray[i].count = 0;
  }// end of initializing these bins


  //cout << "test" << endl; 
  int b = 0;
  for(int i=0;i<N;i++)    {	
    if(X[i]>0) {// Note that we don't consider isolated nodes (X[i]=0)
      b = (int)( floor( log(X[i]/Xmin)/log(ratio) ) ); // get the index of the areabin
      Binarray[b].count++; // for this bin, update it with new data point 
    }
  }


  int Dmax = 0;
  int Dsum = 0;
  for(int b=0;b<Nbins;b++)   {
    Dsum += Binarray[b].count;
    if(Binarray[b].count>Dmax)
      Dmax = Binarray[b].count;
  }
  //cout << "\n Dsum = " << Dsum;
  //cout << "\n Dmax = " << Dmax << endl;

  // save the distribution to file
  char filename[256];
  sprintf(filename,"%s", fname);
  ofstream fout(filename, ios_base::out);
    
  for(int b=0; b<Nbins; b++)  {
    Binarray[b].normalized_count = Binarray[b].count/(double)N;

    fout << b << ' ';
    fout.width(8);   fout << Binarray[b].below << ' ';
    fout.width(8);   fout << Binarray[b].above << ' ';
    fout.width(8);   fout << Binarray[b].count << ' ';
    fout.width(8);   fout << Binarray[b].normalized_count << ' ';
    fout << endl;
    
  }// end of writing data for each bin
  fout.close();

  delete [] Binarray;


}
///////////////////////////////////////////////////////////////////////////




///////////////////////////////////////////////////////////////////////////
// log bin for discrete distribution
void GetHistogram(const vector<int>& X, double ratio, char* fname, double& Kmean, double& K2mean, double& Kmax)
{
  int N = X.size();
    
  double Xmax, Xmin, Xave, Xvar, Xerror;
  Statistics(X, Xmax, Xmin, Xave, Xvar, Xerror);

  Kmax   = Xmax;
  Kmean  = Xave;
  K2mean = Xvar + Xave*Xave;
    
  //cout << "\n Kmean = " << Kmean << endl; //test

  if(Xmin==0) {
    cout << "Note that Xmin = 0.\n";
    Xmin=1;
  }

  int Nbins = (int)(floor(log(Xmax/Xmin)/log(ratio))) + 1; 
  Bin* Binarray = new Bin [Nbins];
  for(int i=0; i<Nbins; i++) {    // initialize these bins
    if(i==0) 
      Binarray[i].below = Xmin;
    else 
      Binarray[i].below = Binarray[i-1].above;
    Binarray[i].above = Binarray[i].below*ratio;
    Binarray[i].count = 0;
  }// end of initializing these bins
  

  //cout << "test" << endl; 
  int b = 0;
  for(int i=0;i<N;i++)  {	
    if(X[i]>0) { // Note that we don't consider isolated nodes (X[i]=0)
      b = (int)( floor( log(X[i]/Xmin)/log(ratio) ) ); // get the index of the areabin
      Binarray[b].count++; // for this bin, update it with new data point 
    }
  }


  int Dmax = 0;
  int Dsum = 0;
  for(int b=0;b<Nbins;b++)   {
    Dsum += Binarray[b].count;
    if(Binarray[b].count>Dmax)
      Dmax = Binarray[b].count;
  }
  //cout << "\n Dsum = " << Dsum;
  //cout << "\n Dmax = " << Dmax << endl;

  // save the distribution to file
  char filename[256];
  sprintf(filename,"%s", fname);
  ofstream fout(filename, ios_base::out);
    
  for(int b=0; b<Nbins; b++)    {
    Binarray[b].normalized_count = Binarray[b].count/(double)N;

    fout << b << ' ';
    fout.width(8);   fout << Binarray[b].below << ' ';
    fout.width(8);   fout << Binarray[b].above << ' ';
    fout.width(8);   fout << Binarray[b].count << ' ';
    fout.width(8);   fout << Binarray[b].normalized_count << ' ';
    fout << endl;
  }// end of writing data for each bin
  fout.close();

  delete [] Binarray;


}
///////////////////////////////////////////////////////////////////////////





///////////////////////////////////////////////////////////////////////////
// linear bin for discrete distribution
void GetHistogram(const vector<int>& X)
{
  int N = X.size();
    
  double Xmax, Xmin, Xave, Xvar, Xerror;
  Statistics(X, Xmax, Xmin, Xave, Xvar, Xerror);

  double Xbinwidth  = 1;//(Xmax-Xmin)/Nbins; // h bin's width
  int Nbins = Xmax-Xmin+1;
  vector<int> D(Nbins,0);                // the bin vector 
  for(int i=0;i<N;i++)    {	
    int b = (int)floor((X[i]-Xmin)/Xbinwidth);
    if(b>=Nbins)      {
      cout << "\n (X[i]-Xmin)/Xbinwidth = " << (X[i]-Xmin)/Xbinwidth  <<  " b = " << b << endl; //test
      b = Nbins-1;
    }
    if(b<0)	{
      cout << "\n (X[i]-Xmin)/Xbinwidth = " << (X[i]-Xmin)/Xbinwidth  <<  " b = " << b << endl; //test
      b = 0;
    }

    D[b]++;    // get the histogram 
  }


  int Dmax = 0;
  int Dsum = 0;
  for(int b=0;b<Nbins;b++)    {
    Dsum += D[b];
    if(D[b]>Dmax)
      Dmax = D[b];
  }
  //cout << "\n Dsum = " << Dsum;
  //cout << "\n Dmax = " << Dmax << endl;


  // save the distribution to file
  char filename[256];
  sprintf(filename,"Dist.dat");
  ofstream fout(filename, ios_base::out);
    
  for(int b=0; b<Nbins; b++)    {
    //fout << Xmin+(b)*Xbinwidth << ' ' << D[b]/(double)Dsum << endl;
    //fout << Xmin+(b+0.5)*Xbinwidth << ' ' << D[b]/(double)N << endl;
    //fout << Xmin+(b+0.5)*Xbinwidth << ' ' << D[b]/(double)(N*Xbinwidth) << endl; // this way we get a normalized histogram, i.e. including a total surface area equal to unity
    fout << Xmin+(b+0.5)*Xbinwidth << ' ' << D[b]/(double)N << ' ' << D[b]/(double)(N*Xbinwidth) << endl; 

  }
  fout.close();

}
///////////////////////////////////////////////////////////////////////////



///////////////////////////////////////////////////////////////////////////
// case 1: linear bin for discrete distribution
void GetHistogram(const vector<int>& X, char* fname)
{
  int N = X.size();
    
  double Xmax, Xmin, Xave, Xvar, Xerror;
  Statistics(X, Xmax, Xmin, Xave, Xvar, Xerror);

  double Xbinwidth  = 1;//(Xmax-Xmin)/Nbins; // h bin's width
  int Nbins = Xmax-Xmin+1;

  vector<int> D(Nbins,0);                // the bin vector 
  for(int i=0;i<N;i++)    {	
    int b = (int)floor((X[i]-Xmin)/Xbinwidth);
    if(b>=Nbins)	{
      cout << "\n (X[i]-Xmin)/Xbinwidth = " << (X[i]-Xmin)/Xbinwidth  <<  " b = " << b << endl; //test
      b = Nbins-1;
    }
    if(b<0)	{
      cout << "\n (X[i]-Xmin)/Xbinwidth = " << (X[i]-Xmin)/Xbinwidth  <<  " b = " << b << endl; //test
      b = 0;
    }

    D[b]++;    // get the histogram 
  }


  int Dmax = 0;
  int Dsum = 0;
  for(int b=0;b<Nbins;b++)    {
    Dsum += D[b];
    if(D[b]>Dmax)
      Dmax = D[b];
  }
  //cout << "\n Dsum = " << Dsum;
  //cout << "\n Dmax = " << Dmax << endl;


  // save the distribution to file
  char filename[256];
  sprintf(filename,"%s", fname);
  ofstream fout(filename, ios_base::out);

  
  for(int b=0; b<Nbins; b++)    {
    //fout << Xmin+(b)*Xbinwidth << ' ' << D[b]/(double)Dsum << endl;
    if(D[b]>0)
      //fout << Xmin+(b)*Xbinwidth << "     " << D[b]/(double)N << "     " << D[b] << endl;
      //fout << Xmin+(b+0.5)*Xbinwidth << ' ' << D[b]/(double)(N*Xbinwidth) << endl; // this way we get a normalized histogram, i.e. including a total surface area equal to unity
    fout << Xmin+(b)*Xbinwidth << "     " << D[b]/(double)N << "     " << D[b]/(double)(N*Xbinwidth) << endl; 
    
  }
  fout.close();

}
///////////////////////////////////////////////////////////////////////////



///////////////////////////////////////////////////////////////////////////
// linear bin for discrete distribution
void GetHistogram(const vector<int>& X, char* fname, double& Kmean, double& K2mean, double& Kmax)
{
  int N = X.size();
    
  double Xmax, Xmin, Xave, Xvar, Xerror;
  Statistics(X, Xmax, Xmin, Xave, Xvar, Xerror);

  Kmax   = Xmax;
  Kmean  = Xave;
  K2mean = Xvar + Xave*Xave;

  double Xbinwidth  = 1;//(Xmax-Xmin)/Nbins; // h bin's width
  int Nbins = Xmax-Xmin+1;

  vector<int> D(Nbins,0);                // the bin vector 
  for(int i=0;i<N;i++)    {	
    int b = (int)floor((X[i]-Xmin)/Xbinwidth);

    if(b>=Nbins)	{
      cout << "\n (X[i]-Xmin)/Xbinwidth = " << (X[i]-Xmin)/Xbinwidth  <<  " b = " << b << endl; //test
      b = Nbins-1;
    }
    if(b<0)	{
      cout << "\n (X[i]-Xmin)/Xbinwidth = " << (X[i]-Xmin)/Xbinwidth  <<  " b = " << b << endl; //test
      b = 0;
    }

    D[b]++;    // get the histogram 
  }


  int Dmax = 0;
  int Dsum = 0;
  for(int b=0;b<Nbins;b++)    {
    Dsum += D[b];
    if(D[b]>Dmax)
      Dmax = D[b];
  }
  //cout << "\n Dsum = " << Dsum;
  //cout << "\n Dmax = " << Dmax << endl;


  // save the distribution to file
  char filename[256];
  sprintf(filename,"%s", fname);
  ofstream fout(filename, ios_base::out);

  
  for(int b=0; b<Nbins; b++)    {
    //fout << Xmin+(b)*Xbinwidth << ' ' << D[b]/(double)Dsum << endl;
    //fout << Xmin+(b)*Xbinwidth << "     " << D[b]/(double)N << endl;
    //fout << Xmin+(b+0.5)*Xbinwidth << ' ' << D[b]/(double)(N*Xbinwidth) << endl; // this way we get a normalized histogram, i.e. including a total surface area equal to unity
    fout << Xmin+(b)*Xbinwidth << "     " << D[b]/(double)N << ' ' << D[b]/(double)(N*Xbinwidth) << endl; 
  }
  fout.close();

}
///////////////////////////////////////////////////////////////////////////



///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// linear bin for discrete distribution, for Erdos-Renyi model. p is the connection probability
void GetHistogram_ER(const vector<int>& X, double p, char* fname)
{
  int N = X.size();
  double Xmax, Xmin, Xave, Xvar, Xerror;
  Statistics(X, Xmax, Xmin, Xave, Xvar, Xerror);

  double Xbinwidth  = 1;//(Xmax-Xmin)/Nbins; // h bin's width
  int Nbins = Xmax-Xmin+1;
  vector<int> D(Nbins,0);                // the bin vector 
  for(int i=0;i<N;i++)    {	
    int b = (int)floor((X[i]-Xmin)/Xbinwidth);

    if(b>=Nbins)	{
      cout << "\n (X[i]-Xmin)/Xbinwidth = " << (X[i]-Xmin)/Xbinwidth  <<  " b = " << b << endl; //test
      b = Nbins-1;
    }
    if(b<0)	{
      cout << "\n (X[i]-Xmin)/Xbinwidth = " << (X[i]-Xmin)/Xbinwidth  <<  " b = " << b << endl; //test
      b = 0;
    }

    D[b]++;    // get the histogram 
  }


  int Dmax = 0;
  int Dsum = 0;
  for(int b=0;b<Nbins;b++)    {
    Dsum += D[b];
    if(D[b]>Dmax)
      Dmax = D[b];
  }
  //cout << "\n Dsum = " << Dsum;
  //cout << "\n Dmax = " << Dmax << endl;


  // save the distribution to file
  //char filename[256];
  //sprintf(filename,"Dist.dat");
  ofstream fout(fname, ios_base::out);
    
  for(int b=0; b<Nbins; b++)    {
    int k = Xmin + b;
    
    //fout << Xmin+(b)*Xbinwidth << ' ' << D[b]/(double)Dsum << endl;
    //fout << Xmin+(b)*Xbinwidth << "     " << D[b]/(double)N << "     " << Binomial(N, p, k) << endl;
    fout << Xmin+(b)*Xbinwidth << "     " 
	 << D[b]/(double)N     << "     " 
	 << Poisson(Xave, k)     << "     " 
	 << Binomial(N, p, k) << endl;
  }
  fout.close();
  
}
///////////////////////////////////////////////////////////////////////////





///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// linear bin for discrete distribution, for Erdos-Renyi model. p is the connection probability
void GetHistogram_ER(const vector<int>& X, double p, char* fname, double& Kmean, double& K2mean, double& Kmax)
{
  int N = X.size();

  double Xmax, Xmin, Xave, Xvar, Xerror;
  Statistics(X, Xmax, Xmin, Xave, Xvar, Xerror);

  Kmax   = Xmax;
  Kmean  = Xave;
  K2mean = Xvar + Xave*Xave;

  double Xbinwidth  = 1;//(Xmax-Xmin)/Nbins; // h bin's width
  int Nbins = Xmax-Xmin+1;
  vector<int> D(Nbins,0);                // the bin vector 
  for(int i=0;i<N;i++)    {	
    int b = (int)floor((X[i]-Xmin)/Xbinwidth);

    if(b>=Nbins)	{
      cout << "\n (X[i]-Xmin)/Xbinwidth = " << (X[i]-Xmin)/Xbinwidth  <<  " b = " << b << endl; //test
      b = Nbins-1;
    }
    if(b<0)	{
      cout << "\n (X[i]-Xmin)/Xbinwidth = " << (X[i]-Xmin)/Xbinwidth  <<  " b = " << b << endl; //test
      b = 0;
    }

    D[b]++;    // get the histogram 
  }


  int Dmax = 0;
  int Dsum = 0;
  for(int b=0;b<Nbins;b++)    {
    Dsum += D[b];
    if(D[b]>Dmax)
      Dmax = D[b];
  }
  //cout << "\n Dsum = " << Dsum;
  //cout << "\n Dmax = " << Dmax << endl;


  // save the distribution to file
  ofstream fout(fname, ios_base::out);
    

  for(int b=0; b<Nbins; b++)    {
    int k = Xmin + b;

    //fout << Xmin+(b)*Xbinwidth << ' ' << D[b]/(double)Dsum << endl;
    //fout << Xmin+(b)*Xbinwidth << "     " << D[b]/(double)N << "     " << Binomial(N, p, k) << endl;
    fout << Xmin+(b)*Xbinwidth << "     " 
	 << D[b]/(double)N     << "     " 
	 << Poisson(Xave, k)     << "     " 
	 << Binomial(N, p, k) << endl;
  }
  fout.close();

}
///////////////////////////////////////////////////////////////////////////





///////////////////////////////////////////////////////////////////////////
// linear bin for continuous distribution
void GetHistogram(const vector<int>& X, int Nbins)
{
  int N = X.size();
  double Xmax, Xmin, Xave, Xvar, Xerror;
  Statistics(X, Xmax, Xmin, Xave, Xvar, Xerror);

  double Xbinwidth  = (Xmax-Xmin)/Nbins; // h bin's width
  vector<int> D(Nbins,0);                // the bin vector 
  for(int i=0;i<N;i++)    {	
    int b = (int)floor((X[i]-Xmin)/Xbinwidth);

    if(b>=Nbins)	{
      cout << "\n (X[i]-Xmin)/Xbinwidth = " << (X[i]-Xmin)/Xbinwidth  <<  " b = " << b << endl; //test
      b = Nbins-1;
    }
    if(b<0)	{
      cout << "\n (X[i]-Xmin)/Xbinwidth = " << (X[i]-Xmin)/Xbinwidth  <<  " b = " << b << endl; //test
      b = 0;
    }

    D[b]++;    // get the histogram 
  }


  int Dmax = 0;
  int Dsum = 0;
  for(int b=0;b<Nbins;b++)    {
    Dsum += D[b];
    if(D[b]>Dmax)
      Dmax = D[b];
  }
  //cout << "\n Dsum = " << Dsum;
  //cout << "\n Dmax = " << Dmax << endl;


  // save the distribution to file
  char filename[256];
  sprintf(filename,"Dist.dat");
  ofstream fout(filename, ios_base::out);
    
  for(int b=0; b<Nbins; b++)    {
    //fout << Xmin+(b+0.5)*Xbinwidth << ' ' << D[b]/(double)Dmax << endl;
    //fout << Xmin+(b+0.5)*Xbinwidth << ' ' << D[b] << endl;
    //fout << Xmin+(b+0.5)*Xbinwidth << ' ' << D[b]/(double)(N*Xbinwidth) << endl; // this way we get a normalized histogram, i.e. including a total surface area equal to unity

    fout << Xmin+(b+0.5)*Xbinwidth << ' ' << D[b] << ' ' << D[b]/(double)(N*Xbinwidth) << endl; 
  }
  fout.close();

}
///////////////////////////////////////////////////////////////////////////




//////////////////////////////////////////////////////////////////////////////////
// case 4: log bin for continous distribution, e.g. the distribution of betweeness
void GetHistogram(const vector<double>& X, double ratio, char* fname)
{
  int N = X.size();
    
  double Xmax, Xmin, Xave, Xvar, Xerror;
  Statistics(X, Xmax, Xmin, Xave, Xvar, Xerror);

  double Xsmin;

  if(Xmin==0) {
    cout << "Note that Xmin = 0.\n";
    Xsmin = GetXsecondmin(X, Xmin);
    cout << "We get the second minimum value: " << Xsmin << endl;
  }
  else {
    Xsmin = Xmin;
  }

  int Nbins = int(floor(log(Xmax/Xsmin)/log(ratio))) + 1; 
  Bin* Binarray = new Bin [Nbins];
  for(int i=0; i<Nbins; i++){     // initialize these bins
    
    if(i==0) 
      Binarray[i].below = Xsmin;
    else 
      Binarray[i].below = Binarray[i-1].above;
    Binarray[i].above = Binarray[i].below*ratio;
    Binarray[i].count = 0;
  }// end of initializing these bins

  vector<double> Bin0;

  int b = 0;
  for(int i=0;i<N;i++)    {	
    if(X[i]>0) {// Note that we don't consider isolated nodes (X[i]=0)
      b = (int)( floor( log(X[i]/Xsmin)/log(ratio) ) ); // get the index of the areabin
      Binarray[b].count++; // for this bin, update it with new data point 
    }
    else {
      Bin0.push_back(X[i]);
    }
  }


  int Dmax = 0;
  int Dsum = 0;
  for(int b=0;b<Nbins;b++)    {
    Dsum += Binarray[b].count;
    if(Binarray[b].count>Dmax)
      Dmax = Binarray[b].count;
  }
  //cout << "\n Dsum = " << Dsum;
  //cout << "\n Dmax = " << Dmax << endl;

  // save the distribution to file
  char filename[256];
  sprintf(filename,"%s", fname);
  ofstream fout(filename, ios_base::out);


  if(Xmin<Xsmin) {
    fout << -1 << ' ';
    fout.width(8);   fout << 0 << ' ';
    fout.width(8);   fout << 0 << ' ';
    fout.width(8);   fout << Bin0.size() << ' ';
    fout.width(8);   fout << Bin0.size()/(double)N << ' ';
    fout << endl;
  }

  for(int b=0; b<Nbins; b++)    {
    Binarray[b].normalized_count = Binarray[b].count/(double)N;

    fout << b << ' ';
    fout.width(8);   fout << Binarray[b].below << ' ';
    fout.width(8);   fout << Binarray[b].above << ' ';
    fout.width(8);   fout << Binarray[b].count << ' ';
    fout.width(8);   fout << Binarray[b].normalized_count << ' ';
    fout << endl;

  }// end of writing data for each bin
  fout.close();

  delete [] Binarray;

}
///////////////////////////////////////////////////////////////////////////



///////////////////////////////////////////////////////////////////////////
// log bin for continous distribution, e.g. the distribution of Cs
void GetHistogram(const vector<double>& X, double ratio, char* fname, double& Xmin)
{
  int N = X.size();
    
  double Xmax, Xave, Xvar, Xerror;
  Statistics(X, Xmax, Xmin, Xave, Xvar, Xerror);
  //cout << "Xmin = " << Xmin << endl;

  int Nbins = int(floor(log(Xmax/Xmin)/log(ratio))) + 1; 
  Bin* Binarray = new Bin [Nbins];
  for(int i=0; i<Nbins; i++) {     // initialize these bins
    
    if(i==0) 
      Binarray[i].below = Xmin;
    else 
      Binarray[i].below = Binarray[i-1].above;
    Binarray[i].above = Binarray[i].below*ratio;
    Binarray[i].count = 0;
  }// end of initializing these bins

  int b = 0;
  for(int i=0;i<N;i++)    {	
    if(X[i]>0) // Note that we don't consider isolated nodes (X[i]=0)
      b = (int)( floor( log(X[i]/Xmin)/log(ratio) ) ); // get the index of the areabin
  
    Binarray[b].count++; // for this bin, update it with new data point 
  }


  int Dmax = 0;
  int Dsum = 0;
  for(int b=0;b<Nbins;b++)    {
    Dsum += Binarray[b].count;
    if(Binarray[b].count>Dmax)
      Dmax = Binarray[b].count;
  }
  //cout << "\n Dsum = " << Dsum;
  //cout << "\n Dmax = " << Dmax << endl;

  // save the distribution to file
  char filename[256];
  sprintf(filename,"%s", fname);
  ofstream fout(filename, ios_base::out);
    
  for(int b=0; b<Nbins; b++)    {
    Binarray[b].normalized_count = Binarray[b].count/(double)N;

    fout << b << ' ';
    fout.width(8);   fout << Binarray[b].below << ' ';
    fout.width(8);   fout << Binarray[b].above << ' ';
    fout.width(8);   fout << Binarray[b].count << ' ';
    fout.width(8);   fout << Binarray[b].normalized_count << ' ';
    fout << endl;

  }// end of writing data for each bin
  fout.close();

  delete [] Binarray;

}
///////////////////////////////////////////////////////////////////////////



/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// case 2: linear bin for continuous distribution, e.g. the eigenvalue distribution or the Csl distribution
void GetHistogram(const vector<double>& X, int Nbins, char* fname)
{
  int N = X.size();
  //cout << "Size= " << N << endl; //test

  double Xmax, Xmin, Xave, Xvar, Xerror;
  Statistics(X, Xmax, Xmin, Xave, Xvar, Xerror);

  // save the distribution to file
  char filename[256];
  sprintf(filename,"%s", fname);
  //ofstream fout(filename, ios_base::app);
  ofstream fout(filename, ios_base::out);

  double Xbinwidth  = (Xmax-Xmin)/Nbins; // h bin's width

  if(Xbinwidth>1e-10) {
    vector<int> D(Nbins,0);                // the bin vector 
    for(int i=0;i<N;i++) {	
      int b = (int)floor((X[i]-Xmin)/Xbinwidth);

      if(b>=Nbins){
	//cout << "\n (X[i]-Xmin)/Xbinwidth = " << (X[i]-Xmin)/Xbinwidth  <<  " b = " << b << endl; //test
	b = Nbins-1;
      }
      if(b<0){
	//cout << "\n (X[i]-Xmin)/Xbinwidth = " << (X[i]-Xmin)/Xbinwidth  <<  " b = " << b << endl; //test
	b = 0;
      }
      D[b]++;    // get the histogram 
    }


    int Dmax = 0;
    int Dsum = 0;
    for(int b=0;b<Nbins;b++){
      Dsum += D[b];
      if(D[b]>Dmax)
	Dmax = D[b];
    }
    //cout << "\n Dsum = " << Dsum;
    //cout << "\n Dmax = " << Dmax << endl;

    for(int b=0; b<Nbins; b++) {
      //fout << Xmin+(b+0.5)*Xbinwidth << ' ' << D[b]/(double)Dmax << endl;
      //fout << Xmin+(b+0.5)*Xbinwidth << ' ' << D[b]/(double)(N) << endl;
      //fout << Xmin+(b+0.5)*Xbinwidth << ' ' << D[b]/(double)(N*Xbinwidth) << endl; // this way we get a normalized histogram, i.e. including a total surface area equal to unity. Be careful if X is small, we will get a P(X) larger than 1. 

      fout << Xmin+(b+0.5)*Xbinwidth << ' ' << D[b]/(double)(N) << ' ' << D[b]/(double)(N*Xbinwidth) << endl; 

    }
  }
  else // if Xmin=Xmax, then this is a delta peak. To show this, we need the following trick.
    {
      fout << Xmin-1e-10 << ' ' << 0 << endl;
      fout << Xmin << ' ' << 1 << endl;
      fout << Xmax+1e-10 << ' ' << 0 << endl;
    }

  fout << endl;
  fout << endl;
 
  fout.close();

}
///////////////////////////////////////////////////////////////////////////






//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// linear bin for continuous distribution, 
// also export the result to two vectors:
// P: each element is a pair value: (w, P(w))
// Q: each element is a pair value: (w, wP(w)/<w>)
void GetPQw(const vector<double>& W, int Nbins, vector<ProbofW>& P, vector<ProbofW>& Q)
{
  int N = W.size();
  //cout << "Size= " << N << endl; //test

  double Wmax, Wmin, Wave, Wvar, Werror;
  Statistics(W, Wmax, Wmin, Wave, Wvar, Werror);

  double Wbinwidth  = (Wmax-Wmin)/Nbins; // h bin's width

  if(Wbinwidth>1e-10) {
    vector<int> D(Nbins,0);                // the bin vector 
    for(int i=0;i<N;i++) {	
      int b = (int)floor((W[i]-Wmin)/Wbinwidth);

      if(b>=Nbins) b = Nbins-1;
      if(b<0)      b = 0;
      D[b]++;    // get the histogram 
    }

    int Dmax = 0;
    int Dsum = 0;
    for(int b=0;b<Nbins;b++){
      Dsum += D[b];
      if(D[b]>Dmax)
	Dmax = D[b];
    }

    for(int b=0; b<Nbins; b++) {
      ProbofW prob;
      prob.w = Wmin+(b+0.5)*Wbinwidth;
      prob.p = D[b]/(double)(N);
      P.push_back(prob);

      prob.p = prob.w * prob.p /Wave;
      Q.push_back(prob);
    }
  }
  else {// if Xmin=Xmax, then this is a delta peak. To show this, we need the following trick.
    
    ProbofW prob;
    //prob.w = Wmin-1e-10;
    //prob.p = 0;
    //P.push_back(prob);
    //Q.push_back(prob);
    
    prob.w = Wmin;
    prob.p = 1;
    P.push_back(prob);
    Q.push_back(prob);

    //prob.w = Wmin+1e-10;
    //prob.p = 0;
    //P.push_back(prob);
    //Q.push_back(prob);
  }
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////






//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// linear bin for continuous distribution, 
// also export the result to one vector:
// P: each element is a pair value: (w, P(w))
void GetPw(const vector<double>& W, int Nbins, vector<ProbofW>& P)
{
  int N = W.size();
  //cout << "Size= " << N << endl; //test

  double Wmax, Wmin, Wave, Wvar, Werror;
  Statistics(W, Wmax, Wmin, Wave, Wvar, Werror);

  // force w in [0,1] (this is used in the edgetyping project, where the inner-module index z is in [0,1])
  Wmin = 0;
  Wmax = 1;

  double Wbinwidth  = (Wmax-Wmin)/Nbins; // h bin's width

  if(Wbinwidth>1e-10) {
    vector<int> D(Nbins,0);                // the bin vector 
    for(int i=0;i<N;i++) {	
      int b = (int)floor((W[i]-Wmin)/Wbinwidth);

      if(b>=Nbins) b = Nbins-1;
      if(b<0)      b = 0;
      D[b]++;    // get the histogram 
    }

    int Dmax = 0;
    int Dsum = 0;
    for(int b=0;b<Nbins;b++){
      Dsum += D[b];
      if(D[b]>Dmax)
	Dmax = D[b];
    }

    for(int b=0; b<Nbins; b++) {
      ProbofW prob;
      prob.w = Wmin + (b+0.5)*Wbinwidth;
      prob.p = D[b]/(double)(N);
      P.push_back(prob);
    }
  }
  else {// if Xmin=Xmax, then this is a delta peak. To show this, we need the following trick.
    
    ProbofW prob;
    //prob.w = Wmin-1e-10;
    //prob.p = 0;
    //P.push_back(prob);
    //Q.push_back(prob);
    
    prob.w = Wmin;
    prob.p = 1;
    P.push_back(prob);

    //prob.w = Wmin+1e-10;
    //prob.p = 0;
    //P.push_back(prob);
    //Q.push_back(prob);
  }
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////








//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// The following is to get the discrete P(X) distribution (especially useful for small networks)
// It is copied from void Network::ClassifyCSL(vector<double>& CSL, 
//			  vector<int>& CSLtype, 
//			  int& NCSLtype, 
//			  vector<double>& CSL_of_type, 
//			  vector<int>& NUM_of_type)

// This template function doesn't work.
// Compiling error:
// Network.cpp:(.text+0x405d): undefined reference to `void GetDiscreteDistribution<int>(std::vector<int, std::allocator<int> >&, char*)'
// why??????????

/*
  template <typename T>
  void GetDiscreteDistribution(vector<T>& X, char* fname)
  {
  vector<int> Xtype;
  int NXtype;
  vector<T> X_of_type;
  vector<int> NUM_of_type;

  int N = X.size();
  vector<T> Xtemp(N);
  copy(X.begin(), X.end(), Xtemp.begin());   // copy X to Xtemp
  sort(Xtemp.begin(), Xtemp.end() );// sort those X values in Xtemp

  // get the number of X types
  X_of_type.push_back(Xtemp[0]);
  for(int i=1; i<N; i++) {
  if(fabs(Xtemp[i]-Xtemp[i-1])>1e-10) 
  X_of_type.push_back(Xtemp[i]);
  }
  NXtype = X_of_type.size();
  NUM_of_type.resize(NXtype, 0);
  
  // get the Xtype for each node 
  Xtype.resize(N);
  for(int i=0; i<N; i++) {
  for(int type=0; type<NXtype; type++) {
  if(fabs(X[i]-X_of_type[type])<=1e-10) { 
  Xtype[i] = type;
  NUM_of_type[type]++;
  break;
  }
  }
  //cout << i << ' ' << Xl[i] << ' ' << Xltype[i] << endl; //debug
  }

  char filename [256];
  sprintf(filename,"%s.discrete", fname);
  ofstream fout(filename, ios_base::out);
  cout.precision(10);
  for(int type=0; type<NXtype; type++) {
  cout << " X["     << type << "] = " << X_of_type[type] 
  << " Num["     << type << "] = " << NUM_of_type[type] 
  << " density[" << type << "] = " << NUM_of_type[type]/(double)N
  << endl; 
  fout << X_of_type[type]  << "    " << NUM_of_type[type]/(double)N << endl;
  }
  fout.close();

  sprintf(filename,"%s.discrete.u", fname);
  ofstream fout2(filename, ios_base::out);
  for(int type=0; type<NXtype; type++) {
  // to show these peaks in a better way, we add some ghost zeros around those peaks
  double eps=1e-10;
  fout2 << X_of_type[type] - eps  << "    " << 0 << endl;
  fout2 << X_of_type[type]  << "    " << NUM_of_type[type]/(double)N << endl;
  fout2 << X_of_type[type] + eps  << "    " << 0 << endl;
  }
  fout2.close();

  }
*/
///////////////////////////////////////////////////////////////////////////



////////////////////////////////////////////////////////////
//case 5: discrete without binning
void GetDiscreteDistribution(vector<double>& X, char* fname)
{
  double Xmax, Xmin, Xave, Xvar, Xerror;
  Statistics(X, Xmax, Xmin, Xave, Xvar, Xerror);


  //////// designed for the calculation of Zscore:= (Xreal- <Xrand>)/ SigmaX      ///////////////////////
  ofstream Xout("X.dat", ios_base::app);
  Xout << Xave << ' ' << sqrt(Xvar) << ' ' << Xerror << endl;
  Xout.close();
  ///////////////////////////////////////////////////////////////////////////////////////////////////////


  vector<int> Xtype;
  int NXtype;
  vector<double> X_of_type;
  vector<int> NUM_of_type;

  int N = X.size();
  vector<double> Xtemp(X.begin(), X.end());
  sort(Xtemp.begin(), Xtemp.end() );// sort those X values in Xtemp

  // get the number of X types
  X_of_type.push_back(Xtemp[0]);
  for(int i=1; i<N; i++) {
    if(fabs(Xtemp[i]-Xtemp[i-1])>1e-10) 
      X_of_type.push_back(Xtemp[i]);
  }
  NXtype = X_of_type.size();
  NUM_of_type.resize(NXtype, 0);
  
  // get the Xtype for each node 
  Xtype.resize(N);
  for(int i=0; i<N; i++) {
    for(int type=0; type<NXtype; type++) {
      if(fabs(X[i]-X_of_type[type])<=1e-10) { 
	Xtype[i] = type;
	NUM_of_type[type]++;
	break;
      }
    }
    //cout << i << ' ' << Xl[i] << ' ' << Xltype[i] << endl; //debug
  }


  // get Pr(X>x) i.e. the Complementary cumulative distribution function 
  vector<double> CCDF1(NXtype); 
  double temp = 0;
  for(int type=0; type<NXtype; type++) {
    temp += NUM_of_type[type];
    CCDF1[type] = N-temp;
  } 
  // Note that if we use this definition, then for the largest Xtype, 
  // CCDF[type] = 0, because Pr(X>xmax)=0
  // Since eventually we want to plot CCDF on log-log plot, y=0 will cause trouble.
  // To solve this issue, we can define CCDF to be Pr(X>=x), as Aaron Clauset did in his paper 
  // "Power-law distribution in emphricial data". If we use this definition we have

  // get Pr(X>=x), another definition of the Complementary cumulative distribution function 
  vector<double> CCDF2(NXtype); 
  CCDF2[0] = N;
  temp = 0;
  for(int type=1; type<NXtype; type++) {
    temp += NUM_of_type[type-1];
    CCDF2[type] = N-temp;
  } 


  char filename [256];
  sprintf(filename,"%s.discrete", fname);
  ofstream fout(filename, ios_base::out);
  for(int type=0; type<NXtype; type++) {
    fout << X_of_type[type]              << "    "  // 1
         << NUM_of_type[type]            << "    "  // 2
         << NUM_of_type[type]/(double)N  << "    "  // 3
         << CCDF1[type]                  << "    "  // 4
         << CCDF1[type]/(double)N        << "    "  // 5
         << CCDF2[type]                  << "    "  // 6
         << CCDF2[type]/(double)N        << endl;   // 7 

  }
  fout.close();


  sprintf(filename,"%s.discrete.u", fname);
  ofstream fout2(filename, ios_base::out);
  //ofstream fout2(filename, ios_base::app);
  fout2 << endl << endl;
  for(int type=0; type<NXtype; type++) {
    // to show these peaks in a better way, we add some ghost zeros around those peaks
    double eps=1e-10;
    fout2 << X_of_type[type] - eps  << "  " << 0 <<  endl;
    fout2 << X_of_type[type]        << "  " << NUM_of_type[type]/(double)N << endl;
    fout2 << X_of_type[type] + eps  << "  " << 0 <<  endl;
  }
  fout2.close();

}
/////////////////////////////////////////////////////////////////////////////////////////////////








/////////////////////////////////////////////////////////////////////////////////////////////////
// comparison, not case sensitive.
bool compare_nocase (string first, string second)
{
  unsigned int i=0;
  while ( (i<first.length()) && (i<second.length()) )
    {
      if (tolower(first[i])<tolower(second[i])) return true;
      else if (tolower(first[i])>tolower(second[i])) return false;
      ++i;
    }
  if (first.length()<second.length()) return true;
  else return false;
}
/////////////////////////////////////////////////////////////////////////////////////////////////




/////////////////////////////////////////////////////////////////////////////////////////////////
void GetDiscreteDistribution(vector<string>& X, char* fname)
{
  int N = X.size();
  cout << "\n There are in total " << N << " items.\n";
  set<string> itemset;
  for(int i=0; i<N; i++) {
    if(X[i]=="  ") X[i]=" Unassigned ";

    if(X[i]!=" subsystem " && X[i]!=" compartment ")
      itemset.insert(X[i]);
  }

  list<string> liststr(itemset.begin(), itemset.end()); 
  liststr.sort(compare_nocase);

  vector<string> basket(liststr.begin(), liststr.end()); 
  int Ntype = basket.size();
  vector<int> count(Ntype,0);

  int Ncount = 0;
  for(int i=0; i<N; i++) {
    for(int j=0; j<Ntype; j++) {
      if(X[i]==basket[j]) { // how can we make sure that "tRNA charging" and "tRNA Charging" are the same 
	count[j]++;
	Ncount++;
      }
    }
  }

  for(int j=0; j<Ntype; j++) {
    if(basket[j]=="  ") {
      cout << j << ',' << basket[j] << endl;
    }
  }
  
  char filename [256];
  sprintf(filename,"%s.discrete", fname);
  ofstream fout(filename, ios_base::out);

  char filename1 [256];
  sprintf(filename1,"%s.discrete.u", fname);
  ofstream fout1(filename1, ios_base::out);

  for(int j=0; j<Ntype; j++) {
    cout  << basket[j] << " ; " << count[j] << " ; " << count[j]/(double)Ncount << endl;
    fout1 << basket[j] << " ; " << count[j] << " ; " << count[j]/(double)Ncount << endl;
    fout << count[j] << "  " << count[j]/(double)Ncount << endl;
  }
  fout.close();
  fout1.close();

  char filename2 [256];
  sprintf(filename2,"%s.discrete.plotscript", fname);
  ofstream gout(filename2, ios_base::out);
  
  gout << "set terminal postscript eps enhanced color solid" << endl;
  gout << "set size 1.0,2.75" << endl;
  gout << "set ytics nomirror" << endl;
  gout << "set boxwidth 0.4" << endl;
  gout << "set yrange [0:*]" << endl;
  gout << "set xlabel ' ' " << endl;
  gout << "set xrange [*:*]" << endl;
  gout << "set xtics nomirror rotate by -90 ( \\" << endl;
  for(int j=0; j<Ntype-1; j++) {
    gout << "\"" << basket[j] << "\"        " << j << ",\\" << endl;
  }
  gout << "\"" << basket[Ntype-1] << "\"        " << Ntype-1 << ")" << endl;    
  gout << "set ylabel '#' font 'Times, 25' " << endl;
  gout << "set style data histogram" << endl;
  gout << "set style histogram rowstacked" << endl;
  gout << "set style fill solid border -1" << endl;
  gout << "set key horizontal " << endl;
  gout << "set output '" << filename << ".eps'" << endl;
  gout << "plot \\" << endl;
  gout << "'"<< filename << "' u 2 noti" << endl;
  gout.close();


}
/////////////////////////////////////////////////////////////////////////////////////////////////








/////////////////////////////////////////////////////////////////////////////////////////////////
void GetPz(vector<int>& K, vector<double>& P, double& z)
{
  int N = K.size();

  //int kmin = 0; // by default
  int kmax = 0;
  double kave = 0;
  for(int i=0; i<N; i++) {
    if(K[i]>kmax)
      kmax = K[i];

    kave += K[i];
  }
  kave /= (double)N;
  
  z = kave;

  P.clear();
  P.resize(kmax+1,0);
  for(int k=0; k<=kmax; k++) {
    P[k] = count(K.begin(), K.end(), k)/(double)N;
    //cout << "P[k=" << k << "]=" << P[k] << endl;
  }

}
/////////////////////////////////////////////////////////////////////////////////////////////////




/////////////////////////////////////////////////////////////////////////////////////////////////
// P[k] : the degree distribution
// Q[k] : the remaining degree distribution
// z: the mean degree 
void GetPQz(vector<int>& K, vector<double>& P, vector<double>& Q, double& z)
{
  int N = K.size();

  //int kmin = 0; // by default
  int kmax = 0;
  double kave = 0;
  for(int i=0; i<N; i++) {
    if(K[i]>kmax)
      kmax = K[i];

    kave += K[i];
  }
  kave /= (double)N;
  
  z = kave;

  P.clear();
  P.resize(kmax+1,0);
  for(int k=0; k<=kmax; k++) {
    P[k] = count(K.begin(), K.end(), k)/(double)N;
    //cout << "P[k=" << k << "]=" << P[k] << endl;
  }

  Q.clear();
  //Q.resize(kmax,0);
  //for(int k=0; k<kmax; k++) {
  //Q[k] = (k+1)*P[k+1]/kave;  
  //}
 
  Q.resize(kmax+1,0);
  for(int k=0; k<=kmax; k++) {
    Q[k] = k*P[k]/kave;  
  }
}
/////////////////////////////////////////////////////////////////////////////////////////////////




/////////////////////////////////////////////////////////////////////////////////////////////////
// Count[k] : the degree histogram
// P[k] : the degree distribution
// Q[k] : the remaining degree distribution
// z: the mean degree 
// v: the variance
void GetPQz(vector<int>& K, vector<int>& Count, vector<double>& P, vector<double>& Q, double& z, double& v)
{
  int N = K.size();

  //int kmin = 0; // by default
  int kmax = 0;
  double kave = 0;
  for(int i=0; i<N; i++) {
    if(K[i]>kmax)
      kmax = K[i];

    kave += K[i];
  }
  kave /= (double)N;
  
  z = kave;



  ////////// variance 
  double ep = 0.0;
  double kvar = 0.0;
  double s = 0.0;
  for (int i=0; i<N; i++) {
    ep += (s=K[i]-kave);
    kvar += s*s;
  }
  kvar=(kvar-ep*ep/N)/(N-1); // Corrected two-pass algorithm
  
  v= kvar;


  //cout << z << ' ' << v << endl; //test

  Count.clear();
  Count.resize(kmax+1,0);

  P.clear();
  P.resize(kmax+1,0);

  for(int k=0; k<=kmax; k++) {
    Count[k] = count(K.begin(), K.end(), k);
    P[k] = Count[k]/(double)N;
    //P[k] = count(K.begin(), K.end(), k)/(double)N;
    //cout << "P[k=" << k << "]=" << P[k] << endl;
  }

  Q.clear();
  //Q.resize(kmax,0);
  //for(int k=0; k<kmax; k++) {
  //Q[k] = (k+1)*P[k+1]/kave;  
  //}

  Q.resize(kmax+1,0);
  //double sum = 0;
  for(int k=0; k<=kmax; k++) {
    Q[k] = k*P[k]/kave;  
    //sum += Q[k];
  }
  //cout << "sumQ[k] = " << sum << endl; //test
 
}
/////////////////////////////////////////////////////////////////////////////////////////////////




/////////////////////////////////////////////////////////////////////////////////////////////////
// P[k] : the degree distribution
// Q[k] : the remaining degree distribution
// z: the mean degree 
// v: the variance
void GetPQz(vector<int>& K, vector<double>& P, vector<double>& Q, double& z, double& v)
{
  int N = K.size();

  //int kmin = 0; // by default
  int kmax = 0;
  double kave = 0;
  for(int i=0; i<N; i++) {
    if(K[i]>kmax)
      kmax = K[i];

    kave += K[i];
  }
  kave /= (double)N;
  
  z = kave;



  ////////// variance 
  double ep = 0.0;
  double kvar = 0.0;
  double s = 0.0;
  for (int i=0; i<N; i++) {
    ep += (s=K[i]-kave);
    kvar += s*s;
  }
  kvar=(kvar-ep*ep/N)/(N-1); // Corrected two-pass algorithm
  
  v= kvar;


  //cout << z << ' ' << v << endl; //test

  P.clear();
  P.resize(kmax+1,0);
  for(int k=0; k<=kmax; k++) {
    P[k] = count(K.begin(), K.end(), k)/(double)N;
    //cout << "P[k=" << k << "]=" << P[k] << endl;
  }

  Q.clear();
  //Q.resize(kmax,0);
  //for(int k=0; k<kmax; k++) {
  //Q[k] = (k+1)*P[k+1]/kave;  
  //}

  Q.resize(kmax+1,0);
  //double sum = 0;
  for(int k=0; k<=kmax; k++) {
    Q[k] = k*P[k]/kave;  
    //sum += Q[k];
  }
  //cout << "sumQ[k] = " << sum << endl; //test
 
}
/////////////////////////////////////////////////////////////////////////////////////////////////


/////////////////////////////////////////////////////////////////////////////////////////////////
// Count[k] : the degree histogram
// P[k] : the degree distribution
// Q[k] : the remaining degree distribution
// z: the mean degree 
// v: the variance
void GetPQz(vector<int>& K, vector<int>& Count, vector<long double>& P, vector<long double>& Q, double& z, double& v)
{
  int N = K.size();

  //int kmin = 0; // by default
  int kmax = 0;
  double kave = 0;
  for(int i=0; i<N; i++) {
    if(K[i]>kmax)
      kmax = K[i];

    kave += K[i];
  }
  kave /= (double)N;
  
  z = kave;



  ////////// variance 
  double ep = 0.0;
  double kvar = 0.0;
  double s = 0.0;
  for (int i=0; i<N; i++) {
    ep += (s=K[i]-kave);
    kvar += s*s;
  }
  kvar=(kvar-ep*ep/N)/(N-1); // Corrected two-pass algorithm
  
  v= kvar;


  //cout << z << ' ' << v << endl; //test

  Count.clear();
  Count.resize(kmax+1,0);

  P.clear();
  P.resize(kmax+1,0);

  for(int k=0; k<=kmax; k++) {
    Count[k] = count(K.begin(), K.end(), k);
    P[k] = Count[k]/(double)N;
    //P[k] = count(K.begin(), K.end(), k)/(double)N;
    //cout << "P[k=" << k << "]=" << P[k] << endl;
  }

  Q.clear();
  //Q.resize(kmax,0);
  //for(int k=0; k<kmax; k++) {
  //Q[k] = (k+1)*P[k+1]/kave;  
  //}

  Q.resize(kmax+1,0);
  //double sum = 0;
  for(int k=0; k<=kmax; k++) {
    Q[k] = k*P[k]/kave;  
    //sum += Q[k];
  }
  //cout << "sumQ[k] = " << sum << endl; //test
 
}
/////////////////////////////////////////////////////////////////////////////////////////////////




/////////////////////////////////////////////////////////////////////////////////////////////////
// P[k] : the degree distribution
// Q[k] : the remaining degree distribution
// z: the mean degree 
// v: the variance
void GetPQz(vector<int>& K, vector<long double>& P, vector<long double>& Q, double& z, double& v)
{
  int N = K.size();

  //int kmin = 0; // by default
  int kmax = 0;
  double kave = 0;
  for(int i=0; i<N; i++) {
    if(K[i]>kmax)
      kmax = K[i];

    kave += K[i];
  }
  kave /= (double)N;
  
  z = kave;



  ////////// variance 
  double ep = 0.0;
  double kvar = 0.0;
  double s = 0.0;
  for (int i=0; i<N; i++) {
    ep += (s=K[i]-kave);
    kvar += s*s;
  }
  kvar=(kvar-ep*ep/N)/(N-1); // Corrected two-pass algorithm
  
  v= kvar;


  //cout << z << ' ' << v << endl; //test

  P.clear();
  P.resize(kmax+1,0);
  for(int k=0; k<=kmax; k++) {
    P[k] = count(K.begin(), K.end(), k)/(double)N;
    //cout << "P[k=" << k << "]=" << P[k] << endl;
  }

  Q.clear();
  //Q.resize(kmax,0);
  //for(int k=0; k<kmax; k++) {
  //Q[k] = (k+1)*P[k+1]/kave;  
  //}

  Q.resize(kmax+1,0);
  //double sum = 0;
  for(int k=0; k<=kmax; k++) {
    Q[k] = k*P[k]/kave;  
    //sum += Q[k];
  }
  //cout << "sumQ[k] = " << sum << endl; //test
 
}
/////////////////////////////////////////////////////////////////////////////////////////////////



/////////////////////////////////////////////////////////////////////////////////////////////////
/*
// P[k] : the degree distribution
// Q[k] : the remaining degree distribution
// z: the mean degree 
void GetPQz(vector<int>& Kin, vector<int>& Kout, vector<double>& P, vector<double>& Q, double& z)
{
Rand rand;
rand.seed(1);

int N = Kin.size();
vector<int> K(N);

int kmin = 0; // by default
int kmax = 0;
double kave = 0;
for(int i=0; i<N; i++) {
int a =  rand.discrete(0,N-1);
int b =  rand.discrete(0,N-1);
int c =  rand.discrete(0,N-1);
int d =  rand.discrete(0,N-1);
K[i] = Kin[a] + Kout[b] + Kin[c] + Kout[d]; // what does this mean??????
//cout << a << ' ' << b << ' ' << c << ' ' << d << ' ' << K[i] << endl; 
if(K[i]>kmax)
kmax = K[i];

kave += K[i];
}
kave /= (double)N;
  
z = kave;

P.clear();
P.resize(kmax+1,0);
for(int k=0; k<=kmax; k++) {
P[k] = count(K.begin(), K.end(), k)/(double)N;
//cout << "P[k=" << k << "]=" << P[k] << endl;
}

Q.clear();
//Q.resize(kmax,0);
//for(int k=0; k<kmax; k++) {
//Q[k] = (k+1)*P[k+1]/kave;  
//}
Q.resize(kmax+1,0);
for(int k=0; k<=kmax; k++) {
Q[k] = k*P[k]/kave;  
}
  
}
*/
 /////////////////////////////////////////////////////////////////////////////////////////////////





