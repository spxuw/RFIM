#include "Statistics.h"


double Beta(double t, double x, double Itot) 
{
  double B = (1e-3) * pow(t,x-1)*pow(1-t,Itot-x-1)/gsl_sf_beta(x,Itot-x);
  if(B<1e-16) B = 1e-16;

  return B;
}


double Gaussian(double x, double mean, double var) 
{
  //return 1.0/sqrt(2*PI*var) * exp(-(x-mean)*(x-mean)/(2*var));
  //return (1e-3) * 1.0/sqrt(2*PI*var) * exp(-(x-mean)*(x-mean)/(2*var));
  //this gives the probability that x is in the range that [x, x+dx], with dx = 1e-3.

  double G = (1e-3) * 1.0/sqrt(2*PI*var) * exp(-(x-mean)*(x-mean)/(2*var));
  if(G<1e-16)  G = 1e-16;
  
  return G;
  
}


long double Poisson(int lambda, int k)
{
  int sign = (k>=lambda) ? 1 : -1;

  int x = sign*(k-lambda);
  int n = x / XMAX;
  int r = x % XMAX;
  long double P = exp(sign*r);

  for(int i=1;i<=k;i++){
    P *= lambda/(EXP*i);
  }

  if(sign==1)
    for(int i=1;i<=n;i++)
      P *= EXPXMAX;
  else
    for(int i=1;i<=n;i++)
      P /= EXPXMAX;
    
  return P;
}


/*inline double Poisson(int lambda, int x)
  {
  double P = 1.0/exp((double)lambda);
  for(int i=1;i<=x;i++)
  P *= (lambda/(double)i);

  return P;
  }
*/



/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Statistics(const vector<double>& data, double& ave, double& squareave, double& var)
{
  int k;
  int n = data.size();

  ////////// average value
  double s = 0.0;
  double s2 = 0.0;
  for (k=0;k<n;k++)    {
    s +=  data[k];
    s2 += data[k]*data[k];
  }
  ave = s/n;
  squareave = s2/n;

  ////////// variance 
  double ep = 0.0;
  var = 0.0;
  for (k=0;k<n;k++) {
    ep += (s=data[k]-ave);
    var += s*s;
  }
  var=(var-ep*ep/n)/(n-1); // Corrected two-pass algorithm
}
////////////////////////////////////////////////////////////////////////////////////////////////


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// find the second minimum value 
double GetXsecondmin(const vector<double>& data, double Xmin)
{
  int n = data.size();
  double MIN=+1e300;
  for (int k=0;k<n;k++)    {
    if(data[k]<MIN && data[k]>Xmin)
      MIN = data[k];
  }
  return MIN;
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Statistics(const vector<double>& data, double& max, double& min, double& ave, double& var, double& error)
{
  int k;
  int n = data.size();

  double MAX=-1e300;
  double MIN=+1e300;
  ////////////////////          Simple Statistics         ///////////////////////////////////
  ////////// average value
  double s = 0.0;
  for (k=0;k<n;k++)
    {
      s += data[k];
 
      if(data[k]>MAX)
	MAX = data[k];
      if(data[k]<MIN)
	MIN = data[k];
    }
  ave = s/n;
  max = MAX;
  min = MIN;

  ////////// variance 
  double ep = 0.0;
  var = 0.0;
  for (k=0;k<n;k++) 
    {
      ep += (s=data[k]-ave);
      var += s*s;
    }
  var=(var-ep*ep/n)/(n-1); // Corrected two-pass algorithm

  //  Finally, the error bar is given by square root of Variance/N_effective
  //  here n_effective = n/Kappa
  double Kappa = 1.0;
  error = sqrt(var*Kappa/n);

    
    cout <<"\n size            = " << n;
    cout <<"\n max             = " << max;
    cout <<"\n min             = " << min;
    cout <<"\n mean            = " << ave;
    cout <<"\n variance        = " << var;
    cout <<"\n sigma           = " << sqrt(var);
    cout <<"\n error of mean   = " << error << endl;
  

}
////////////////////////////////////////////////////////////////////////////////////////////////






/*
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Statistics(const vector<double>& data, double& max, double& min, double& ave, double& var, double& error, double& ci)
{
  int k;
  int n = data.size();

  double MAX=-1e300;
  double MIN=+1e300;
  ////////////////////          Simple Statistics         ///////////////////////////////////
  ////////// average value
  double s = 0.0;
  for (k=0;k<n;k++)
    {
      s += data[k];
 
      if(data[k]>MAX)
	MAX = data[k];
      if(data[k]<MIN)
	MIN = data[k];
    }
  ave = s/n;
  max = MAX;
  min = MIN;

  ////////// variance 
  double ep = 0.0;
  var = 0.0;
  for (k=0;k<n;k++) 
    {
      ep += (s=data[k]-ave);
      var += s*s;
    }
  var=(var-ep*ep/n)/(n-1); // Corrected two-pass algorithm

  //  Finally, the error bar is given by square root of Variance/N_effective
  //  here n_effective = n/Kappa
  double Kappa = 1.0;
  error = sqrt(var*Kappa/n);

  if(n<=3)
    ci = error*4.0;
  else
    ci = error*2.0;
  
  
    cout <<"\n size            = " << n;
    cout <<"\n max             = " << max;
    cout <<"\n min             = " << min;
    cout <<"\n mean            = " << ave;
    cout <<"\n variance        = " << var;
    cout <<"\n sigma           = " << sqrt(var);
    cout <<"\n error of mean   = " << error; 
    cout <<"\n confidence interval (95%) = " << ci;
	 << endl;
  

}
////////////////////////////////////////////////////////////////////////////////////////////////
*/








/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Statistics(const vector<int>& data, double& max, double& min, double& ave, double& var, double& error)
{
  int k;
  int n = data.size();

  double MAX=-1e300;
  double MIN=+1e300;
  ////////////////////          Simple Statistics         ///////////////////////////////////
  ////////// average value
  double s = 0.0;
  for (k=0;k<n;k++)
    {
      s += data[k];

      if(data[k]>MAX)
	MAX = data[k];
      if(data[k]<MIN)
	MIN = data[k];
    }
  ave = s/n;
  max = MAX;
  min = MIN;

  ////////// variance 
  double ep = 0.0;
  var = 0.0;
  for (k=0;k<n;k++) 
    {
      ep += (s=data[k]-ave);
      var += s*s;
    }
  var=(var-ep*ep/n)/(n-1); // Corrected two-pass algorithm

  //  Finally, the error bar is given by square root of Variance/N_effective
  //  here n_effective = n/Kappa
  double Kappa = 1.0;
  error = sqrt(var*Kappa/n);

  
  cout <<"\n size            = " << n;
  cout <<"\n [min, mean, max]= [" << min << " , " <<  ave << " , " << max << ']';
  cout <<"\n variance        = " << var;
  cout <<"\n sigma           = " << sqrt(var);
  cout <<"\n error of mean   = " << error << endl;
  

  cout << n << "   " << min << "  " <<  ave << "  " << max << "  " << sqrt(var) << "  " << error << endl;

}
////////////////////////////////////////////////////////////////////////////////////////////////





/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Statistics(const vector<int>& data, double& xmax, double& xmin, double& xave, double& xvar, double& xerror, double& x2ave)
{
  int k;
  int n = data.size();

  double MAX=-1e300;
  double MIN=+1e300;
  ////////////////////          Simple Statistics         ///////////////////////////////////
  ////////// average value
  double s = 0.0;
  double s2 = 0.0;
  for (k=0;k<n;k++)
    {
      s += data[k];
      s2 += data[k]*data[k];

      if(data[k]>MAX)
	MAX = data[k];
      if(data[k]<MIN)
	MIN = data[k];
    }
  xave  = s/n;
  x2ave = s2/n;
  xmax  = MAX;
  xmin  = MIN;

  ////////// variance 
  double ep = 0.0;
  xvar = 0.0;
  for (k=0;k<n;k++) 
    {
      ep += (s=data[k]-xave);
      xvar += s*s;
    }
  xvar=(xvar-ep*ep/n)/(n-1); // Corrected two-pass algorithm

  //  Finally, the error bar is given by square root of Variance/N_effective
  //  here n_effective = n/Kappa
  double Kappa = 1.0;
  xerror = sqrt(xvar*Kappa/n);

    
  cout <<"\n size            = " << n;
  cout <<"\n xmax            = " << xmax;
  cout <<"\n xmin            = " << xmin;
  cout <<"\n xmean           = " << xave;
  cout <<"\n x2mean          = " << x2ave;
  cout <<"\n variance        = " << xvar;
  cout <<"\n sigma           = " << sqrt(xvar);
  cout <<"\n error of mean   = " << xerror << endl;
    
}
////////////////////////////////////////////////////////////////////////////////////////////////





/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Statistics(const vector<int>& data, double& Xmin, double& Xmax, 
		double& X1inverseave, //  < 1/(X+1)>
		double& logX1ave,     //  < log(X+1) >
		double& Xsrave,       //  < X^0.5 > 
		double& Xave,         //  < X > 
		double& X2ave )       //  < x^2 > 
{
  int k;
  int n = data.size();

  double MAX=-1e300;
  double MIN=+1e300;
  ////////////////////          Simple Statistics         ///////////////////////////////////
  ////////// average value
  double s1 = 0.0;
  double s2 = 0.0;
  double s3 = 0.0; 
  double s4 = 0.0;
  double s5 = 0.0; 

  for (k=0;k<n;k++) {
    s1 += 1.0/(data[k]+1.);
    s2 += log(data[k]+1.);
    s3 += sqrt(data[k]);
    s4 += data[k];
    s5 += data[k]*data[k];
    
    if(data[k]>MAX)
      MAX = data[k];
    if(data[k]<MIN)
      MIN = data[k];
  }
  X1inverseave = s1/n;
  logX1ave     = s2/n;
  Xsrave       = s3/n;
  Xave         = s4/n;
  X2ave        = s5/n;

  Xmax  = MAX;
  Xmin  = MIN;

  ////////// variance 
  double s = 0.0;
  double ep = 0.0;
  double xvar = 0.0;
  for (k=0;k<n;k++) 
    {
      ep += (s=data[k]-Xave);
      xvar += s*s;
    }
  xvar=(xvar-ep*ep/n)/(n-1); // Corrected two-pass algorithm

  //  Finally, the error bar is given by square root of Variance/N_effective
  //  here n_effective = n/Kappa
  double Kappa = 1.0;
  double xerror = sqrt(xvar*Kappa/n);

  
  cout <<"\n size            = " << n;
  cout <<"\n xmin            = " << Xmin;
  cout <<"\n xmax            = " << Xmax;
  cout <<"\n <1/(X+1)>       = " << X1inverseave;
  cout <<"\n <log X>         = " << logX1ave;
  cout <<"\n <X^0.5>         = " << Xsrave;
  cout <<"\n <X>             = " << Xave;
  cout <<"\n <X^2>           = " << X2ave;
  cout <<"\n variance        = " << xvar << endl;
  cout <<"\n sigma           = " << sqrt(xvar);
  cout <<"\n error of mean   = " << xerror << endl;
  
}
////////////////////////////////////////////////////////////////////////////////////////////////






// rewritten from NR
void moment(const vector<int>& data, double &ave, double &adev, double &sdev, double &var, double &skew, double &curt)
{
  int j;
  double ep=0.0,s,p;

  int n=data.size();
  if (n <= 1) {cout << "n must be at least 2 in moment"; exit(0);}

  s=0.0;
  for (j=0;j<n;j++) 
    s += data[j];
  ave=s/n;

  adev=var=skew=curt=0.0;
  for (j=0;j<n;j++) {
    adev += fabs(s=data[j]-ave);
    ep += s;
    var += (p=s*s);
    skew += (p *= s);
    curt += (p *= s);
  }
  adev /= n;
  var=(var-ep*ep/n)/(n-1);
  sdev=sqrt(var);
  
  if (var != 0.0) {
    skew /= (n*var*sdev);
    curt=curt/(n*var*var)-3.0;
  } 
  else {
    cout << "No skew/kurtosis when variance = 0 (in moment)"; 
    skew=-100000;
    curt=-100000;
  }
  
}




int InverseCDF(vector<double>& P, double x)
{
  int k;
  double sum = 0;
  for(k=0; sum <= x; k++)  {
    sum += P[k];
  }
  if(sum>x){
    sum -= P[k];
    k--;
  }
  
  return k;
}


///////////////////////////////////////////////////////////////////////////////////////////////////////
void StatisticalDispersion(vector<int>& X, double& H7, double& H8, double& H9, double& H10)
{
  int N = X.size();
  vector<int> Xtemp(X.begin(), X.end());
  sort(Xtemp.begin(), Xtemp.end() );// sort those X values in Xtemp
  int Xmin = Xtemp[0];
  int Xmax = Xtemp[N-1];
  cout << Xmin << ' ' << Xmax << endl; // test

  // get the number of X types
  vector<int> X_of_type;
  X_of_type.push_back(Xtemp[0]);
  for(int i=1; i<N; i++) {
    if(Xtemp[i] != Xtemp[i-1]) 
      X_of_type.push_back(Xtemp[i]);
  }
  int NXtype = X_of_type.size();
  vector<int> NUM_of_type(NXtype, 0);
  
  // get the Xtype for each node 
  vector<int> Xtype(N);
  for(int i=0; i<N; i++) {
    for(int type=0; type<NXtype; type++) {
      if( X[i] == X_of_type[type]) { 
	Xtype[i] = type;
	NUM_of_type[type]++;
	break;
      }
    }
    //cout << i << ' ' << Xl[i] << ' ' << Xltype[i] << endl; //debug
  }

  //for(int type=0; type<NXtype; type++) {
  //cout << X_of_type[type]  << "    " << NUM_of_type[type] << "    " << NUM_of_type[type]/(double)N << endl;
  //}


  ///////////    Now we get the discrete distribution P[k]   ////////////////////////////////////////////////////
  vector<double> P(Xmax+1,0);
  double Kmean = 0;
  for(int k=0; k<=Xmax; k++) {
    for(int type=0; type<NXtype; type++) {
      if(k == X_of_type[type]) 
	P[k] = NUM_of_type[type]/(double)N;
    }
    //cout << k << ' ' << P[k] << endl;
    Kmean += k*P[k];
  }
  //cout << "Xmean = " << Kmean << endl;
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
  


  // Now we calculate the statistical dispersion measure
  // H7 = <|K/<K> - 1|> = 2 \sum_{k=0}^z (1-k/z) P(k) 
  double sum = 0.0;
  for(int k=0; k<=Kmean; k++)  {
    sum += (1. - k/Kmean) * P[k];
  }
  H7 = 2 * sum;


  
  // Shannon entropy 
  // H8 = -<logP(k)> = - sum_k  P(k) logP(k)
  sum = 0.0;
  for(int k=0; k<=Xmax ; k++)  {
    if(P[k]>0) {
      sum -= P[k] * log(P[k]);
    }
  }
  H8 = sum;



  // H9 = RMD : Relative mean difference = MD / <k>
  //    MD = 2 sum_i sum_{j<i} (i-j) P[i] P[j] 
  sum = 0.0;
  for(int i=0; i<=Xmax ; i++)  {
    double term = 0;
    for(int j=0; j<i ; j++)  {
      term += (i-j)*P[j];
    }
    sum += term * P[i];
  }
  H9 = (2 * sum)/Kmean;
  


  //H10 = QCD : Quartile coefficient of dispersion
  //    = (Q3-Q1)/(Q3+Q1)
  double Q1 = InverseCDF(P, 0.25);
  double Q3 = InverseCDF(P, 0.75);
  H10 = (Q3-Q1)/(Q3+Q1);

  
}
/////////////////////////////////////////////////////////////////////////////////////////////////



/////////////////////////////////////////////////////////////////////////////////////////////////////////////
// calculate the value v_P, of the P-th percentile of an ascending ordered dataset containing N elements with
// values v1 <= v2 <= ... <= vN
double Percentile(const vector<double>& v, int P){
  int N = v.size();

  // calculate rank
  double n = (double)P*(N-1)/100. + 1; // this is used by Microsoft Excel
  double k, d; // n = k + d = integer component + decimal component
  d = modf(n, &k);
  
  if((int)k==1)
    return v[0];
  else if((int)k==N)
    return v[N-1];
  else {
    return v[(int)k-1] + d*(v[(int)k]-v[(int)k-1]);
  }
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////



/////////////////////////////////////////////////////////////////////////////////////////////////////////////
// bin data into B bins with (almost) equal # of data points
// If we want B=2 bins, the binning follows:  [min, median); [median,max]
// then the bin size= median-min; max-median
// If we want B=4 bins, the binning follows:  [min, Q1), [Q1, Q2=median), [Q2,Q3), [Q3,max]
// then the bin size= Q1-min; Q2-Q1; Q3-Q2; max-Q4
void Binning(const vector<double>& data, int B, vector<double>& Binsize, vector<double>& Q)
{
  Binsize.clear();
  Binsize.resize(B,0);
  Q.clear();
  Q.resize(B+1,0);

  int N = data.size();
  cout <<"\n size = " << N;
  vector<double> temp(data.begin(), data.end());
  sort(temp.begin(), temp.end());
  
  //get percentiles
  //Q0 = Xtemp[0]; // min
  //Q1 = Percentile(Xtemp, 25);
  //Q2 = Percentile(Xtemp, 50);
  //Q3 = Percentile(Xtemp, 75);
  //Q4 = Xtemp[N-1]; // max
  
  Q[0] = temp[0];   // min
  cout <<"\n Q[0] = " << Q[0];
  for(int b=1; b<B; b++) {
    Q[b] = Percentile(temp, b*100/B);
    cout <<"\n Q[" << b << "] = " << Q[b] << "  " << b*100/B;
    Binsize[b-1] = Q[b]-Q[b-1];
  }
  Q[B] = temp[N-1]; // max
  cout <<"\n Q[" << B << "] = " << Q[B] << endl;
  Binsize[B-1] = Q[B]-Q[B-1];

}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////



/////////////////////////////////////////////////////////////////////////////////////////////////////////////
int GetBinIndex(const double& x, const vector<double>& Q)
{
  int B = Q.size()-1;
  if(x<Q[0] || x>Q[B]) {
    cout << "error!\n";
    exit(0);
  }

  int index = 0;
  
  // in case Q[0]=Q[1], we consider all points with x=Q[0]=Q[1] belong to the 1st bin
  // only for points Q[0]=Q[1] < x < Q[2] belong to the 2nd bin
  if(fabs(Q[0]-Q[1])<1e-6) { 
    if(fabs(x-Q[0])<1e-6)
      index = 0;
    else if(x>Q[0] && x<=Q[2])
      index = 1;
    else {
      for(int b=2; b<B; b++) {
	if(x>Q[b] && x<=Q[b+1])
	  index = b;
      }
    }
  }
  else {

    for(int b=0; b<B; b++) {
      if(x>Q[b] && x<=Q[b+1])
	index = b;
    }

    if(fabs(x-Q[0])<1e-6) // if x=Q[0]
      index = 0;
  }
  
  return index;
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////



/////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Statistics(const vector<double>& data, 
		double& min, double& Q1, double& Q2, double& Q3, double& max, 
		double& ave, double& var, double& error,
		ofstream& fout1, ofstream& fout2)
{
  int n = data.size();

  vector<double> Xtemp(data.begin(), data.end());
  sort(Xtemp.begin(), Xtemp.end());
  
  // get percentiles
  min = Xtemp[0];
  Q1 = Percentile(Xtemp, 25);
  Q2 = Percentile(Xtemp, 50);
  Q3 = Percentile(Xtemp, 75);
  max = Xtemp[n-1];
  fout1 <<"\n size            = " << n;
  fout1 <<"\n min             = " << min;
  fout1 <<"\n Q1              = " << Q1;
  fout1 <<"\n Q2              = " << Q2;
  fout1 <<"\n Q3              = " << Q3;
  fout1 <<"\n max             = " << max;


  // calcluate the cumulative distribution function (CDF)_X (x) = P(X<=x) 
  /*// Method-I
  // a. get the number of X types, i.e. the unique # of data values
  vector<double> X_of_type;
  X_of_type.push_back(Xtemp[0]);
  for(int i=1; i<n; i++) {
    if(fabs(Xtemp[i]-Xtemp[i-1])>1e-10) 
      X_of_type.push_back(Xtemp[i]);
  }
  int NXtype = X_of_type.size();
  vector<int> NUM_of_type(NXtype, 0);
  
  // b. get the Xtype for each data point
  vector<int> Xtype(n);
  for(int i=0; i<n; i++) {
    for(int type=0; type<NXtype; type++) {
      if(fabs(Xtemp[i]-X_of_type[type])<=1e-10) { 
	Xtype[i] = type;
	NUM_of_type[type]++;
	break;
      }
    }
  }

  // c. get CDF(x) = Pr(X<=x) i.e. the cumulative distribution function 
  vector<double> CDF(NXtype); 
  double temp = 0;
  for(int type=0; type<NXtype; type++) {
    temp += NUM_of_type[type];
    CDF[type] = temp;
  } 
  */

  
  //Method-II (equivlant to Method-I, but more efficient)
  // a. get the number of X types, i.e. the unique # of data values
  vector<double> X_of_type;
  X_of_type.push_back(Xtemp[0]);
  vector<int> NUM_of_type;
  NUM_of_type.push_back(1);
  vector<int> Xtype(n);
  Xtype[0] = 0;

  for(int i=1; i<n; i++) {
    if(fabs(Xtemp[i]-Xtemp[i-1])>1e-10) { // find a different X value
      X_of_type.push_back(Xtemp[i]);

      int Ntype = NUM_of_type.size();
      NUM_of_type.push_back(1);
      Xtype[i] = Ntype-1;
    }
    else {
      int Ntype = NUM_of_type.size();
      NUM_of_type[Ntype-1]++;
      Xtype[i] = Ntype;
    }
  }
  int NXtype = X_of_type.size(); // # of X types, i.e. the # of unique X values

  // b. get CDF(x) = Pr(X<=x) i.e. the cumulative distribution function 
  vector<double> CDF(NXtype); 
  double temp = 0;
  for(int type=0; type<NXtype; type++) {
    temp += NUM_of_type[type];
    CDF[type] = temp;
  } 


  for(int type=0; type<NXtype; type++) {
    fout2 << X_of_type[type]              << "    "  // 1
          << CDF[type]/(double)n          << "    "  // 2
          << NUM_of_type[type]/(double)n  << endl;   // 3

  }
  fout2.close();


  ////////////////////          Simple Statistics         ///////////////////////////////////
  ////////// average value
  double s = 0.0;
  for(int k=0;k<n;k++) {
    s += data[k];
  }
  ave = s/n;

  ////////// variance 
  double ep = 0.0;
  var = 0.0;
  for(int k=0;k<n;k++) {
    ep += (s=data[k]-ave);
    var += s*s;
  }
  var=(var-ep*ep/n)/(n-1); // Corrected two-pass algorithm

  //  Finally, the error bar is given by square root of Variance/N_effective
  //  here n_effective = n/Kappa
  double Kappa = 1.0;
  error = sqrt(var*Kappa/n);


  cout <<"\n size            = " << n;
  cout <<"\n min             = " << min;
  cout <<"\n Q1              = " << Q1;
  cout <<"\n Q2              = " << Q2;
  cout <<"\n Q3              = " << Q3;
  cout <<"\n max             = " << max;
  cout <<"\n mean            = " << ave;
  cout <<"\n variance        = " << var;
  cout <<"\n sigma           = " << sqrt(var);
  cout <<"\n error of mean   = " << error << endl;

  fout1 <<"\n mean            = " << ave;
  fout1 <<"\n variance        = " << var;
  fout1 <<"\n sigma           = " << sqrt(var);
  fout1 <<"\n error of mean   = " << error << endl;
  fout1.close();

}
////////////////////////////////////////////////////////////////////////////////////////////////





