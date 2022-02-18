#include "State.h"

using namespace std;


//////////////////////////////////////////////////////////////////////////////////////////////////////
void State::SetRandomField(int dim, int length, double disorder, int rngseed)
{
  Rand rand;     	    // set up the random number generator(RNG)    
  rand.seed(rngseed);   

  int num = 1;
  for(int i=0;i<dim;i++) 
    num *= length;     
  int z=2*dim;

  Vec randomField(num);
  double abs_min = 10000;
  double abs_max = -10000;

  double dE  = 10000;
  double min = 1e300;
  double max = -1e300;
  double sum = 0;
  for(int i=0; i<num; i++) 
    {

      switch(::dist)
	{
	case 0: // 0. Bimodal distributionw -R or +R
	  randomField[i] = rand.bimodal(disorder); 
	  break;

	case 1: // 1. Gaussian distributionw
	  randomField[i] = rand.gaussian(disorder); 
	  break;

	case 2: // 2. rectangular distribution [-R,R]
	  randomField[i] = rand.rectangular(disorder); 
	  break;

	case 3: // 3. Lorentzian distribution
	  randomField[i] = rand.lorentzian(disorder); 
	  break;
		
	case 4: // 4. Parabolic distribution
	  randomField[i] = rand.parabolic(disorder); 
	  break;
		
	default: 
	  cout <<"\n Distribution type is wrong!\n";
	}
	    

      if(randomField[i] > max) max = randomField[i];
      if(randomField[i] < min) min = randomField[i];

      if(fabs(randomField[i]) < abs_min) 
	abs_min = fabs(randomField[i]);
      
      if(fabs(randomField[i]) > abs_max) 
	abs_max = fabs(randomField[i]);

     
      for(int nsame=0; nsame<z+1; nsame++)
	{
	  double tmp= fabs(4*(nsame-dim) + 2*randomField[i]);
	  if( tmp < dE)
	    dE = tmp;
	}

      sum += randomField[i];
    }

  
  State::h.resize(num);
  State::h        = randomField; 
  State::abshmin  = abs_min;
  State::abshmax  = abs_max;
 
  State::dEmin    = dE; 
  State::hmin     = min; 
  State::hmax     = max; 
  State::havg     = sum/(double)num;
  State::hsum     = sum;

  cout << "hmin= " << hmin << endl;
  cout << "hmax= " << hmax << endl;

}

int BisectionSearch(vector<double>& xx, const double x)
{
  int ju,jm,jl;
  bool ascnd;

  int n=xx.size();
  jl=-1; // initialize lower 
  ju=n;  // and upper limits
  ascnd=(xx[n-1] >= xx[0]); // true if ascending order of table, false otherwise

  while (ju-jl > 1) {
    jm=(ju+jl) >> 1; // compute a midpoint
    if ((x >= xx[jm]) == ascnd) jl=jm;       // replace either the lower limit
    else        	          ju=jm;       // or the upper limit, as appropriate       
  }

  int j;
  if (x == xx[0]) j=0;
  else if (x == xx[n-1]) j=n-1;
  else j=jl+1;

  return j;
}
////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////
#include <gsl/gsl_sf_erf.h>
double NormalCDF(double z)
{
  return 0.5* (1 + gsl_sf_erf(z/sqrt(2.0)) );
}
////////////////////////////////////////////////////////////////////////



//////////////////////////////////////////////////////////////////////////////////////////////////////

void State::Set_RandomField_RandomBond(int num, double kmean, double gamma, double disorder, int rngseed)
{
  Rand rand;
  rand.seed(rngseed);

  double alpha = 1.0/(gamma-1.0); // gamma = 1 + 1./alpha;
  vector<double> w(num);
  double Sum_w = 0;
  for(int i=0;i<num;i++) {
    w[i]= pow((double)(i+1),-alpha); // note that (i+1)^(-alpha) since we actually index node from 1 to N
    Sum_w += w[i];
  }

  vector<double> Cum_w(num); // cumulant of w[i]
  double temp = 0;
  for(int i=0;i<num;i++){
    temp    += w[i]/Sum_w;
    Cum_w[i] = temp; 
  }

  int numedge = 0.5* kmean * num; 


  vector<int> Degree(num,0);
  vector<Nbl> Adj(num);
  vector<Link> EdgeList;  
  map<string,int>  IndexMAP; // the map between (i,j) and the edge index e
  map<string,double>  WeightMAP; // the map between (i,j) and the edge weight w_ij
   
  int edge = 0;
  double wsum = 0;


  for(edge=1; edge<=numedge; ) {
    int i = BisectionSearch(Cum_w, rand.ran1());
    int j = BisectionSearch(Cum_w, rand.ran1());
  
    if(!find(Adj[i],j) && i!=j ) {
 
      double wij = rand.ran1();
      

      Link link(i, j, edge, wij); 
      EdgeList.push_back(link);
      wsum += wij;

      stringstream sst1; sst1 << i << ">" << j;  
      stringstream sst2; sst2 << j << ">" << i;
      IndexMAP[sst1.str()] = edge;
      IndexMAP[sst2.str()] = edge;
      WeightMAP[sst1.str()] = wij;
      WeightMAP[sst2.str()] = wij;

      edge++;
       
      Adj[i].push_back(j); 
      Adj[j].push_back(i); 
	
      Degree[i]++;
      Degree[j]++;
    }
  }
  cout << "# of edges E= " << numedge << endl; //test

  State::A.resize(num);
  State::K.resize(num);
  State::LINK.resize(numedge);
  State::NodeName.resize(num);

  State::A      = Adj;
  State::K      = Degree;
  State::LINK   = EdgeList;
  State::MAP    = IndexMAP;
  State::Weight = WeightMAP; 
  State::Wsum   = wsum;
  State::c      = kmean;
 
  Vec randomField(num);
  double abs_min = 10000;
  double abs_max = -10000;

  double dE  = 10000;
  double min = 1e300;
  double max = -1e300;
  double sum = 0;
  double tmp = 0;


  for(int i=0; i<num; i++) 
    {
      NodeName[i] = Int2String(i);

      switch(::dist)
	{
	case 0: // 0. Bimodal distributionw -R or +R
	  randomField[i] = rand.bimodal(disorder); 
	  break;

	case 1: // 1. Gaussian distributionw
	  randomField[i] = rand.gaussian(disorder); 
	  break;

	case 2: // 2. rectangular distribution [-R,R]
	  randomField[i] = rand.rectangular(disorder); 
	  break;

	case 3: // 3. Lorentzian distribution
	  randomField[i] = rand.lorentzian(disorder); 
	  break;
		
	case 4: // 4. Parabolic distribution
	  randomField[i] = rand.parabolic(disorder); 
	  break;
		
	default: 
	  cout <<"\n Distribution type is wrong!\n";
	}

      if(randomField[i] > max) max = randomField[i];
      if(randomField[i] < min) min = randomField[i];

      if(fabs(randomField[i]) < abs_min) 
	abs_min = fabs(randomField[i]);

      if(fabs(randomField[i]) > abs_max) 
	abs_max = fabs(randomField[i]);

      for(int nsame=0; nsame<Degree[i]+1; nsame++)
	{
	  tmp= fabs(4*nsame-2*Degree[i] + 2*randomField[i]);
	  if(tmp < dE) dE = tmp;
	  tmp= fabs(4*nsame-2*Degree[i] - 2*randomField[i]);
	  if(tmp < dE) dE = tmp;
	}

      sum += randomField[i];
    }

  State::h.resize(num);
  State::h        = randomField; 
  State::abshmin  = abs_min;
  State::abshmax  = abs_max;



  State::dEmin    = dE; 
  State::hmin     = min; 
  State::hmax     = max; 
  State::havg     = sum/(double)num;
  State::hsum     = sum;

  cout << "hmin = " << min << endl; //test
  cout << "hmax = " << max << endl; //test
  cout << "havg = " << State::havg << endl; //test

  cout << "wsum = " << State::Wsum << endl; //test
  
 
}
///////////////////////////////////////////////////////////////////////////////


void State::Set_RandomPvalue_RandomBond(int num, double kmean, double gamma, int rngseed)
{
  Rand rand;
  rand.seed(rngseed);
  
  // construct the random network using the static model  
  double alpha = 1.0/(gamma-1.0); // gamma = 1 + 1./alpha;
  vector<double> w(num);
  double Sum_w = 0;
  for(int i=0;i<num;i++) {
    w[i]= pow((double)(i+1),-alpha); // note that (i+1)^(-alpha) since we actually index node from 1 to N
    Sum_w += w[i];
  }

  vector<double> Cum_w(num); // cumulant of w[i]
  double temp = 0;
  for(int i=0;i<num;i++){
    temp    += w[i]/Sum_w;
    Cum_w[i] = temp; 
  }

  int numedge = 0.5* kmean * num; 

  vector<int> Degree(num,0);
  vector<Nbl> Adj(num);
  vector<Link> EdgeList;  
  map<string,int>  IndexMAP; // the map between (i,j) and the edge index e
  map<string,double>  WeightMAP; // the map between (i,j) and the edge weight w_ij
   
  int edge = 0;
  double wsum = 0;


  char fname[256];
  sprintf(fname,"./data/N%d-c%e-gamma%.3lf-seed%d.edgelist", num, kmean, gamma, rngseed);
  ofstream fout1(fname, ios_base::out);

  for(edge=1; edge<=numedge; ) {
    int i = BisectionSearch(Cum_w, rand.ran1());
    int j = BisectionSearch(Cum_w, rand.ran1());
  
    if(!find(Adj[i],j) && i!=j ) {
 
      double wij = rand.ran1();

      
      fout1 << i << ' ' << j << ' ' << wij << endl; 

      Link link(i, j, edge, wij); 
      EdgeList.push_back(link);
      wsum += wij;

      // build a mapping between the edge (i,j) and the edge index e
      stringstream sst1; sst1 << i << ">" << j;  
      stringstream sst2; sst2 << j << ">" << i;
      IndexMAP[sst1.str()] = edge;
      IndexMAP[sst2.str()] = edge;
      WeightMAP[sst1.str()] = wij;
      WeightMAP[sst2.str()] = wij;

      edge++;
       
      Adj[i].push_back(j); 
      Adj[j].push_back(i); 
	
      Degree[i]++;
      Degree[j]++;
    }
  }
  fout1.close();
  cout << "# of edges E= " << numedge << endl; //test

  State::A.resize(num);
  State::K.resize(num);
  State::LINK.resize(numedge);
  State::NodeName.resize(num);

  State::A      = Adj;
  State::K      = Degree;
  State::LINK   = EdgeList;
  State::MAP    = IndexMAP;
  State::Weight = WeightMAP; 
  State::Wsum   = wsum;
  State::c      = kmean;
 
  vector<double> P(num); 

  sprintf(fname,"./data/N%d-c%e-gamma%.3lf-seed%d.nodeweight", num, kmean, gamma, rngseed);
  ofstream fout2(fname, ios_base::out);
  for(int i=0; i<num; i++) {
    P[i] = rand.ran1();
    fout2 << i << ' ' << P[i] << endl;
  }
  fout2.close();


  State::PVALUE.resize(num);
  State::PVALUE   = P;
 
}
///////////////////////////////////////////////////////////////////////////////






//////////////////////////////////////////////////////////////////////////////////////////////////////
void CalCoords(int loc, int D, int* stride, int* newCoords)
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



/////////////////////////////////////////////////////////////////////////////////////////////////
// undirected
// construct a lattice from given Dimension, Linear size and boundary conditions; but store the lattice as a undirected network
void State::Set_RandomField_RandomBond_Lattice(int dim, int a, double disorder, int seed, bool PBC)
{
  Rand rand;     	    // set up the random number generator(RNG)    
  rand.seed(seed);   

  int num=1;
  for(int i=0;i<dim;i++) 
    num*=a;
  int z=2*dim;

  Vec randomField(num);
  double abs_min = 10000;
  double abs_max = -10000;

  double dE  = 10000;
  double min = 1e300;
  double max = -1e300;
  double sum = 0;
  for(int i=0; i<num; i++) {
    randomField[i] = rand.gaussian(disorder); 

    if(randomField[i] > max) max = randomField[i];
    if(randomField[i] < min) min = randomField[i];

   if(fabs(randomField[i]) < abs_min) 
      abs_min = fabs(randomField[i]);

    if(fabs(randomField[i]) > abs_max) 
      abs_max = fabs(randomField[i]);


    for(int nsame=0; nsame<z+1; nsame++)
      {
	double tmp= fabs(4*(nsame-dim) + 2*randomField[i]);
	if( tmp < dE)
	  dE = tmp;
      }

    sum += randomField[i];
  }

  State::h.resize(num);
  State::h        = randomField; 
  State::abshmin  = abs_min;
  State::abshmax  = abs_max;
 
  State::dEmin    = dE; 
  State::hmin     = min; 
  State::hmax     = max; 
  State::havg     = sum/(double)num;
  State::hsum     = sum;

  cout << "hmin= " << hmin << endl;
  cout << "hmax= " << hmax << endl;


  vector<int> Degree(num,0);
  vector<Nbl> Adj(num);
  vector<Link> EdgeList;  
  vector<string> namelist;

  map<string,int>  IndexMAP; // the map between (i,j) and the edge index e
  map<string,double>  WeightMAP; // the map between (i,j) and the edge weight w_ij
   
  int numedge = 0;
  double wsum = 0;
 

  int* Stride = new int[dim];
  Stride[dim-1] = 1;
  for(int i=dim-2;i>=0;i--) {
    Stride[i] = Stride[i+1]*a;
  }
  
  int* NeighborLocs = new int[z];
  int* Coords = new int [dim];  

  if(PBC)
    {
      for(int Loci=0;Loci<num;Loci++)
	{

	  CalCoords(Loci, dim, Stride, Coords); 


	  for(int i=0;i<dim;i++) {
	    for(int j=0;j<dim;j++) {
	      if(j==i) {
		if(Coords[j]==a-1)
		  {
		    NeighborLocs[2*j]=Loci-(a-1)*Stride[j];
		    NeighborLocs[2*j+1]=Loci-Stride[j];
		  }
		else if(Coords[j]==0)
		  {
		    NeighborLocs[2*j]=Loci+Stride[j];
		    NeighborLocs[2*j+1]=Loci+(a-1)*Stride[j];
		  }
		else
		  {
		    NeighborLocs[2*j]=Loci+Stride[j];
		    NeighborLocs[2*j+1]=Loci-Stride[j];
		  }
	      }
	    }
	  }
	  
	  for(int k=0;k<z;k++)
	    {
	      int Locj = NeighborLocs[k];
	      //cout << loci << ' ' << locj << endl; //test

	      if(!find(Adj[Loci], Locj)) {
		//double wij = rand.ran1();  // random bond
		double wij = 1.0;            // identical bond
		Link link(Loci,Locj, numedge, wij);
		EdgeList.push_back(link);
		wsum += wij;
	
		// build a mapping between the edge (i,j) and the edge index e
		stringstream sst1; sst1 << Loci << ">" << Locj;  
		stringstream sst2; sst2 << Locj << ">" << Loci;
		IndexMAP[sst1.str()] = numedge;
		IndexMAP[sst2.str()] = numedge;
		WeightMAP[sst1.str()] = wij;
		WeightMAP[sst2.str()] = wij;
		numedge++;

		Adj[Loci].push_back(Locj);
		Adj[Locj].push_back(Loci);

		Degree[Loci]++; 
		Degree[Locj]++;
	      }
	    }
	}
    }
  else
    {
      for(int Loci=0;Loci<num;Loci++)
	{
	  CalCoords(Loci, dim, Stride, Coords); 

	  for(int i=0;i<dim;i++) {
	    for(int j=0;j<dim;j++) {
	      if(j==i) {
		if(Coords[j]==a-1)
		  {
		    NeighborLocs[2*j]=-1;
		    NeighborLocs[2*j+1]=Loci-Stride[j];
		  }
		else if(Coords[j]==0)
		  {
		    NeighborLocs[2*j]=Loci+Stride[j];
		    NeighborLocs[2*j+1]=num+1;
		  }
		else
		  {
		    NeighborLocs[2*j]=Loci+Stride[j];
		    NeighborLocs[2*j+1]=Loci-Stride[j];
		  }
	      }
	    }
	  }

	  for(int k=0;k<z;k++)
	    {
	      int Locj = NeighborLocs[k];
	      //cout << Loci << ' ' << Locj << endl; //test		

	      if((Locj<num && Locj>0))
		{
		  if(!find(Adj[Loci], Locj)) {

		    //double wij = rand.ran1(); 
		    double wij = 1.0 ;//


		    /*
		    double scale = 0.2;
		    if(Loci <0.6*num && Loci>0.4*num && Locj <0.6*num && Locj>0.4*num)
		      wij *= 1;
		    else
		      wij *= scale;
		    */


		    Link link(Loci,Locj, numedge, wij);
		    EdgeList.push_back(link);
		    wsum += wij;
		    
		    // build a mapping between the edge (i,j) and the edge index e
		    stringstream sst1; sst1 << Loci << ">" << Locj;  
		    stringstream sst2; sst2 << Locj << ">" << Loci;
		    IndexMAP[sst1.str()] = numedge;
		    IndexMAP[sst2.str()] = numedge;
		    WeightMAP[sst1.str()] = wij;
		    WeightMAP[sst2.str()] = wij;
		    numedge++;
		    
		    Adj[Loci].push_back(Locj);
		    Adj[Locj].push_back(Loci);
		    
		    Degree[Loci]++; 
		    Degree[Locj]++;
		  }
		}
	    }
	}
      
    }

  delete [] NeighborLocs;
  delete [] Coords;
 
 cout << "# of edge= " << numedge << endl; //test
 
  State::A.resize(num);
  State::K.resize(num);
  State::LINK.resize(numedge);
  State::NodeName.resize(num);
  

  State::A      = Adj;
  State::K      = Degree;
  State::LINK   = EdgeList;
  State::MAP    = IndexMAP;
  State::Weight = WeightMAP; 
  State::Wsum   = wsum;
  State::c      = z; 

  for(int i=0;i<num;i++)
    NodeName[i] = Int2String(i);
}
/////////////////////////////////////////////////////////////////////////////////////////////////





#include "inversenormalCDF.h"
#include "find.h"

/* read the standard network file 
 #N E pmin p1min pmax wmin wmax wave
  N E pmin p1min pmax wmin wmax wave
 #nodeindex nodename nodeweight
  .
  .
  .
 #source target weight
  .
  .
  .

 Here   
  N      # of nodes
  E      # of edges
  pmin   minimum p-value 
  p1min  minimum non-zero p-value 
  pmax   maximum p-value
  wmin   minimum edge weight
  wmax   maximum edge weight
  wave   average edge weight 

  nodeweight = pvalue
  eightweight = confidence level
*/

///////////////////////////////////////////////////////////////////////////////////////////////// 
void State::Read_Fields_Bonds(string file)
{
  ifstream fin(file.c_str(), ios_base::in);
  if(!fin) {cout << "Cannot open " << file  << " for read."; exit(0);}

  int num, numedge;
  double pmin, p1min, pmax, wmin, wmax, wave;
  string s;
  getline(fin, s, '\n');
  fin >> num >> numedge >> pmin >> p1min >> pmax >> wmin >> wmax >> wave;
  getline(fin, s, '\n'); // to finish reading this line ended by '\n'
  cout << "N         E        pmin        p1min       pmax        wmin        wmax        wave " << endl;
  cout << num << ' ' << numedge << ' ' << pmin << ' ' << p1min << ' ' << pmax << ' ' << wmin << ' ' << wmax << ' ' << wave << endl;
 
  int index;
  string name;
  double pvalue;
  vector<string> namelist;
  map<string,int> nodemap;
  vector<double> P;//(num);
  getline(fin, s, '\n');
  //cout << s << endl;
  for(int i=0;i<num;i++) {
    fin >> index >> name >> pvalue;
    getline(fin, s, '\n'); // to finish reading this line ended by '\n'
    //cout << index << ' ' << name << ' ' << pvalue << endl; //test
    nodemap[name] = index;
    P.push_back(pvalue);//P[i] = pvalue;
    namelist.push_back(name);
  }
 
  vector<int> Degree(num,0);
  vector<Nbl> Adj(num);
  vector<Link> EdgeList;  
  map<string,int>  IndexMAP; // the map between (i,j) and the edge index e
  map<string,double>  WeightMAP; // the map between (i,j) and the edge weight w_ij
   
  double wsum = 0;
  int i,j;
  double wij;

  getline(fin, s, '\n');
  //cout << s << endl;
  for(int e=1;e<=numedge;e++) {
    fin >> i >> j >> wij;
    getline(fin, s, '\n'); // to finish reading this line ended by '\n'
    //cout << i << ' ' << j << ' ' << wij << endl; //test

    stringstream sst1; sst1 << i << ">" << j;  
    stringstream sst2; sst2 << j << ">" << i;
    
    wij /= wave;  // normalize the edge weight such that <wij>=1. YYL 04/11/2016


    if(IndexMAP[sst1.str()]==0) { // avoid reading repeated edges 
      Link link(i, j, e, wij); 
      EdgeList.push_back(link);
      wsum += wij;

      // build a mapping between the edge (i,j) and the edge index e
      IndexMAP[sst1.str()]  = e;
      IndexMAP[sst2.str()]  = e;
      WeightMAP[sst1.str()] = wij;
      WeightMAP[sst2.str()] = wij;
       
      Adj[i].push_back(j); 
      Adj[j].push_back(i); 
      
      Degree[i]++;
      Degree[j]++;
    }
    
  }
  
  fin.close();
  numedge = EdgeList.size(); // update the number of edges 
  cout << "\nThe true number of edges (after removing repeated edges) = " << numedge << endl;  //test
  
  State::A.resize(num);
  State::K.resize(num);
  State::LINK.resize(numedge);
  State::NodeName.resize(num);
  State::PVALUE.resize(num);

  State::A      = Adj;
  State::K      = Degree;
  State::NodeName = namelist;
  State::PVALUE   = P;

  State::LINK = EdgeList;
  State::NodeMAP = nodemap;
  State::MAP    = IndexMAP;
  State::Weight = WeightMAP; 
  State::Wsum   = wsum;
  State::c      = 2*numedge/(double)num;

  Vec randomField(num);
  double abs_min = 10000;
  double abs_max = -10000;

  double dE  = 10000;
  double min = 1e300;
  double max = -1e300;
  double sum = 0;
  double tmp = 0;

  double zscore;

  int N0 = 0;   // # of isolated spins 

  for(int i=0;i<num;i++) {
    zscore = (P[i]==0) ? NormalCDFInverse(1-0.01*p1min) : NormalCDFInverse(1-P[i]); // 1-p[i] = NormalCDF(z[i])
    randomField[i] = zscore;

    if(randomField[i] > max) 
      max = randomField[i];
    if(randomField[i] < min) 
      min = randomField[i];

    if(fabs(randomField[i]) < abs_min) 
      abs_min = fabs(randomField[i]);

    if(fabs(randomField[i]) > abs_max) 
      abs_max = fabs(randomField[i]);

    
    for(int nsame=0; nsame<Degree[i]+1; nsame++) {
      tmp= fabs(4*nsame-2*Degree[i] + 2*randomField[i]);
      if(tmp < dE) 
	dE = tmp;
      tmp= fabs(4*nsame-2*Degree[i] - 2*randomField[i]);
      if(tmp < dE) 
	dE = tmp;
    }

    sum += randomField[i];

    if(Degree[i]==0) 
      N0++;
    
    
  }
  // Note that for those isolated nodes, their spin values will be solely determined by their local field 
  cout << "# of isolated spins   = " << N0  << endl; 

  State::h.resize(num);
  State::h        = randomField; 
  State::abshmin  = abs_min;
  State::abshmax  = abs_max;
 
  State::dEmin    = dE; 
  State::hmin     = min; 
  State::hmax     = max; 
  State::havg     = sum/(double)num;
  State::hsum     = sum;

  cout << "hmin = " << min << endl; //test
  cout << "hmax = " << max << endl; //test
  cout << "havg = " << State::havg << endl; //test
  cout << "Finish reading fields and bonds from file!\n\n";
  
  
}
/////////////////////////////////////////////////////////////////////////////////////////////////


///////////////////////////////////////////////////////////////////////////////////////////////// 
void State::Randomize_Fields(int seed)
{
  int num = h.size();
  cout << num << endl; 

  /*
  for(int i=0;i<20;i++) 
    cout << h[i] << ',';
  cout << endl;
  */

  // This way, those genes without assigned p-values will also be shuffled. 
  random_shuffle(h.begin(), h.end());

  /*
  for(int i=0;i<20;i++) 
    cout << h[i] << ',';
  cout << endl;
  */

}
///////////////////////////////////////////////////////////////////////////////////////////////// 



///////////////////////////////////////////////////////////////////////////////////////////////// 
void State::Randomize_Bonds(int seed)
{
  int num = K.size();
  int numedge = LINK.size();
  double kmean = 2.0*numedge/(double)num;

  vector<int> Degree(num,0);
  vector<Nbl> Adj(num);
  vector<Link> EdgeList;  
  map<string,int>  IndexMAP; // the map between (i,j) and the edge index e
  map<string,double>  WeightMAP; // the map between (i,j) and the edge weight w_ij
   
  Rand RNG;
  RNG.seed(seed);

  double wsum = 0;
  int i,j;
  for(int edge=1; edge<=numedge; ) {
    RNG.discrete(0, num-1, i, j); 
    //cout << i << ' ' << j << endl; //test 
    if(!find(Adj[i],j) && i!=j) {
      double wij = RNG.ran1();
     
      Link link(i, j, edge, wij); 
      EdgeList.push_back(link);
      wsum += wij;

      // build a mapping between the edge (i,j) and the edge index e
      stringstream sst1; sst1 << i << ">" << j;  
      stringstream sst2; sst2 << j << ">" << i;
      IndexMAP[sst1.str()] = edge;
      IndexMAP[sst2.str()] = edge;
      WeightMAP[sst1.str()] = wij;
      WeightMAP[sst2.str()] = wij;

      edge++;
       
      Adj[i].push_back(j); 
      Adj[j].push_back(i); 
	
      Degree[i]++;
      Degree[j]++;
    }
  }
  cout << "# of edges E= " << numedge << endl; //test

  State::A.clear();
  State::K.clear();
  State::LINK.clear();
  State::NodeName.clear();

  State::A.resize(num);
  State::K.resize(num);
  State::LINK.resize(numedge);
  State::NodeName.resize(num);

  State::A      = Adj;
  State::K      = Degree;
  State::LINK   = EdgeList;
  State::MAP    = IndexMAP;
  State::Weight = WeightMAP; 
  State::Wsum   = wsum;
  State::c      = kmean;
 

}
///////////////////////////////////////////////////////////////////////////////////////////////// 


/////////////////////////////////////////////////////////////////////////////////////////////////
void Get_N_E_NameList(string file, int& num, int& numedge, vector<string>& namelist)  
{
  ifstream fin(file.c_str(), ios_base::in);
  if(!fin) {cout << "Cannot open " << file  << " for read."; exit(0);}

  //int num, numedge;
  double pmin, p1min, pmax, wmin, wmax;
  string s;
  getline(fin, s, '\n');
  fin >> num >> numedge >> pmin >> p1min >> pmax >> wmin >> wmax;
  getline(fin, s, '\n'); // to finish reading this line ended by '\n'
  cout << "N      E      pmin     p1min     pmax    wmin     wmax " << endl;
  cout << num << ' ' << numedge << ' ' << pmin << ' ' << p1min << ' ' << pmax << ' ' << wmin << ' ' << wmax << endl;
 
  int index;
  string name;
  double pvalue;
  //vector<string> namelist;
  getline(fin, s, '\n');
  for(int i=0;i<num;i++) {
    fin >> index >> name >> pvalue;
    getline(fin, s, '\n'); // to finish reading this line ended by '\n'
    namelist.push_back(name);
  }
  fin.close();

}
///////////////////////////////////////////////////////////////////////////////////////////////// 
 







//////////////////////////////////////////////////////////////////////////////////////////////////////
inline hType rounding(double x)  // round the randomField to an integer (for hi_pr algorithm) 
{return  (hType)floor(RESOLUTION*x);}

void State::SetEffectiveField(double H)
{
  // If H is too positice or too negative, we may have overflow problem. (Because we have defined RESOLUTION to be a very big number. And hType = long long int, which has a clear upper limit.
  // YYL 04/17/2014
  // This problem will be serious if we just care about the ground state at a given H.
  // If we calculate the whole M-H curve, this problem typically will not show up, because we have taken care of the hmin and hmax values into account. The cross field Hx should always be well behaved.
  // unless there is something wrong with the energy calculation. See my debug note 04/17/2014
 
  /*
  if(H > -hmin) {
    //cout << "H > -hmin. The external field is too positive! We can reset it to be H= -hmin * 1.01.\n";
    H = -hmin * 1.01;
  }
  else if (H < -hmax) {
    //cout << "H < -hmax. The external field is too negative! We can reset it to be H= -hmax * 1.01.\n";
    H = -hmax * 1.01;
  }
  */

  Hext = H;
  
  // type-I
  for(int i=0; i<N; i++) 
    heff[i] = rounding(Hext + h[i]);     
}
//////////////////////////////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////////////////////////////
void State::SetEffectiveField_typeII(double lambda)
{
  Hext = 0;

  // type-II
  for(int i=0; i<N; i++) 
    heff[i] = rounding(lambda*h[i]);     
  
}
//////////////////////////////////////////////////////////////////////////////////////////////////////
