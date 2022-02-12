//////////////////////////////////////////////////////////////////////////////////////////////////////
// init Window by taking central chunk of size W, L, H   
Window::Window( istream& f)
{
    int tmpstore;
    f >> W >> L >> H;   // note that D <= 3. The Droplet code cannot be applied to higher dimension. 
    N = W * L * H;

    if(H==1)
	if(L==1) D = 1;
        else     D = 2;
    else D = 3;

    Z = 2*D;            // number of nearest neighbors

    size = new int[D];
    stride = new int[D];
    neighborLocs = new int[Z];
    neighborCoords = new int* [Z];
    pseudoneighborCoords = new int* [Z];

    switch(D)
    {
	case 1: 
	    size[0]   = W;    stride[0] = 1;  //x
	    break;

	case 2:
	    size[1]   = W;    stride[1] = 1;  //x
	    size[0]   = L;    stride[0] = W;  //y
	    break;

	case 3: 
	    size[2]   = W;    stride[2] = 1;  //x
	    size[1]   = L;    stride[1] = W;  //y
	    size[0]   = H;    stride[0] = W*L;//z
	    break;
    }
	
    for(int i=0;i<Z;i++)
    {
	neighborCoords[i] = new int[D];
	pseudoneighborCoords[i] = new int[D];
    }

    spin = new char [N];
    label = NULL;

    int mark=0;
    for (int k = 0; k < H; ++k)
	for (int j = 0; j < L; ++j)//for (int j = L-1; j >=0; j--) //
	    for (int i = 0; i < W; ++i)
	    {
		f >> tmpstore;
		spin[i + W * (j + L * k)] = tmpstore;
		if(tmpstore==1 && mark==0)
		{
		    startingSpinLoc = i + W * (j + L * k);
		    mark = 1;
		}
	    }

    //Note that in the data file, the lower line corresponds
    //to higher value of L. 
}
//////////////////////////////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////////////////////////////
// init Window by taking the valued cluster(a separated avalanche)                     //Y.L. 06/15/05
Window::Window(int dimension, int length, list<int>& spinList)
{
    D= dimension;

    switch(D)
    {
	case 1:
            W = length; L = H = 1;
	    break;
	case 2:    
	    W = L = length;  H = 1;	    
	    break;
	case 3:    
	    W = L = H = length;	    
	    break;
    }
    
    N = W * L * H;
    Z = 2*D;

    size = new int[D];
    stride = new int[D];
    neighborLocs = new int[Z];
    neighborCoords = new int* [Z];
    pseudoneighborCoords = new int* [Z];

    switch(D)
    {
	case 1: 
	    size[0]   = W;    stride[0] = 1;  //x
	    break;

	case 2:
	    size[1]   = W;    stride[1] = 1;  //x
	    size[0]   = L;    stride[0] = W;  //y
	    break;

	case 3: 
	    size[2]   = W;    stride[2] = 1;  //x
	    size[1]   = L;    stride[1] = W;  //y
	    size[0]   = H;    stride[0] = W*L;//z
	    break;
    }
	
    for(int i=0;i<Z;i++)
    {
	neighborCoords[i] = new int[D];
	pseudoneighborCoords[i] = new int[D];
    }

    spin = new char [N];
    label = NULL;

    for(int i=0; i<N; i++) 
	spin[i] = SPINDN;

    list<int>::iterator p = spinList.begin(); //the starting spin when scanning the avalanche  
    startingSpinLoc = (*p);
    
    for(; p != spinList.end(); p++)
	spin[(*p)] = SPINUP;

}
//////////////////////////////////////////////////////////////////////////////////////////////////////




//////////////////////////////////////////////////////////////////////////////////////////////////////
// init Window by taking the valued array(a separated avalanche)                      // Y.L. 06/15/05
Window::Window(int dimension, int length, char* spinCluster)
{
    D= dimension;

    switch(D)
    {
	case 1:
            W = length; L = H = 1;
	    break;
	case 2:    
	    W = L = length;  H = 1;	    
	    break;
	case 3:    
	    W = L = H = length;	    
	    break;
    }

    N = W * L * H;

    Z = 2*D;

    size = new int[D];
    stride = new int[D];
    neighborLocs = new int[Z];
    neighborCoords = new int* [Z];
    pseudoneighborCoords = new int* [Z];


    switch(D)
    {
	case 1: 
	    size[0]   = W;    stride[0] = 1;  //x
	    break;

	case 2:
	    size[1]   = W;    stride[1] = 1;  //x
	    size[0]   = L;    stride[0] = W;  //y
	    break;

	case 3: 
	    size[2]   = W;    stride[2] = 1;  //x
	    size[1]   = L;    stride[1] = W;  //y
	    size[0]   = H;    stride[0] = W*L;//z
	    break;
    }


    for(int i=0;i<Z;i++)
    {
	neighborCoords[i] = new int[D];
	pseudoneighborCoords[i] = new int[D];
    }

    spin = new char [N];
    label = NULL;

    int mark=0;
    for(int i=0; i<N; i++) 
    {
	spin[i] = spinCluster[i];
	if(spin[i]==1 && mark==0)
	{
	    startingSpinLoc = i;
	    mark = 1;
	}
    }
}
//////////////////////////////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////////////////////////////
// constructor: init window from another window 
Window::Window(const Window& x)                       
{
    
    W = x.W;
    L = x.L;
    H = x.H;
    N = x.N;
    D = x.D;
    startingSpinLoc = x.startingSpinLoc;


    Z = 2*D;

    size = new int[D];
    stride = new int[D];
    neighborLocs = new int[Z];
    neighborCoords = new int* [Z];
    pseudoneighborCoords = new int* [Z];

    switch(D)
    {
	case 1: 
	    size[0]   = W;    stride[0] = 1;  //x
	    break;

	case 2:
	    size[1]   = W;    stride[1] = 1;  //x
	    size[0]   = L;    stride[0] = W;  //y
	    break;

	case 3: 
	    size[2]   = W;    stride[2] = 1;  //x
	    size[1]   = L;    stride[1] = W;  //y
	    size[0]   = H;    stride[0] = W*L;//z
	    break;
    }

    for(int i=0;i<Z;i++)
    {
	neighborCoords[i] = new int[D];
	pseudoneighborCoords[i] = new int[D];
    }

    spin = new char[N];
    for (int i = 0; i < N; ++i)
	spin[i] = x.spin[i];
}
////////////////////////////////////////////////////////////////////////////////////////////////////// 
   



//////////////////////////////////////////////////////////////////////////////////////////////////////
// constructor: init window from a state
Window::Window(const State& x)                       
{
    int length=x.GetL();
    D = x.GetD();
    N = x.GetN();
    Z = x.GetZ();

    switch(D)
    {
	case 1:
            W = length; L = H = 1;
	    break;
	case 2:    
	    W = L = length;  H = 1;	    
	    break;
	case 3:    
	    W = L = H = length;	    
	    break;
    }


    size = new int[D];
    stride = new int[D];
    neighborLocs = new int[Z];
    neighborCoords = new int* [Z];
    pseudoneighborCoords = new int* [Z];


    switch(D)
    {
	case 1: 
	    size[0]   = W;    stride[0] = 1;  //x
	    break;

	case 2:
	    size[1]   = W;    stride[1] = 1;  //x
	    size[0]   = L;    stride[0] = W;  //y
	    break;

	case 3: 
	    size[2]   = W;    stride[2] = 1;  //x
	    size[1]   = L;    stride[1] = W;  //y
	    size[0]   = H;    stride[0] = W*L;//z
	    break;
    }


    for(int i=0;i<Z;i++)
    {
	neighborCoords[i] = new int[D];
	pseudoneighborCoords[i] = new int[D];
    }

    spin = new char [N];
    label = NULL;

    int mark=0;
    for(int i=0; i<N; i++) 
    {
	spin[i] = x.GetSpin(i);
	if(spin[i]==1 && mark==0)
	{
	    startingSpinLoc = i;
	    mark = 1;
	}
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////// 




//////////////////////////////////////////////////////////////////////////////////////////////////////
void Window::printme()
{

    cout << "\nSpin configuration:  \n";  //06/16/05 by Yang Liu

    for (int k = 0; k < H; ++k) {
	for (int j = 0; j < L; ++j) {//for (int j = L-1; j >=0; --j) {	//for (int j = 0; j < L; ++j) {
	    for (int i = 0; i < W; ++i) {
		if (spin[i + W * (j + L * k)] == SPINUP)
		    cout << '1';
		else
		    cout << '0';
	    }
	    cout << "\n";
	}
	cout << "\n";
    }
    //Note that in the data file, the lower line corresponds to higher value of L. 
    cout << "Dimension       = " << D << " and StartingSpinLoc = " << startingSpinLoc << endl;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////
