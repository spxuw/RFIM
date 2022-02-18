

Window::Window( istream& f)
{
    int tmpstore;
    f >> W >> L >> H;   
    N = W * L * H;

    if(H==1)
	if(L==1) D = 1;
        else     D = 2;
    else D = 3;

    Z = 2*D;            

    size = new int[D];
    stride = new int[D];
    neighborLocs = new int[Z];
    neighborCoords = new int* [Z];
    pseudoneighborCoords = new int* [Z];

    switch(D)
    {
	case 1: 
	    size[0]   = W;    stride[0] = 1;  
	    break;

	case 2:
	    size[1]   = W;    stride[1] = 1;  
	    size[0]   = L;    stride[0] = W;  
	    break;

	case 3: 
	    size[2]   = W;    stride[2] = 1;  
	    size[1]   = L;    stride[1] = W;  
	    size[0]   = H;    stride[0] = W*L;
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
	for (int j = 0; j < L; ++j)
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

    
    
}





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
	    size[0]   = W;    stride[0] = 1;  
	    break;

	case 2:
	    size[1]   = W;    stride[1] = 1;  
	    size[0]   = L;    stride[0] = W;  
	    break;

	case 3: 
	    size[2]   = W;    stride[2] = 1;  
	    size[1]   = L;    stride[1] = W;  
	    size[0]   = H;    stride[0] = W*L;
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

    list<int>::iterator p = spinList.begin(); 
    startingSpinLoc = (*p);
    
    for(; p != spinList.end(); p++)
	spin[(*p)] = SPINUP;

}







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
	    size[0]   = W;    stride[0] = 1;  
	    break;

	case 2:
	    size[1]   = W;    stride[1] = 1;  
	    size[0]   = L;    stride[0] = W;  
	    break;

	case 3: 
	    size[2]   = W;    stride[2] = 1;  
	    size[1]   = L;    stride[1] = W;  
	    size[0]   = H;    stride[0] = W*L;
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
	    size[0]   = W;    stride[0] = 1;  
	    break;

	case 2:
	    size[1]   = W;    stride[1] = 1;  
	    size[0]   = L;    stride[0] = W;  
	    break;

	case 3: 
	    size[2]   = W;    stride[2] = 1;  
	    size[1]   = L;    stride[1] = W;  
	    size[0]   = H;    stride[0] = W*L;
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
	    size[0]   = W;    stride[0] = 1;  
	    break;

	case 2:
	    size[1]   = W;    stride[1] = 1;  
	    size[0]   = L;    stride[0] = W;  
	    break;

	case 3: 
	    size[2]   = W;    stride[2] = 1;  
	    size[1]   = L;    stride[1] = W;  
	    size[0]   = H;    stride[0] = W*L;
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







void Window::printme()
{

    cout << "\nSpin configuration:  \n";  

    for (int k = 0; k < H; ++k) {
	for (int j = 0; j < L; ++j) {
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
    
    cout << "Dimension       = " << D << " and StartingSpinLoc = " << startingSpinLoc << endl;
}

