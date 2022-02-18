





int Window::LY_labelCC(int loc, int tolabel, char s, 
		       int* surface, int* boxvol, int* maxdim)
{
    int locn;

    int* max_dim = new int[D];     
    int** mark_del = new int* [D]; 
    for(int i=0;i<D;i++)
	mark_del[i] = new int[size[i]];

    for(int i=0;i<D;i++)
    {
	max_dim[i] = 0;
	for(int j=0;j<size[i];j++)
	    mark_del[i][j] = 0;
    }

    int headQ;                     
    queue<int> Queue;           
    Queue.push(loc);              
    
    label[loc] = tolabel;
    spin[loc] |= GRAY;	             

    int* startingcoords = new int[D];
    GetCoords(loc, startingcoords);   

    *surface = 0;
    int volume = 0;

    int* headcoords = new int[D];   
 
    while (!Queue.empty())     
    {
	headQ = Queue.front();  
	Queue.pop();         

	GetCoords(headQ, headcoords); 
	GetNeighbors(headQ, headcoords);  

	
	for(int i=D-1;i>=0;i--)
	{
	    
	    
            
	    int temp = ((headcoords[i] - startingcoords[i] + size[i] + size[i]/2) % size[i]); 
	    
	    if(mark_del[i][temp]==0) 
	    {
		max_dim[i]++; 
		mark_del[i][temp]=1;
	    }
	}

	
	for(int i=0;i<Z;i++) 
	{
	    locn = neighborLocs[i]; 
	    if ((spin[locn] & 1) != s)  
		(*surface)++;
	    else if((spin[locn] & NOTWHITE) == 0)
	    {
		spin[locn] |= GRAY; 
		label[locn] = tolabel; 
		Queue.push(locn);    
	    }
	}


	spin[headQ] -= GRAY; 
	spin[headQ] |= BLACK;
	volume++;
    }

    *maxdim = max_dim[D-1]; 
    for (int i=D-2;i>=0;i--)
	if(max_dim[i] > *maxdim) 
	    *maxdim = max_dim[i];

    *boxvol = max_dim[D-1];
    for (int i=D-2;i>=0;i--)
	*boxvol *= max_dim[i];  


    delete [] max_dim;
    for(int i=0;i<D;i++)
	delete [] mark_del[i];
    delete [] mark_del;
    delete [] startingcoords;
    delete [] headcoords;

 
    return volume;
}







int Window::LY_labelCC(int loc, int tolabel, char s, 
		       int* surface, int* boxvol, int* maxdim, 
		       bool* BoundaryTouching)
{
    int locn;

    int* max_dim = new int[D];     
    int** mark_del = new int* [D]; 
    for(int i=0;i<D;i++)
	mark_del[i] = new int[size[i]];

    for(int i=0;i<D;i++)
    {
	max_dim[i] = 0;
	for(int j=0;j<size[i];j++)
	    mark_del[i][j] = 0;
    }

    int headQ;                     
    queue<int> Queue;           
    Queue.push(loc);              
    
    label[loc] = tolabel;
    spin[loc] |= GRAY;	             

    int* startingcoords = new int[D];
    GetCoords(loc, startingcoords);   

    *surface = 0;
    int volume = 0;

    int* headcoords = new int[D];   
    bool TouchingB=false;
 
    while (!Queue.empty())     
    {
	headQ = Queue.front();  
	Queue.pop();         

	GetCoords(headQ, headcoords); 
	GetNeighbors(headQ, headcoords);  

	
	for(int i=D-1;i>=0;i--)
	{
	    
	    
            
	    int temp = ((headcoords[i] - startingcoords[i] + size[i] + size[i]/2) % size[i]); 
	    
	    if(mark_del[i][temp]==0) 
	    {
		max_dim[i]++; 
		mark_del[i][temp]=1;
	    }
	}

	
	for(int i=0;i<Z;i++) 
	{
	    locn = neighborLocs[i]; 
	    if ((spin[locn] & 1) != s)  
		(*surface)++;
	    else if((spin[locn] & NOTWHITE) == 0)
	    {
		spin[locn] |= GRAY; 
		label[locn] = tolabel; 
		Queue.push(locn);    
	    }
	}


	spin[headQ] -= GRAY; 
	spin[headQ] |= BLACK;
	volume++;

	
        
	
	if(!TouchingB)
	    for(int i=D-1;i>=0;i--)
		if( (headcoords[i]==0)||(headcoords[i]==size[i]-1) )
		    TouchingB = true;


    }

    *maxdim = max_dim[D-1]; 
    for (int i=D-2;i>=0;i--)
	if(max_dim[i] > *maxdim) 
	    *maxdim = max_dim[i];

    *boxvol = max_dim[D-1];
    for (int i=D-2;i>=0;i--)
	*boxvol *= max_dim[i];  

    *BoundaryTouching = TouchingB;


    delete [] max_dim;
    for(int i=0;i<D;i++)
	delete [] mark_del[i];
    delete [] mark_del;
    delete [] startingcoords;
    delete [] headcoords;

 
    return volume;
}









int Window::LY_labelCC_getQ(int loc, int tolabel, char s, 
			    int* surface, int* boxvol, int* maxdim, 
			    Vec_DP* EVinfo)
{
    int locn;	  

    int* max_dim = new int[D];     
    int** mark_del = new int* [D]; 
    for(int i=0;i<D;i++)
	mark_del[i] = new int[size[i]];

    for(int i=0;i<D;i++)
    {
	max_dim[i] = 0;
	for(int j=0;j<size[i];j++)
	    mark_del[i][j] = 0;
    }
  

    int headQ;
    queue<int> locQueue;
    locQueue.push(loc);              

    label[loc] = tolabel;
    spin[loc] |= GRAY;		

    int* startingcoords = new int[D];
    GetCoords(loc, startingcoords);   

    *surface = 0;
    int volume = 0;

    list<Point>  pointList;   
    queue<Point> pointQueue;  
    int* coords = new int[D];
    GetCoords(loc, coords);             
    Point point(D, coords);                
    pointList.push_back(point);    
    pointQueue.push(point);        
    
    int* headcoords = new int[D];   

    while (!locQueue.empty())     
    {
	headQ = locQueue.front();              
	locQueue.pop();                     
	GetCoords(headQ, headcoords); 
	GetNeighbors(headQ, headcoords);  

	Point point = pointQueue.front();  
 	pointQueue.pop();                  
	GetPseudoNeighbors(point);         


	
	for(int i=D-1;i>=0;i--)
	{
	    
	    
            
	    int temp = ((headcoords[i] - startingcoords[i] + size[i] + size[i]/2) % size[i]); 
	    
	    if(mark_del[i][temp]==0) 
	    {
		max_dim[i]++; 
		mark_del[i][temp]=1;
	    }
	}

	
	
	
	for(int i=0;i<Z;i++) 
	{
	    
	    locn = neighborLocs[i]; 

	    if ((spin[locn] & 1) != s)  
		(*surface)++;
	    else if((spin[locn] & NOTWHITE) == 0)
	    {
		spin[locn] |= GRAY; 
		label[locn] = tolabel; 
		locQueue.push(locn);    

		Point point(D, pseudoneighborCoords[i]);           
		pointList.push_back(point);                        
		pointQueue.push(point);
	    }
	}

	spin[headQ] -= GRAY; 
	spin[headQ] |= BLACK;
	volume++;

    }

    *maxdim = max_dim[D-1]; 
    for (int i=D-2;i>=0;i--)
	if(max_dim[i] > *maxdim) 
	    *maxdim = max_dim[i];

    *boxvol = max_dim[D-1];
    for (int i=D-2;i>=0;i--)
	*boxvol *= max_dim[i];  

    
    *EVinfo = get_eigenvalue(pointList);

    delete [] max_dim;
    for(int i=0;i<D;i++)
	delete [] mark_del[i];
    delete [] mark_del;
    delete [] startingcoords;
    delete [] headcoords;
    delete [] coords;
    return volume;
}













int Window::LY_labelCC_getQ(int loc, int tolabel, char s, 
			    int* surface, int* boxvol, int* maxdim, 
			    Vec_DP* EVinfo, 
			    bool* BoundaryTouching)
{
    int locn;	  

    int* max_dim = new int[D];     
    int** mark_del = new int* [D]; 
    for(int i=0;i<D;i++)
	mark_del[i] = new int[size[i]];

    for(int i=0;i<D;i++)
    {
	max_dim[i] = 0;
	for(int j=0;j<size[i];j++)
	    mark_del[i][j] = 0;
    }
  

    int headQ;
    queue<int> locQueue;
    locQueue.push(loc);              

    label[loc] = tolabel;
    spin[loc] |= GRAY;		

    int* startingcoords = new int[D];
    GetCoords(loc, startingcoords);   

    *surface = 0;
    int volume = 0;

    list<Point>  pointList;   
    queue<Point> pointQueue;  
    int* coords = new int[D];
    GetCoords(loc, coords);             
    Point point(D, coords);                
    pointList.push_back(point);    
    pointQueue.push(point);        
    
    int* headcoords = new int[D];   

    bool TouchingB=false;

    while (!locQueue.empty())     
    {
	headQ = locQueue.front();              
	locQueue.pop();                     
	GetCoords(headQ, headcoords); 
	GetNeighbors(headQ, headcoords);  

	Point point = pointQueue.front();  
 	pointQueue.pop();                  
	GetPseudoNeighbors(point);         


	
	for(int i=D-1;i>=0;i--)
	{
	    
	    
            
	    int temp = ((headcoords[i] - startingcoords[i] + size[i] + size[i]/2) % size[i]); 
	    
	    if(mark_del[i][temp]==0) 
	    {
		max_dim[i]++; 
		mark_del[i][temp]=1;
	    }
	}

	
	
	
	for(int i=0;i<Z;i++) 
	{
	    
	    locn = neighborLocs[i]; 

	    if ((spin[locn] & 1) != s)  
		(*surface)++;
	    else if((spin[locn] & NOTWHITE) == 0)
	    {
		spin[locn] |= GRAY; 
		label[locn] = tolabel; 
		locQueue.push(locn);    

		Point point(D, pseudoneighborCoords[i]);           
		pointList.push_back(point);                        
		pointQueue.push(point);
	    }
	}

	spin[headQ] -= GRAY; 
	spin[headQ] |= BLACK;
	volume++;

	
        
	
	if(!TouchingB)
	    for(int i=D-1;i>=0;i--)
		if( (headcoords[i]==0)||(headcoords[i]==size[i]-1) )
		    TouchingB = true;

    }

    *maxdim = max_dim[D-1]; 
    for (int i=D-2;i>=0;i--)
	if(max_dim[i] > *maxdim) 
	    *maxdim = max_dim[i];

    *boxvol = max_dim[D-1];
    for (int i=D-2;i>=0;i--)
	*boxvol *= max_dim[i];  

    
    
    
    *EVinfo = get_eigenvalue(pointList,volume); 

    *BoundaryTouching = TouchingB;

    delete [] max_dim;
    for(int i=0;i<D;i++)
	delete [] mark_del[i];
    delete [] mark_del;
    delete [] startingcoords;
    delete [] headcoords;
    delete [] coords;
    return volume;
}











int Window::LY_labelCC(int loc, int tolabel, char s, bool* BoundaryTouching)
{
    int locn;	  

    int headQ;
    queue<int> locQueue;
    locQueue.push(loc);              

    label[loc] = tolabel;
    spin[loc] |= GRAY;		

    int volume = 0;
    int* headcoords = new int[D];   
    bool TouchingB=false;

    while (!locQueue.empty())     
    {
	headQ = locQueue.front();              
	locQueue.pop();                     
	GetCoords(headQ, headcoords); 
	GetNeighbors(headQ);  

	
	
	
	for(int i=0;i<Z;i++) 
	{
	    
	    locn = neighborLocs[i]; 

	    if ((spin[locn] & 1) != s)  
		continue; 
	    else if((spin[locn] & NOTWHITE) == 0)
	    {
		spin[locn] |= GRAY; 
		label[locn] = tolabel; 
		locQueue.push(locn);    
	    }
	}

	spin[headQ] -= GRAY; 
	spin[headQ] |= BLACK;
	volume++;

	
        
	
	if(!TouchingB)
	    for(int i=D-1;i>=0;i--)
		if( (headcoords[i]==0)||(headcoords[i]==size[i]-1) )
		    TouchingB = true;

    }

    *BoundaryTouching = TouchingB;

    delete [] headcoords;
    return volume;
}









int Window::LY_labelCC_getQ(int loc, int tolabel, char s,
			    int* surface, int* boxvol, int* maxdim, 
			    Vec_DP* EVinfo, bool* BoundaryTouching, 
			    list<Point>* newList, vector<int>* newSize)
{
    int locn;	  

    vector<int> min_coord(D);     
    vector<int> max_dim(D);       

    int** mark_del = new int* [D]; 
    for(int i=0;i<D;i++)
	mark_del[i] = new int[size[i]];

    for(int i=0;i<D;i++)
    {
	max_dim[i] = 0;
	for(int j=0;j<size[i];j++)
	    mark_del[i][j] = 0;
    }
  

    int headQ;
    queue<int> locQueue;
    locQueue.push(loc);              

    label[loc] = tolabel;
    spin[loc] |= GRAY;		

    int* startingcoords = new int[D];
    GetCoords(loc, startingcoords);   

    *surface = 0;
    int volume = 0;

    list<Point>  pointList;   
    queue<Point> pointQueue;  

    int* coords = new int[D];
    GetCoords(loc, coords);             

    for(int i=0;i<D;i++)                
	min_coord[i] = coords[i];

    Point point(D, coords);                
    pointList.push_back(point);    
    pointQueue.push(point);        

    int* headcoords = new int[D];   

    bool TouchingB=false;

    while (!locQueue.empty())     
    {
	headQ = locQueue.front();              
	locQueue.pop();                     
	GetCoords(headQ, headcoords); 
	GetNeighbors(headQ, headcoords);  

	Point point = pointQueue.front();  
 	pointQueue.pop();                  
	GetPseudoNeighbors(point);         


	
	for(int i=D-1;i>=0;i--)
	{
	    
	    
            
	    int u = ((headcoords[i] - startingcoords[i] + size[i] + size[i]/2) % size[i]); 
	    
	    if(mark_del[i][u]==0) 
	    {
		max_dim[i]++; 
		mark_del[i][u]=1;
	    }

	    int v = point.coord(D,i);
	    if(min_coord[i]>v)
		min_coord[i]=v;

	}

	
	
	
	for(int i=0;i<Z;i++) 
	{
	    
	    locn = neighborLocs[i]; 

	    if ((spin[locn] & 1) != s)  
		(*surface)++;
	    else if((spin[locn] & NOTWHITE) == 0)
	    {
		spin[locn] |= GRAY; 
		label[locn] = tolabel; 
		locQueue.push(locn);    

		Point point(D, pseudoneighborCoords[i]);           
		pointList.push_back(point);                        
		pointQueue.push(point);


	    }
	}

	spin[headQ] -= GRAY; 
	spin[headQ] |= BLACK;
	volume++;

	
        
	
	if(!TouchingB)
	    for(int i=D-1;i>=0;i--)
		if( (headcoords[i]==0)||(headcoords[i]==size[i]-1) )
		    TouchingB = true;

    }

    *maxdim = max_dim[D-1]; 
    for (int i=D-2;i>=0;i--)
	if(max_dim[i] > *maxdim) 
	    *maxdim = max_dim[i];

    *boxvol = max_dim[D-1];
    for (int i=D-2;i>=0;i--)
	*boxvol *= max_dim[i];  

    
    *EVinfo = get_eigenvalue(pointList);
    *BoundaryTouching = TouchingB;


    
    
    
    
    
    vector<int> deltCoord(D);
    vector<int> sideSize(D);
    for(int i=0;i<D;i++){
	deltCoord[i] = min_coord[i] - 1;
	sideSize[i] = max_dim[i] + 2; 
    }

    for(list<Point>::iterator p = pointList.begin(); p != pointList.end(); p++)
	(*p).Shift(D, deltCoord);

    *newList= pointList; 
    *newSize= sideSize;


    for(int i=0;i<D;i++)
	delete [] mark_del[i];
    delete [] mark_del;
    delete [] startingcoords;
    delete [] headcoords;
    delete [] coords;
    return volume;
}

