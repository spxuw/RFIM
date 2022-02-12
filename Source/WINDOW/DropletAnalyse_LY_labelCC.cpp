
//////////////////////////////////////////////////////////////////////////////////////////////////////
//This function will return volume of the connected cluster. the starting point with location loc, has spin s, 
//index idx, and finally all spins in this connected cluster will be labelled with tolabel! 
//The surfacearea, boxvolume and maxdim will also be passed back via pointers. Y.L.02/27/06

int Window::LY_labelCC(int loc, int tolabel, char s, 
		       int* surface, int* boxvol, int* maxdim)
{
    int locn;

    int* max_dim = new int[D];     // max dim of in x,y,z... directions (D dimensional)
    int** mark_del = new int* [D]; // check whether in a certain direction, a certain point has been counted
    for(int i=0;i<D;i++)
	mark_del[i] = new int[size[i]];

    for(int i=0;i<D;i++)
    {
	max_dim[i] = 0;
	for(int j=0;j<size[i];j++)
	    mark_del[i][j] = 0;
    }

    int headQ;                     // head (location) of the queue
    queue<int> Queue;           
    Queue.push(loc);              // push it into the real loc list
    
    label[loc] = tolabel;
    spin[loc] |= GRAY;	             // "gray"

    int* startingcoords = new int[D];
    GetCoords(loc, startingcoords);   

    *surface = 0;
    int volume = 0;

    int* headcoords = new int[D];   
 
    while (!Queue.empty())     
    {
	headQ = Queue.front();  // get the real loc from the list
	Queue.pop();         // pop it from the list

	GetCoords(headQ, headcoords); // get its coordinates stored in coords[i] 
	GetNeighbors(headQ, headcoords);  // get its real neighbors' locs stored in neighborLocs[i] 

	//   calculate the distance between coords and startingcoords
	for(int i=D-1;i>=0;i--)
	{
	    // del[i]= ((headcoords[i] - startingcoords[i] + size[i] + size[i]/2) % size[i]) - size[i]/2
	    // since del[i] (the relative_coords) is in the range [-size[i]/2, size[i]/2)
            // so temp= (del[i] + size[i]/2) is in the range [0,size[i]) 
	    int temp = ((headcoords[i] - startingcoords[i] + size[i] + size[i]/2) % size[i]); //relative_coords
	    
	    if(mark_del[i][temp]==0) 
	    {
		max_dim[i]++; 
		mark_del[i][temp]=1;
	    }
	}

	// now look at nearest neighbor: if white, same spin, then color gray and queue
	for(int i=0;i<Z;i++) 
	{
	    locn = neighborLocs[i]; 
	    if ((spin[locn] & 1) != s)  //don't include!
		(*surface)++;//but increment surface
	    else if((spin[locn] & NOTWHITE) == 0)
	    {
		spin[locn] |= GRAY; 
		label[locn] = tolabel; 
		Queue.push(locn);    
	    }
	}// end of for each neighbor loop


	spin[headQ] -= GRAY; // better have this bit set! (should by algorithm)
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
//////////////////////////////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////////////////////////////
// The following function is overloaded from the 
// LY_labelCC(int loc, int tolabel, char s, int* surface, int* boxvol, int* maxdim)
// The only difference is that it can tell whether the connected cluster is boundary touching. 03/14/06
int Window::LY_labelCC(int loc, int tolabel, char s, 
		       int* surface, int* boxvol, int* maxdim, 
		       bool* BoundaryTouching)
{
    int locn;

    int* max_dim = new int[D];     // max dim of in x,y,z... directions (D dimensional)
    int** mark_del = new int* [D]; // check whether in a certain direction, a certain point has been counted
    for(int i=0;i<D;i++)
	mark_del[i] = new int[size[i]];

    for(int i=0;i<D;i++)
    {
	max_dim[i] = 0;
	for(int j=0;j<size[i];j++)
	    mark_del[i][j] = 0;
    }

    int headQ;                     // head (location) of the queue
    queue<int> Queue;           
    Queue.push(loc);              // push it into the real loc list
    
    label[loc] = tolabel;
    spin[loc] |= GRAY;	             // "gray"

    int* startingcoords = new int[D];
    GetCoords(loc, startingcoords);   

    *surface = 0;
    int volume = 0;

    int* headcoords = new int[D];   
    bool TouchingB=false;
 
    while (!Queue.empty())     
    {
	headQ = Queue.front();  // get the real loc from the list
	Queue.pop();         // pop it from the list

	GetCoords(headQ, headcoords); // get its coordinates stored in coords[i] 
	GetNeighbors(headQ, headcoords);  // get its real neighbors' locs stored in neighborLocs[i] 

	//   calculate the distance between coords and startingcoords
	for(int i=D-1;i>=0;i--)
	{
	    // del[i]= ((headcoords[i] - startingcoords[i] + size[i] + size[i]/2) % size[i]) - size[i]/2
	    // since del[i] (the relative_coords) is in the range [-size[i]/2, size[i]/2)
            // so temp= (del[i] + size[i]/2) is in the range [0,size[i]) 
	    int temp = ((headcoords[i] - startingcoords[i] + size[i] + size[i]/2) % size[i]); //relative_coords
	    
	    if(mark_del[i][temp]==0) 
	    {
		max_dim[i]++; 
		mark_del[i][temp]=1;
	    }
	}

	// now look at nearest neighbor: if white, same spin, then color gray and queue
	for(int i=0;i<Z;i++) 
	{
	    locn = neighborLocs[i]; 
	    if ((spin[locn] & 1) != s)  //don't include!
		(*surface)++;//but increment surface
	    else if((spin[locn] & NOTWHITE) == 0)
	    {
		spin[locn] |= GRAY; 
		label[locn] = tolabel; 
		Queue.push(locn);    
	    }
	}// end of for each neighbor loop


	spin[headQ] -= GRAY; // better have this bit set! (should by algorithm)
	spin[headQ] |= BLACK;
	volume++;

	// check whether this point is a boundary point
        // if it was, then this connected cluster is defined as 
	// a boundary touching cluster!
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
//////////////////////////////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////////////////////////////
// The following function is almost the same as the function 
//        int Window::LY_labelCC(int idx, int tolabel, char s, 
//                               int* surface, int* boxvol, int* maxdim)
// The only difference is that it can pass back information about the shape tensor. Y.L 02/25/06
//////////////////////////////////////////////////////////////////////////////////////////////////////
int Window::LY_labelCC_getQ(int loc, int tolabel, char s, 
			    int* surface, int* boxvol, int* maxdim, 
			    Vec_DP* EVinfo)
{
    int locn;	  // location(index) of nearest neighbor

    int* max_dim = new int[D];     // max dim of in x,y,z... directions (D dimensional)
    int** mark_del = new int* [D]; // check whether in a certain direction, a certain point has been counted
    for(int i=0;i<D;i++)
	mark_del[i] = new int[size[i]];

    for(int i=0;i<D;i++)
    {
	max_dim[i] = 0;
	for(int j=0;j<size[i];j++)
	    mark_del[i][j] = 0;
    }
  

    int headQ;// head (location) of the queue
    queue<int> locQueue;
    locQueue.push(loc);              // push it into the real loc list

    label[loc] = tolabel;
    spin[loc] |= GRAY;		// "gray"

    int* startingcoords = new int[D];
    GetCoords(loc, startingcoords);   

    *surface = 0;
    int volume = 0;

    list<Point>  pointList;   // pseudo cluster's coordinates
    queue<Point> pointQueue;  // pseudo cluster's coordinates
    int* coords = new int[D];
    GetCoords(loc, coords);             // get the coordinates of the starting spin with the function 
    Point point(D, coords);                
    pointList.push_back(point);    // push it into the virtual point queue 
    pointQueue.push(point);        // and the virtual point list   
    
    int* headcoords = new int[D];   

    while (!locQueue.empty())     //while(!pointQueue.empty())         
    {
	headQ = locQueue.front();              // get the real loc from the list
	locQueue.pop();                     // pop it from the list
	GetCoords(headQ, headcoords); // get its coordinates stored in coords[i] 
	GetNeighbors(headQ, headcoords);  // get its real neighbors' locs stored in neighborLocs[i] 

	Point point = pointQueue.front();  // get the virtual point from the queue 
 	pointQueue.pop();                  // pop it from the queue
	GetPseudoNeighbors(point);         // get its pseudo neighbors  


	//   calculate the distance between coords and startingcoords
	for(int i=D-1;i>=0;i--)
	{
	    // del[i]= ((headcoords[i] - startingcoords[i] + size[i] + size[i]/2) % size[i]) - size[i]/2
	    // since del[i] (the relative_coords) is in the range [-size[i]/2, size[i]/2)
            // so temp= (del[i] + size[i]/2) is in the range [0,size[i]) 
	    int temp = ((headcoords[i] - startingcoords[i] + size[i] + size[i]/2) % size[i]); //relative_coords
	    
	    if(mark_del[i][temp]==0) 
	    {
		max_dim[i]++; 
		mark_del[i][temp]=1;
	    }
	}

	// loop over this spin's pseudo neighbors and record them 
	// if the neighbor is also spin-up and has not been pushed into the queue 		
	// also check the nearest neighbor: if white, same spin, then color gray and queue
	for(int i=0;i<Z;i++) 
	{
	    //locn = GetLoc_PBC(pseudoneighborCoords[i]); // get pseudo neighbors' real locations. 
	    locn = neighborLocs[i]; // equivalent to the result GetLoc_PBC(pseudoneighborCoords[i])

	    if ((spin[locn] & 1) != s)  //don't include!
		(*surface)++;//but increment surface
	    else if((spin[locn] & NOTWHITE) == 0)
	    {
		spin[locn] |= GRAY; 
		label[locn] = tolabel; 
		locQueue.push(locn);    

		Point point(D, pseudoneighborCoords[i]);           
		pointList.push_back(point);                        // and push it into the queue and the list. 
		pointQueue.push(point);
	    }
	}// end of for each neighbor loop

	spin[headQ] -= GRAY; // better have this bit set! (should by algorithm)
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

    // Once we get the point list of the cluster/avalanche, the following work is trivial. !!!
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
//////////////////////////////////////////////////////////////////////////////////////////////////////






//////////////////////////////////////////////////////////////////////////////////////////////////////
// The following function is overloaded from 
// int Window::LY_labelCC_getQ(int loc, int tolabel, char s, 
//                             int* surface, int* boxvol, int* maxdim, Vec_DP* EVinfo)
// It can tell whether the connected cluster is boundary touching.
//////////////////////////////////////////////////////////////////////////////////////////////////////
int Window::LY_labelCC_getQ(int loc, int tolabel, char s, 
			    int* surface, int* boxvol, int* maxdim, 
			    Vec_DP* EVinfo, 
			    bool* BoundaryTouching)
{
    int locn;	  // location(index) of nearest neighbor

    int* max_dim = new int[D];     // max dim of in x,y,z... directions (D dimensional)
    int** mark_del = new int* [D]; // check whether in a certain direction, a certain point has been counted
    for(int i=0;i<D;i++)
	mark_del[i] = new int[size[i]];

    for(int i=0;i<D;i++)
    {
	max_dim[i] = 0;
	for(int j=0;j<size[i];j++)
	    mark_del[i][j] = 0;
    }
  

    int headQ;// head (location) of the queue
    queue<int> locQueue;
    locQueue.push(loc);              // push it into the real loc list

    label[loc] = tolabel;
    spin[loc] |= GRAY;		// "gray"

    int* startingcoords = new int[D];
    GetCoords(loc, startingcoords);   

    *surface = 0;
    int volume = 0;

    list<Point>  pointList;   // pseudo cluster's coordinates
    queue<Point> pointQueue;  // pseudo cluster's coordinates
    int* coords = new int[D];
    GetCoords(loc, coords);             // get the coordinates of the starting spin with the function 
    Point point(D, coords);                
    pointList.push_back(point);    // push it into the virtual point queue 
    pointQueue.push(point);        // and the virtual point list   
    
    int* headcoords = new int[D];   

    bool TouchingB=false;

    while (!locQueue.empty())     //while(!pointQueue.empty())         
    {
	headQ = locQueue.front();              // get the real loc from the list
	locQueue.pop();                     // pop it from the list
	GetCoords(headQ, headcoords); // get its coordinates stored in coords[i] 
	GetNeighbors(headQ, headcoords);  // get its real neighbors' locs stored in neighborLocs[i] 

	Point point = pointQueue.front();  // get the virtual point from the queue 
 	pointQueue.pop();                  // pop it from the queue
	GetPseudoNeighbors(point);         // get its pseudo neighbors  


	//   calculate the distance between coords and startingcoords
	for(int i=D-1;i>=0;i--)
	{
	    // del[i]= ((headcoords[i] - startingcoords[i] + size[i] + size[i]/2) % size[i]) - size[i]/2
	    // since del[i] (the relative_coords) is in the range [-size[i]/2, size[i]/2)
            // so temp= (del[i] + size[i]/2) is in the range [0,size[i]) 
	    int temp = ((headcoords[i] - startingcoords[i] + size[i] + size[i]/2) % size[i]); //relative_coords
	    
	    if(mark_del[i][temp]==0) 
	    {
		max_dim[i]++; 
		mark_del[i][temp]=1;
	    }
	}

	// loop over this spin's pseudo neighbors and record them 
	// if the neighbor is also spin-up and has not been pushed into the queue 		
	// also check the nearest neighbor: if white, same spin, then color gray and queue
	for(int i=0;i<Z;i++) 
	{
	    //locn = GetLoc_PBC(pseudoneighborCoords[i]); // get pseudo neighbors' real locations. 
	    locn = neighborLocs[i]; // equivalent to the result GetLoc_PBC(pseudoneighborCoords[i])

	    if ((spin[locn] & 1) != s)  //don't include!
		(*surface)++;//but increment surface
	    else if((spin[locn] & NOTWHITE) == 0)
	    {
		spin[locn] |= GRAY; 
		label[locn] = tolabel; 
		locQueue.push(locn);    

		Point point(D, pseudoneighborCoords[i]);           
		pointList.push_back(point);                        // and push it into the queue and the list. 
		pointQueue.push(point);
	    }
	}// end of for each neighbor loop

	spin[headQ] -= GRAY; // better have this bit set! (should by algorithm)
	spin[headQ] |= BLACK;
	volume++;

	// check whether this point is a boundary point
        // if it was, then this connected cluster is defined as 
	// a boundary touching cluster!
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

    // Once we get the point list of the cluster/avalanche, the following work is trivial. !!!
    // Here if volume <3, we can use pre-determined eigenvalues without doing any calculation to save time
    // Yang Liu 06/25/06
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
//////////////////////////////////////////////////////////////////////////////////////////////////////




//////////////////////////////////////////////////////////////////////////////////////////////////////
// The following function is modified from 
// int Window::LY_labelCC_getQ(int loc, int tolabel, char s, 
//                             int* surface, int* boxvol, int* maxdim, Vec_DP* EVinfo, bool* BoundaryTouching)
// It is used for those spin_DOWN connected clusters.
//////////////////////////////////////////////////////////////////////////////////////////////////////
int Window::LY_labelCC(int loc, int tolabel, char s, bool* BoundaryTouching)
{
    int locn;	  // location(index) of nearest neighbor

    int headQ;// head (location) of the queue
    queue<int> locQueue;
    locQueue.push(loc);              // push it into the real loc list

    label[loc] = tolabel;
    spin[loc] |= GRAY;		// "gray"

    int volume = 0;
    int* headcoords = new int[D];   
    bool TouchingB=false;

    while (!locQueue.empty())     //while(!pointQueue.empty())         
    {
	headQ = locQueue.front();              // get the real loc from the list
	locQueue.pop();                     // pop it from the list
	GetCoords(headQ, headcoords); // get its coordinates stored in coords[i] 
	GetNeighbors(headQ);  // get its real neighbors' locs stored in neighborLocs[i] 

	// loop over this spin's neighbors and record them 
	// if the neighbor is also spin-up and has not been pushed into the queue 		
	// also check the nearest neighbor: if white, same spin, then color gray and queue
	for(int i=0;i<Z;i++) 
	{
	    //locn = GetLoc_PBC(pseudoneighborCoords[i]); // get pseudo neighbors' real locations. 
	    locn = neighborLocs[i]; // equivalent to the result GetLoc_PBC(pseudoneighborCoords[i])

	    if ((spin[locn] & 1) != s)  //don't include!
		continue; //(*surface)++;//but increment surface
	    else if((spin[locn] & NOTWHITE) == 0)
	    {
		spin[locn] |= GRAY; 
		label[locn] = tolabel; 
		locQueue.push(locn);    
	    }
	}// end of for each neighbor loop

	spin[headQ] -= GRAY; // better have this bit set! (should by algorithm)
	spin[headQ] |= BLACK;
	volume++;

	// check whether this point is a boundary point
        // if it was, then this connected cluster is defined as 
	// a boundary touching cluster!
	if(!TouchingB)
	    for(int i=D-1;i>=0;i--)
		if( (headcoords[i]==0)||(headcoords[i]==size[i]-1) )
		    TouchingB = true;

    }

    *BoundaryTouching = TouchingB;

    delete [] headcoords;
    return volume;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////////////////////////////
// The following function is overloaded from 
// int Window::LY_labelCC_getQ(int loc, int tolabel, char s, 
//           		       int* surface, int* boxvol, int* maxdim, Vec_DP* EVinfo, bool* BoundaryTouching)
// It can pass back the minmum virtual coordinates in D directions. 03/13/06
//////////////////////////////////////////////////////////////////////////////////////////////////////
int Window::LY_labelCC_getQ(int loc, int tolabel, char s,
			    int* surface, int* boxvol, int* maxdim, 
			    Vec_DP* EVinfo, bool* BoundaryTouching, 
			    list<Point>* newList, vector<int>* newSize)
{
    int locn;	  // location(index) of nearest neighbor

    vector<int> min_coord(D);     // the minmum (virtual) coordinates in D directions, could be negative Y.L. 03/13/06
    vector<int> max_dim(D);       // max dim of in x,y,z... directions (D dimensional)

    int** mark_del = new int* [D]; // check whether in a certain direction, a certain point has been counted
    for(int i=0;i<D;i++)
	mark_del[i] = new int[size[i]];

    for(int i=0;i<D;i++)
    {
	max_dim[i] = 0;
	for(int j=0;j<size[i];j++)
	    mark_del[i][j] = 0;
    }
  

    int headQ;// head (location) of the queue
    queue<int> locQueue;
    locQueue.push(loc);              // push it into the real loc list

    label[loc] = tolabel;
    spin[loc] |= GRAY;		// "gray"

    int* startingcoords = new int[D];
    GetCoords(loc, startingcoords);   

    *surface = 0;
    int volume = 0;

    list<Point>  pointList;   // pseudo cluster's coordinates
    queue<Point> pointQueue;  // pseudo cluster's coordinates

    int* coords = new int[D];
    GetCoords(loc, coords);             // get the coordinates of the starting spin 

    for(int i=0;i<D;i++)                // initialize the min_cord[i]
	min_coord[i] = coords[i];

    Point point(D, coords);                
    pointList.push_back(point);    // push it into the virtual point queue 
    pointQueue.push(point);        // and the virtual point list   

    int* headcoords = new int[D];   

    bool TouchingB=false;

    while (!locQueue.empty())     //while(!pointQueue.empty())         
    {
	headQ = locQueue.front();              // get the real loc from the list
	locQueue.pop();                     // pop it from the list
	GetCoords(headQ, headcoords); // get its coordinates stored in coords[i] 
	GetNeighbors(headQ, headcoords);  // get its real neighbors' locs stored in neighborLocs[i] 

	Point point = pointQueue.front();  // get the virtual point from the queue 
 	pointQueue.pop();                  // pop it from the queue
	GetPseudoNeighbors(point);         // get its pseudo neighbors  


	//   calculate the distance between coords and startingcoords
	for(int i=D-1;i>=0;i--)
	{
	    // del[i]= ((headcoords[i] - startingcoords[i] + size[i] + size[i]/2) % size[i]) - size[i]/2
	    // since del[i] (the relative_coords) is in the range [-size[i]/2, size[i]/2)
            // so u= (del[i] + size[i]/2) is in the range [0,size[i]) 
	    int u = ((headcoords[i] - startingcoords[i] + size[i] + size[i]/2) % size[i]); //relative_coords
	    
	    if(mark_del[i][u]==0) 
	    {
		max_dim[i]++; 
		mark_del[i][u]=1;
	    }

	    int v = point.coord(D,i);
	    if(min_coord[i]>v)
		min_coord[i]=v;

	}

	// loop over this spin's pseudo neighbors and record them 
	// if the neighbor is also spin-up and has not been pushed into the queue 		
	// also check the nearest neighbor: if white, same spin, then color gray and queue
	for(int i=0;i<Z;i++) 
	{
	    //locn = GetLoc_PBC(pseudoneighborCoords[i]); // get pseudo neighbors' real locations. 
	    locn = neighborLocs[i]; // equivalent to the result GetLoc_PBC(pseudoneighborCoords[i])

	    if ((spin[locn] & 1) != s)  //don't include!
		(*surface)++;//but increment surface
	    else if((spin[locn] & NOTWHITE) == 0)
	    {
		spin[locn] |= GRAY; 
		label[locn] = tolabel; 
		locQueue.push(locn);    

		Point point(D, pseudoneighborCoords[i]);           
		pointList.push_back(point);                        // and push it into the queue and the list. 
		pointQueue.push(point);


	    }
	}// end of for each neighbor loop

	spin[headQ] -= GRAY; // better have this bit set! (should by algorithm)
	spin[headQ] |= BLACK;
	volume++;

	// check whether this point is a boundary point
        // if it was, then this connected cluster is defined as 
	// a boundary touching cluster!
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

    // Once we get the point list of the cluster/avalanche, the calculation of shape tensor is trivial. !!!
    *EVinfo = get_eigenvalue(pointList);
    *BoundaryTouching = TouchingB;


    // We get the point list of the cluster/avalanche, we also know the minimum (virtual) coordinates, we can 
    // then shift the cluster/avalanche in this way: (x_min,y_min,z_min) ---> (1,1,1)  
    // And the sizes of smallest fitting Window than can contain the cluster/avalanche are
    // (max_dim_x + 2, max_dim_y + 2, max_dim_z + 2)
    // Y.L.03/13/06
    vector<int> deltCoord(D);
    vector<int> sideSize(D);
    for(int i=0;i<D;i++){
	deltCoord[i] = min_coord[i] - 1;
	sideSize[i] = max_dim[i] + 2; 
    }

    for(list<Point>::iterator p = pointList.begin(); p != pointList.end(); p++)
	(*p).Shift(D, deltCoord);

    *newList= pointList; // The shiffted pointList is pass back to make a new Window.
    *newSize= sideSize;


    for(int i=0;i<D;i++)
	delete [] mark_del[i];
    delete [] mark_del;
    delete [] startingcoords;
    delete [] headcoords;
    delete [] coords;
    return volume;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////
