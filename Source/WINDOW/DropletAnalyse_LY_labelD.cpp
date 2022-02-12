//////////////////////////////////////////////////////////////////////////////////////////////////////
// The following function will return volume of the domain. It is much similar with the function: 
// int Window::LY_label_connected_cluster(int idx, int tolabel, char s, 
//                                        int* surface, int* boxvol, int* maxdim)
// But here, we just concern those spins with index label are "notlabel", we don't care whether it is 
// spin UP or DOWN. So in this way we get the domain information. The surfacearea, boxvolume and maxdim 
// will be passed back via pointers. Y.L.02/27/06
//////////////////////////////////////////////////////////////////////////////////////////////////////
int Window::LY_domain_with_notlabel(int loc, int notlabel, int *surface, int *boxvol, int *maxdim)
{
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

    *surface = 0;
    int volume = 0;
   
    int* startingcoords = new int[D];
    GetCoords(loc, startingcoords);   

    spin[loc] |= GRAY;	             // "gray"

    Queue.push(loc);              // push it into the real loc list   

    int* coords = new int[D];   
 
    while (!Queue.empty())     
    {
	headQ = Queue.front();  // get the real loc from the list
	Queue.pop();         // pop it from the list

	GetCoords(headQ, coords); 	// get its coordinates stored in coords[i] 
	GetNeighbors(headQ, coords);  	//  get its real neighbors' locs stored in neighborLocs[i] 

	//   calculate the distance between coords and startingcoords
	for(int i=D-1;i>=0;i--)
	{
	    int temp = ((coords[i] - startingcoords[i] + size[i] + size[i]/2) % size[i]); //relative_coords
	    
	    if(mark_del[i][temp]==0) 
	    {
		max_dim[i]++; 
		mark_del[i][temp]=1;
	    }
	}

	// now look at nn: if white, not labeled notlabel, then color gray and queue
	for(int i=0;i<Z;i++) 
	{
	    int locn = neighborLocs[i]; 
	    if (label[locn] == notlabel)  //don't include! but increment surface
		(*surface)++;
	    else if ((spin[locn] & NOTWHITE) == 0) // if "white"
	    {
		spin[locn] |= GRAY; Queue.push(locn);
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
    delete [] coords;

 
    return volume;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////





//////////////////////////////////////////////////////////////////////////////////////////////////////
// The following function is overloaded from 
// int Window::LY_domain_with_notlabel(int loc, int notlabel, 
//				       int *surface, int *boxvol, int *maxdim)
// The only difference is that it will NOT calculate the maxdim and boxvol, since for a single avalanche,
// its maxdim and boxvol have already been obtained in the function: LY_labelCC_getQ(...)
//////////////////////////////////////////////////////////////////////////////////////////////////////
int Window::LY_domain_with_notlabel(int loc, int notlabel, int *surface)
{
    int headQ;                     // head (location) of the queue
    queue<int> Queue;           

    *surface = 0;
    int volume = 0;
   
    int* startingcoords = new int[D];
    GetCoords(loc, startingcoords);   

    spin[loc] |= GRAY;	             // "gray"

    Queue.push(loc);              // push it into the real loc list   

    int* coords = new int[D];   
 
    while (!Queue.empty())     
    {
	headQ = Queue.front();  // get the real loc from the list
	Queue.pop();         // pop it from the list

	GetCoords(headQ, coords); 	// get its coordinates stored in coords[i] 
	GetNeighbors(headQ, coords);  	//  get its real neighbors' locs stored in neighborLocs[i] 

	// now look at nn: if white, not labeled notlabel, then color gray and queue
	for(int i=0;i<Z;i++) 
	{
	    int locn = neighborLocs[i]; 
	    if (label[locn] == notlabel)  //don't include! but increment surface
		(*surface)++;
	    else if ((spin[locn] & NOTWHITE) == 0) // if "white"
	    {
		spin[locn] |= GRAY; Queue.push(locn);
	    }
	}// end of for each neighbor loop


	spin[headQ] -= GRAY; // better have this bit set! (should by algorithm)
	spin[headQ] |= BLACK;
	volume++;
    }

    delete [] startingcoords;
    delete [] coords;

 
    return volume;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////



//////////////////////////////////////////////////////////////////////////////////////////////////////
// The following function is overloaded from 
// int Window::LY_domain_with_notlabel(int loc, int notlabel, int *surface)
// The only difference is that it will NOT calculate the surface.
// This is for the spin-DOWN domain.
//////////////////////////////////////////////////////////////////////////////////////////////////////
int Window::LY_domain_with_notlabel(int loc, int notlabel)
{
    int headQ;                     // head (location) of the queue
    queue<int> Queue;           

    int volume = 0;

    spin[loc] |= GRAY;	             // "gray"

    Queue.push(loc);              // push it into the real loc list   

    while (!Queue.empty())     
    {
	headQ = Queue.front();  // get the real loc from the list
	Queue.pop();         // pop it from the list

	GetNeighbors(headQ);  	//  get its real neighbors' locs stored in neighborLocs[i] 

	// now look at nn: if white, not labeled notlabel, then color gray and queue
	for(int i=0;i<Z;i++) 
	{
	    int locn = neighborLocs[i]; 
	    if (label[locn] == notlabel)  //don't include! but increment surface
		continue; //(*surface)++;
	    else if ((spin[locn] & NOTWHITE) == 0) // if "white"
	    {
		spin[locn] |= GRAY; Queue.push(locn);
	    }
	}// end of for each neighbor loop


	spin[headQ] -= GRAY; // better have this bit set! (should by algorithm)
	spin[headQ] |= BLACK;
	volume++;
    }

    return volume;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////




