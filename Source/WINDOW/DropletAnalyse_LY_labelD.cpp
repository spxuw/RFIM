







int Window::LY_domain_with_notlabel(int loc, int notlabel, int *surface, int *boxvol, int *maxdim)
{
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

    *surface = 0;
    int volume = 0;
   
    int* startingcoords = new int[D];
    GetCoords(loc, startingcoords);   

    spin[loc] |= GRAY;	             

    Queue.push(loc);              

    int* coords = new int[D];   
 
    while (!Queue.empty())     
    {
	headQ = Queue.front();  
	Queue.pop();         

	GetCoords(headQ, coords); 	
	GetNeighbors(headQ, coords);  	

	
	for(int i=D-1;i>=0;i--)
	{
	    int temp = ((coords[i] - startingcoords[i] + size[i] + size[i]/2) % size[i]); 
	    
	    if(mark_del[i][temp]==0) 
	    {
		max_dim[i]++; 
		mark_del[i][temp]=1;
	    }
	}

	
	for(int i=0;i<Z;i++) 
	{
	    int locn = neighborLocs[i]; 
	    if (label[locn] == notlabel)  
		(*surface)++;
	    else if ((spin[locn] & NOTWHITE) == 0) 
	    {
		spin[locn] |= GRAY; Queue.push(locn);
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
    delete [] coords;

 
    return volume;
}













int Window::LY_domain_with_notlabel(int loc, int notlabel, int *surface)
{
    int headQ;                     
    queue<int> Queue;           

    *surface = 0;
    int volume = 0;
   
    int* startingcoords = new int[D];
    GetCoords(loc, startingcoords);   

    spin[loc] |= GRAY;	             

    Queue.push(loc);              

    int* coords = new int[D];   
 
    while (!Queue.empty())     
    {
	headQ = Queue.front();  
	Queue.pop();         

	GetCoords(headQ, coords); 	
	GetNeighbors(headQ, coords);  	

	
	for(int i=0;i<Z;i++) 
	{
	    int locn = neighborLocs[i]; 
	    if (label[locn] == notlabel)  
		(*surface)++;
	    else if ((spin[locn] & NOTWHITE) == 0) 
	    {
		spin[locn] |= GRAY; Queue.push(locn);
	    }
	}


	spin[headQ] -= GRAY; 
	spin[headQ] |= BLACK;
	volume++;
    }

    delete [] startingcoords;
    delete [] coords;

 
    return volume;
}










int Window::LY_domain_with_notlabel(int loc, int notlabel)
{
    int headQ;                     
    queue<int> Queue;           

    int volume = 0;

    spin[loc] |= GRAY;	             

    Queue.push(loc);              

    while (!Queue.empty())     
    {
	headQ = Queue.front();  
	Queue.pop();         

	GetNeighbors(headQ);  	

	
	for(int i=0;i<Z;i++) 
	{
	    int locn = neighborLocs[i]; 
	    if (label[locn] == notlabel)  
		continue; 
	    else if ((spin[locn] & NOTWHITE) == 0) 
	    {
		spin[locn] |= GRAY; Queue.push(locn);
	    }
	}


	spin[headQ] -= GRAY; 
	spin[headQ] |= BLACK;
	volume++;
    }

    return volume;
}





