int Window::label_connected_cluster(int idx, int tolabel, char s, 
                                    int* surface, int* boxvol, int* maxdim)
{
    int idxn;		
    int i, j, k, in, jn, kn;
    int del_i, del_j, del_k;
    int maxdim_i=0, maxdim_j=0, maxdim_k=0;
    vector<int> mark_del_i(W,0),mark_del_j(L,0),mark_del_k(H,0); 
    int istart, jstart, kstart;
    int headQ;
    queue<int> Queue;

    *surface = 0;
    int volume = 0;
    
    
    kstart = idx / (W * L);
    jstart = (idx - kstart * W * L)/W;
    istart = idx % W;

    
    spin[idx] |= GRAY;		

    Queue.push(idx);
    label[idx] = tolabel;


    while (!Queue.empty()) 
    {
	headQ = Queue.front(); 
	Queue.pop();



	k = kn = headQ / L / W;
	j = jn = (headQ / W) % L;
	i = in = headQ % W;
	
        
	del_i = ((i - istart + W + W/2) % W) - W/2; 
	del_j = ((j - jstart + L + L/2) % L) - L/2; 
	del_k = ((k - kstart + H + H/2) % H) - H/2; 
 	if(mark_del_i[del_i+W/2]==0) {maxdim_i++; mark_del_i[del_i+W/2]=1;}
 	if(mark_del_j[del_j+L/2]==0) {maxdim_j++; mark_del_j[del_j+L/2]=1;}
 	if(mark_del_k[del_k+H/2]==0) {maxdim_k++; mark_del_k[del_k+H/2]=1;}



	
	in = (i + 1) % W;
	idxn = in + W * (jn + L * kn);
	if ((spin[idxn] & 1) != s)  
	    (*surface)++;
	else if((spin[idxn] & NOTWHITE) == 0) 
	{
	    spin[idxn] |= GRAY; label[idxn] = tolabel; Queue.push(idxn);
	}

	in = (i + W - 1) % W;
	idxn = in + W * (jn + L * kn);
	if ((spin[idxn] & 1) != s) 
	    (*surface)++;
        else if ((spin[idxn] & NOTWHITE) == 0)
	{
	    spin[idxn] |= GRAY; label[idxn] = tolabel; Queue.push(idxn);
	}
	in = i;

	
	jn = (j + L - 1) % L;
	idxn = in + W * (jn + L * kn);
	if ((spin[idxn] & 1)!=s) 
	    (*surface)++;
	else if ((spin[idxn] & NOTWHITE) == 0)
	{
	    spin[idxn] |= GRAY; label[idxn] = tolabel; Queue.push(idxn);
	}
	
	jn = (j + 1) % L;
	idxn = in + W * (jn + L * kn);
	if ((spin[idxn] & 1)!=s) 
	    (*surface)++;
	else if ((spin[idxn] & NOTWHITE) == 0)
	{
	    spin[idxn] |= GRAY; label[idxn] = tolabel; Queue.push(idxn);
	}
	jn = j;
    
  	
	kn = (k + H - 1) % H;
	idxn = in + W * (jn + L * kn);
	if ((spin[idxn] & 1)!=s) 
	    (*surface)++;
	else if((spin[idxn] & NOTWHITE) == 0)
	{
	    spin[idxn] |= GRAY; label[idxn] = tolabel; Queue.push(idxn);
	}
	
	kn = (k + 1) % H;
	idxn = in + W * (jn + L * kn);
	if((spin[idxn] & 1)!=s) 
	    (*surface)++;
	else if((spin[idxn] & NOTWHITE) == 0)
	{
	    spin[idxn] |= GRAY; label[idxn]=tolabel; Queue.push(idxn);
	}
	kn = k;
    
	spin[headQ] -= GRAY; 
	spin[headQ] |= BLACK;
	volume++;
    }
    
    *maxdim = maxdim_i;
    if (maxdim_j > maxdim_i) *maxdim = maxdim_j;
    if (maxdim_k > maxdim_j) *maxdim = maxdim_k;
    *boxvol = maxdim_i * maxdim_j * maxdim_k;  



    return volume;
}

int Window::domain_with_notlabel(int idx, int notlabel, 
				       int *surface, int *boxvol, int *maxdim)
{
    int idxn;		
    int i, j, k, in, jn, kn;
    int del_i, del_j, del_k;
    int maxdim_i=0,maxdim_j=0, maxdim_k=0;
    vector<int> mark_del_i(W,0),mark_del_j(L,0),mark_del_k(H,0); 
    int istart, jstart, kstart;
    int headQ;
    queue <int> Queue;

    *surface = 0;
    int volume = 0;

    
    kstart = idx / (W * L);
    jstart = (idx - kstart * W * L)/W;
    istart = idx % W;

    
    spin[idx] |= GRAY;		

    Queue.push(idx);


    while (!Queue.empty()) 
    {
	headQ = Queue.front(); 
	Queue.pop();

	k = kn = headQ / L / W;
	j = jn = (headQ / W) % L;
	i = in = headQ % W;
	
	del_i = ((i - istart + W + W/2) % W) - W/2; 
	del_j = ((j - jstart + L + L/2) % L) - L/2; 
	del_k = ((k - kstart + H + H/2) % H) - H/2; 

 	if(mark_del_i[del_i+W/2]==0) {maxdim_i++; mark_del_i[del_i+W/2]=1;}
 	if(mark_del_j[del_j+L/2]==0) {maxdim_j++; mark_del_j[del_j+L/2]=1;}
 	if(mark_del_k[del_k+H/2]==0) {maxdim_k++; mark_del_k[del_k+H/2]=1;}

	
	
	in = (i + 1) % W;
	idxn = in + W * (jn + L * kn);
	if (label[idxn] == notlabel)  
	    (*surface)++;
	else if ((spin[idxn] & NOTWHITE) == 0) 
	{
	    spin[idxn] |= GRAY; Queue.push(idxn);
	}
	
	in = (i + W - 1) % W;
	idxn = in + W * (jn + L * kn);
	if (label[idxn] == notlabel)
	    (*surface)++;
	else if ((spin[idxn] & NOTWHITE) == 0)
	{
	    spin[idxn] |= GRAY; Queue.push(idxn);
	}
	in = i;
  
	
	jn = (j + L - 1) % L;
	idxn = in + W * (jn + L * kn);
	if (label[idxn] == notlabel)
	    (*surface)++;
	else if ((spin[idxn] & NOTWHITE) == 0)
	{
	    spin[idxn] |= GRAY; Queue.push(idxn);
	}

	jn = (j + 1) % L;
	idxn = in + W * (jn + L * kn);
	if (label[idxn] == notlabel)
	    (*surface)++;
	else if ((spin[idxn] & NOTWHITE) == 0)
	{
	    spin[idxn] |= GRAY; Queue.push(idxn);
	}
	jn = j;

        
	kn = (k + H - 1) % H;
	idxn = in + W * (jn + L * kn);
	if (label[idxn] == notlabel) {
	    (*surface)++;
	} else
	    if ((spin[idxn] & NOTWHITE) == 0)
	    {spin[idxn] |= GRAY; Queue.push(idxn);}
	

	kn = (k + 1) % H;
	idxn = in + W * (jn + L * kn);
	if (label[idxn] == notlabel)
	    (*surface)++;
	else if ((spin[idxn] & NOTWHITE) == 0)
	{
	    spin[idxn] |= GRAY; Queue.push(idxn);
	}
	kn = k;
    
	spin[headQ] -= GRAY; 
	spin[headQ] |= BLACK;
	volume++;
    }

    *maxdim = maxdim_i;
    if (maxdim_j > maxdim_i) *maxdim = maxdim_j;
    if (maxdim_k > maxdim_j) *maxdim = maxdim_k;
    *boxvol = maxdim_i * maxdim_j * maxdim_k;  

    return volume;
}









void Window::report_connected_clusters(int* bigidx, int* bigvol)
{  
    int volume, surface, boxvol, maxdim;
    int start, numstart = W * L * H;
  
    int cluster_count = 1; 
  
    char s;

    *bigidx = -1;
    *bigvol = 0;

    for (start = 0; start < numstart; ++start) 
    {
	
	if (spin[start] & BLACK) continue;
	
        
	s = spin[start];
	surface = 0;
	volume = label_connected_cluster(start, start, s, &surface, &boxvol, &maxdim);
	if (volume > *bigvol) {*bigvol = volume; *bigidx = start;}
	cout << cluster_count++ << ":  " <<  volume << "/" << surface << "/" << boxvol << "/" << maxdim << "\n";
    }
    cout << "\n";

    highbitzerotwos();   

}












void Window::report_domains_with_notlabel(int notlabel)
{
    int volume, surface, boxvol, maxdim;
    int numstart = W * L * H;
	
    int domain_count = 1; 


    for (int start = 0; start < numstart; ++start) 
    {
	
	if ((label[start] == notlabel) || (spin[start] & BLACK)) 
	    continue;
			
	
	surface = 0;
	
	volume = domain_with_notlabel(start, notlabel, &surface, &boxvol, &maxdim);
		
	cout << domain_count++ << ":  " << volume << "/" << surface << "/" << boxvol << "/" << maxdim << "\n";

    }
    
    highbitzerotwos(); 
}








void Window::print_bubble_surfaces() 
{
    int start, dir, nbr, numstart = W * L * H;
    int surf, vol, maxdim;

    if (label == NULL) label = new int[W * L * H];

    int Q;
    int Qvol;

    cout << "\n cluster_count: volume/surface/boxvol/maxdim \n";
    
    report_connected_clusters(&Q, &Qvol);	

    
    for (start = 0; start < numstart; ++start)
	if (label[start] == Q)
	    label[start] = -1;

    int depth = 0;     
    cout << "\n Domain_count: volume/surface/boxvol/maxdim \n";
  
    while (Qvol != W * L * H) 
    {       
	
        
        
        
        
        
	
	report_domains_with_notlabel(-1);  
	
    
	for (start = 0; start < numstart; ++start) 
	{
	    if (label[start] == -1)
		for (dir = 0; dir < 6; ++dir)
		{
		    nbr = getnbr(start, dir);
		    if (label[nbr] != -1 && label[nbr] != Q)
			
			
			
                        
                        
			Qvol += label_connected_cluster(nbr, Q, spin[nbr], &surf, &vol, &maxdim);
		        
                        
		}
	}
	
        
	for (start = 0; start < numstart; ++start)
	    if (label[start] == Q)
		label[start] = -1;
    
	highbitzerotwos();
	depth++;
    }

    cout << "\nDepth= " << depth << endl << endl;  
}


