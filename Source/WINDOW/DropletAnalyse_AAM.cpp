/****************************************************************************************************
Y.L. 02/28/06
The following is AAM's original code with the bug of maxdim removed by Y.L.

2 functions about connected clusters: 
   int Window::label_connected_cluster(int idx, int tolabel, char s, 
                                       int* surface, int* boxvol, int* maxdim)
   void Window::report_connected_clusters(int* bigidx, int* bigvol)

2 functions about domains: 
   int Window::domain_with_notlabel(int idx, int notlabel, 
	    			 int *surface, int *boxvol, int *maxdim)
   void Window::report_domains_with(int notlabel)

The last function gives a spin configuration's information of all connected clusters and domains.
   void Window::print_bubble_surfaces() 
*****************************************************************************************************/


//////////////////////////////////////////////////////////////////////////////////////////////////////
// return volume of the connected cluster. the starting point has spin s, index idx, and 
// finally all spins in this connected cluster will be labelled with tolabel!
// the surfacearea, boxvolume and maxdim will also be passed back via pointers. Y.L.02/27/06
//////////////////////////////////////////////////////////////////////////////////////////////////////
int Window::label_connected_cluster(int idx, int tolabel, char s, 
                                    int* surface, int* boxvol, int* maxdim)
{
    int idxn;		// index of spin, nbr
    int i, j, k, in, jn, kn;
    int del_i, del_j, del_k;
    int maxdim_i=0, maxdim_j=0, maxdim_k=0;
    vector<int> mark_del_i(W,0),mark_del_j(L,0),mark_del_k(H,0); 
    int istart, jstart, kstart;
    int headQ;
    queue<int> Queue;

    *surface = 0;
    int volume = 0;
    
    // get coordinates directly
    kstart = idx / (W * L);
    jstart = (idx - kstart * W * L)/W;
    istart = idx % W;

    // assume all vertices are "white" - have 0 in bits 4 and 5
    spin[idx] |= GRAY;		// GRAY the staring spin

    Queue.push(idx);
    label[idx] = tolabel;


    while (!Queue.empty()) 
    {
	headQ = Queue.front(); 
	Queue.pop();

	/********** for debug *************
	cout << "headQ=" << headQ << ","; 
	 ********* for debug *************/

	k = kn = headQ / L / W;
	j = jn = (headQ / W) % L;
	i = in = headQ % W;
	//   calculate the distance between i and istart, j and jstart, k and kstart. ////	
        //   pay attention to the part of getting the maxdim                    //Y.L.02/25/06
	del_i = ((i - istart + W + W/2) % W) - W/2; 
	del_j = ((j - jstart + L + L/2) % L) - L/2; 
	del_k = ((k - kstart + H + H/2) % H) - H/2; 
 	if(mark_del_i[del_i+W/2]==0) {maxdim_i++; mark_del_i[del_i+W/2]=1;}
 	if(mark_del_j[del_j+L/2]==0) {maxdim_j++; mark_del_j[del_j+L/2]=1;}
 	if(mark_del_k[del_k+H/2]==0) {maxdim_k++; mark_del_k[del_k+H/2]=1;}

	/********** for debug *************
	cout << "\n in=" << in << " jn=" << jn << " kn=" << kn << endl; 
	cout << "\n istart=" << istart << " jstart=" << jstart << " kstart=" << kstart << endl; 
	********** for debug *************/

	// now look at nearest neighbor: if white, same spin, then color gray and queue
	// Y.L. comment: this part is written wordily. 02/25/06
        // Y.L. So I write another concise fucntion :
        // int Window::label_connected_cluster_LY(int idx, int tolabel, char s, 
       //                                         int* surface, int* boxvol, int* maxdim)

	//////////////////////////      i     //////////////////////////////////
	in = (i + 1) % W;
	idxn = in + W * (jn + L * kn);
	if ((spin[idxn] & 1) != s)  //different spin, don't include! but increment surface
	    (*surface)++;
	else if((spin[idxn] & NOTWHITE) == 0) // if "white"
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

	//////////////////////////      j     //////////////////////////////////
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
    
  	//////////////////////////      k     //////////////////////////////////
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
    
	spin[headQ] -= GRAY; // better have this bit set! (should by algorithm)
	spin[headQ] |= BLACK;// BLACK this spin 
	volume++;
    }
    
    *maxdim = maxdim_i;
    if (maxdim_j > maxdim_i) *maxdim = maxdim_j;
    if (maxdim_k > maxdim_j) *maxdim = maxdim_k;
    *boxvol = maxdim_i * maxdim_j * maxdim_k;  

    /****************** debug
    cout << " maxdim_i=" << maxdim_i 
	 << " maxdim_j=" << maxdim_j 
	 << " maxdim_k=" << maxdim_k << endl; 

    cout << " volume=" << volume 
	 << " maxdim=" << *maxdim 
	 << " boxvol=" << *boxvol 
	 << endl;
     debug***************************/ 

    return volume;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////////////////////////////
// return volume of the domain
// this function is much similar with the above one: 
// int Window::label_connected_cluster(int idx, int tolabel, char s, 
//                                    int* surface, int* boxvol, int* maxdim)
// But here, we just concern those spins with index label are "notlabel"
// we don't care whether it is spin up or spin down
// so in this way we get the domain information.
// the surfacearea, boxvolume and maxdim will be passed back via pointers. Y.L.02/27/06
//////////////////////////////////////////////////////////////////////////////////////////////////////
int Window::domain_with_notlabel(int idx, int notlabel, 
				       int *surface, int *boxvol, int *maxdim)
{
    int idxn;		// index of spin, nbr
    int i, j, k, in, jn, kn;
    int del_i, del_j, del_k;
    int maxdim_i=0,maxdim_j=0, maxdim_k=0;
    vector<int> mark_del_i(W,0),mark_del_j(L,0),mark_del_k(H,0); 
    int istart, jstart, kstart;
    int headQ;
    queue <int> Queue;

    *surface = 0;
    int volume = 0;

    // get coordinates directly
    kstart = idx / (W * L);
    jstart = (idx - kstart * W * L)/W;
    istart = idx % W;

    // assume all vertices are "white" - have 0 in bits 4 and 5
    spin[idx] |= GRAY;		// "gray"

    Queue.push(idx);


    while (!Queue.empty()) 
    {
	headQ = Queue.front(); 
	Queue.pop();

	k = kn = headQ / L / W;
	j = jn = (headQ / W) % L;
	i = in = headQ % W;
	///##########   calculate the distance between i and istart, j and jstart, k and kstart.  #######////
	del_i = ((i - istart + W + W/2) % W) - W/2; 
	del_j = ((j - jstart + L + L/2) % L) - L/2; 
	del_k = ((k - kstart + H + H/2) % H) - H/2; 

 	if(mark_del_i[del_i+W/2]==0) {maxdim_i++; mark_del_i[del_i+W/2]=1;}
 	if(mark_del_j[del_j+L/2]==0) {maxdim_j++; mark_del_j[del_j+L/2]=1;}
 	if(mark_del_k[del_k+H/2]==0) {maxdim_k++; mark_del_k[del_k+H/2]=1;}

	// now look at nn: if white, not labeled notlabel, then color gray and queue
	//////////////////////////      i     //////////////////////////////////
	in = (i + 1) % W;
	idxn = in + W * (jn + L * kn);
	if (label[idxn] == notlabel)  //don't include! but increment surface
	    (*surface)++;
	else if ((spin[idxn] & NOTWHITE) == 0) // if "white"
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
  
	//////////////////////////      j     //////////////////////////////////
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

        //////////////////////////      k     //////////////////////////////////    
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
    
	spin[headQ] -= GRAY; // better have this bit set! (should by algorithm)
	spin[headQ] |= BLACK;// BLACK this spin 
	volume++;
    }

    *maxdim = maxdim_i;
    if (maxdim_j > maxdim_i) *maxdim = maxdim_j;
    if (maxdim_k > maxdim_j) *maxdim = maxdim_k;
    *boxvol = maxdim_i * maxdim_j * maxdim_k;  

    return volume;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////



//////////////////////////////////////////////////////////////////////////////////////////////////////
// determine and label all connected spin clusters. [Trace out droplet surfaces] 
// pass back index and volume of the largest connected cluster
//void Window::label_connected_clusters(int* bigidx, int* bigvol)
//////////////////////////////////////////////////////////////////////////////////////////////////////
void Window::report_connected_clusters(int* bigidx, int* bigvol)
{  
    int volume, surface, boxvol, maxdim;
    int start, numstart = W * L * H;
  
    int cluster_count = 1; // by Liu Yang 07/21/04
  
    char s;

    *bigidx = -1;
    *bigvol = 0;

    for (start = 0; start < numstart; ++start) 
    {
	// skip if already visited 
	if (spin[start] & BLACK) continue;
	
        // init statistics of this droplet. 
	s = spin[start];
	surface = 0;
	volume = label_connected_cluster(start, start, s, &surface, &boxvol, &maxdim);
	if (volume > *bigvol) {*bigvol = volume; *bigidx = start;}
	cout << cluster_count++ << ":  " <<  volume << "/" << surface << "/" << boxvol << "/" << maxdim << "\n";
    }
    cout << "\n";

    highbitzerotwos();   // all 4th bits are now zeroed

}
//////////////////////////////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////////////////////////////
// find all the connected clusters with index "notlabel", i.e. those clusters which are  
// touching with the cluster with index "label".
// here, "notlabel" means labels which are different from "label"
//
// determine their volume and surface area - report these 
// in this way, we get the domain volume and domain surface in this depth
//void report_connected_clusters_not(int notlabel);
//////////////////////////////////////////////////////////////////////////////////////////////////////
void Window::report_domains_with_notlabel(int notlabel)
{
    int volume, surface, boxvol, maxdim;
    int numstart = W * L * H;
	
    int domain_count = 1; // by Liu Yang 07/21/04
//    ofstream outf("test_output_Domains",  ios_base::out);

    for (int start = 0; start < numstart; ++start) 
    {
	// skip if the spin has label == "notlabel" or already visited 
	if ((label[start] == notlabel) || (spin[start] & BLACK)) 
	    continue;
			
	// init statistics of this droplet.
	surface = 0;
	//volume = notlabel_connected_cluster(start, notlabel, &surface, &boxvol, &maxdim);
	volume = domain_with_notlabel(start, notlabel, &surface, &boxvol, &maxdim);
		
	cout << domain_count++ << ":  " << volume << "/" << surface << "/" << boxvol << "/" << maxdim << "\n";
//	outf << domain_count++ << " " << volume << " " << surface << " " << boxvol << " " << maxdim << "\n";
    }
    
    highbitzerotwos(); // all 4th bits are now zeroed
}
//////////////////////////////////////////////////////////////////////////////////////////////////////






/*////////////////////////////////////////////////////////////////////////////////////////////////////
statistics on domains of connected spins, but here concerned about nesting, rather than clusters alone.
     Basic idea: 
    1, Determine and label all connected spin clusters (trace out droplet surfaces)
    2, Start from the largest cluster with index Q and volume Qvol. Relabel it with index -1.
    3, While  Qvol(number of spins with index -1)  is less than W*L*H (the system volume) 
      1)  find all the clusters with indices are NOT -1. i.e. those clusters are  
           touching with the cluster of index -1, but they could have holes
          (different spin values) inside themselves, in this sense, they are
          NOT connected clusters but domains.
      2)  determine their volume and surface area  
           ---- report these (in this way, we get the domain volume and domain surface in this depth)
      3)  in this depth, find all the connected clusters (same spin value)
          which are touching with the cluster of index -1.
          add their volume to Qvol (though not same spin value with the cluster of index
          -1) and relabel those connected clusters with index -1, 
      4)  depth++
      5)  do 1) until Qvol==W*L*H, i.e. until we get the deepest level
////////////////////////////////////////////////////////////////////////////////////////////////////*/
void Window::print_bubble_surfaces() 
{
    int start, dir, nbr, numstart = W * L * H;
    int surf, vol, maxdim;

    if (label == NULL) label = new int[W * L * H];

    int Q;
    int Qvol;

    cout << "\n cluster_count: volume/surface/boxvol/maxdim \n";
    //label_connected_clusters(&Q, &Qvol);	// cluster is labelled with the origin spin
    report_connected_clusters(&Q, &Qvol);	// cluster is labelled with the origin spin

    // the largest connected cluster's label is now Q; but specialize: set all Q spins to have label -1
    for (start = 0; start < numstart; ++start)
	if (label[start] == Q)
	    label[start] = -1;

    int depth = 0;     // depth is how deep the droplets are nested.
    cout << "\n Domain_count: volume/surface/boxvol/maxdim \n";
  
    while (Qvol != W * L * H) 
    {       
	////////////////////////////////////////////////////////////////////////////////////////////
        // tell surface & volume of clusters (connected by not-1 labels)
        // this will give all the depth=0 domains, actually this is enough for our problem!!!!
        // because we don't care how deep the cluster is and how many subdomain ("holes") there are.
        // so. when inserted in the Hysteresis code, we can ignore the following part to get the 
        // subdomians.  Y.L.02/27/06
	//report_connected_clusters_not(-1);  
	report_domains_with_notlabel(-1);  
	////////////////////////////////////////////////////////////////////////////////////////////
    
	for (start = 0; start < numstart; ++start) 
	{
	    if (label[start] == -1)
		for (dir = 0; dir < 6; ++dir)
		{
		    nbr = getnbr(start, dir);
		    if (label[nbr] != -1 && label[nbr] != Q)
			// ^ don't relabel with Q if already done
			// have found a nbr spin to -1 cluster, if not done yet,
			// give this cluster the label of Q, return volume of cluster 
                        // to increment total volume associated with Q (though not same spin)
                        //  Cool! Y.L.02/27/06  
			Qvol += label_connected_cluster(nbr, Q, spin[nbr], &surf, &vol, &maxdim);
		        //Qvol will be equal to the system volume when we get to the deepest level
                        //of the domain
		}
	}
	// now take all of the Q clusters and relabel as -1
        // in this sense, we dig into a deeper level //Y.L.02/27/06
	for (start = 0; start < numstart; ++start)
	    if (label[start] == Q)
		label[start] = -1;
    
	highbitzerotwos();
	depth++;
    }

    cout << "\nDepth= " << depth << endl << endl;  //02/24/06 by Yang Liu
}
//////////////////////////////////////////////////////////////////////////////////////////////////////

