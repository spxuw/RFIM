void Window::LY_reportD_with_notlabel(int notlabel)  
{
    int volume, surface, boxvol, maxdim;
    int nums = N;
	
    int domain_count = 1; 
		
    for (int start = 0; start < nums; start++) 
    {	
        // skip if already visited 
	if ((label[start] == notlabel) || (spin[start] & BLACK)) 
	    continue;
			
	// init statistics of this droplet.
	surface = 0;
	volume = LY_domain_with_notlabel(start, notlabel, &surface, &boxvol, &maxdim);

	cout << domain_count++ << ":  (" 
	    //<< (int)(spin[start]!=BLACK) << ") "
	     << (spin[start] & 1) << ")  "    
	     << volume << "/" << surface << "/" 
	     << boxvol << "/" << maxdim << "\n";
	
    }
  
    highbitzerotwos();   // all 4th bits are now zeroed
}

void Window::LY_reportD_with_notlabel(int depth, int notlabel, int* v_D, int* s_D, int* bv_D, int* md_D)  
{
    int volume, surface, boxvol, maxdim;
    int nums = N;
	
    //int domain_count = 1; 
		
    for (int start = 0; start < nums; start++) 
    {	
        // skip if already visited 
	if ((label[start] == notlabel) || (spin[start] & BLACK)) 
	    continue;
			
	// init statistics of this droplet.
	surface = 0;
	if(spin[start]&1)
	    volume = LY_domain_with_notlabel(start, notlabel, &surface, &boxvol, &maxdim);
	else
	    volume = LY_domain_with_notlabel(start, notlabel);

	// if depth=0 and spin-UP, this is the only one avalanche!!
	if ( (depth==0) && (spin[start] & 1) ) {
	    *v_D = volume;   *s_D = surface;	    
	    *bv_D= boxvol;  *md_D= maxdim;
	}


    }
  
    highbitzerotwos();   // all 4th bits are now zeroed
}

void Window::LY_reportD_with_notlabel(int depth, int notlabel, int* v_D, int* a_D)  
{
    int volume, surface;
    int nums = N;
	
    //int domain_count = 1; 
		
    for (int start = 0; start < nums; start++) 
    {	
        // skip if already visited 
	if ((label[start] == notlabel) || (spin[start] & BLACK)) 
	    continue;
			
	// init statistics of this droplet.
	surface = 0;
	if(spin[start]&1)
	    volume = LY_domain_with_notlabel(start, notlabel, &surface);
	else
	    volume = LY_domain_with_notlabel(start, notlabel);


	// if depth=0 and spin-UP, this is the only one avalanche!!
	if ( (depth==0) && (spin[start] & 1) ) {
	    *v_D = volume;   *a_D = surface;
	}

	
    }
  
    highbitzerotwos();   // all 4th bits are now zeroed
}

void Window::LY_reportD_with_notlabel(int depth, int notlabel, char* filename)  
{
    ofstream outfile(filename, ios_base::app); // data will be appended to this file 
    if(depth==0)
	outfile <<"//     v_D        a_D       bv_D       md_D      depth \n";

    int volume, surface, boxvol, maxdim;
    int nums = N;
	
    int domain_count = 1; 
		
    for (int start = 0; start < nums; start++) 
    {	
        // skip if already visited 
	if ((label[start] == notlabel) || (spin[start] & BLACK)) 
	    continue;
			
	// init statistics of this droplet.
	surface = 0;
	if(spin[start]&1)
	    volume = LY_domain_with_notlabel(start, notlabel, &surface, &boxvol, &maxdim);
	else
	    volume = LY_domain_with_notlabel(start, notlabel);


	if(spin[start]&1)
	{
	    outfile.width(10); outfile << volume          << " ";
	    outfile.width(10); outfile << surface         << " ";
	    outfile.width(10); outfile << boxvol          << " ";
	    outfile.width(10); outfile << maxdim          << " ";
	    outfile.width(10); outfile << depth           << " ";
	    outfile << endl;
	}

	domain_count++;
    }
  
    highbitzerotwos();   // all 4th bits are now zeroed

    outfile.close();

}

void Window::LY_reportD_with_notlabel(int depth, int notlabel, bool mark, int seed, char* filename, double H)  
{
    ofstream outfile(filename, ios_base::app); // append mode
    if(mark)
	if(depth==0)
	    outfile <<"/ seed= "<< seed << " v_D        a_D       bv_D       md_D      depth     H\n";

    int volume, surface, boxvol, maxdim;
    int nums = N;
	
    int domain_count = 1; 
		
    for (int start = 0; start < nums; start++) 
    {	
        // skip if already visited 
	if ((label[start] == notlabel) || (spin[start] & BLACK)) 
	    continue;
			
	// init statistics of this droplet.
	surface = 0;

	if(spin[start]&1)
	    volume = LY_domain_with_notlabel(start, notlabel, &surface, &boxvol, &maxdim);
	else
	    volume = LY_domain_with_notlabel(start, notlabel, &surface);

	//if((spin[start]&1) && (volume >3)) //Y.L. 03/06/06
	if(spin[start] & 1)  //Y.L. 06/16/06
	{
	    outfile.width(10); outfile << volume          << " ";
	    outfile.width(10); outfile << surface         << " ";
	    outfile.width(10); outfile << boxvol          << " ";
	    outfile.width(10); outfile << maxdim          << " ";
	    outfile.width(10); outfile << depth           << " ";
	    outfile.width(10); outfile << H               << endl;
	}

	domain_count++;
    }
  
    highbitzerotwos();   // all 4th bits are now zeroed

    outfile.close();

}
//////////////////////////////////////////////////////////////////////////////////////////////////////
