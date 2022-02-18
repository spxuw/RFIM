void Window::LY_print_bubble_surfaces()
{
    int start, locn;
    int nums = N;

    int surf, vol, maxdim;

    if (label == NULL) 
	label = new int[nums];

    int Q;
    int Qvol;

	
    cout << "\n cluster_count: (up/down) volume/surface/boxvol/maxdim/Asphericity/BoundaryTouching\n";
   
    LY_modified_reportCC(&Q, &Qvol);

    cout << "Q=" << Q << "  Qvol=" << Qvol << endl;


    for (start = 0; start < nums; ++start)
	if (label[start] == Q)
	    label[start] = -1;


    int depth = 0 ;
    cout << "\n Depth:Domain_count: (up/down) volume/surface/boxvol/maxdim \n";

    while (Qvol != nums) 
    {
	LY_reportD_with_notlabel(-1); 

	for (start = 0; start < nums; ++start) {
	    if (label[start] == -1) {
		GetNeighbors(start);
		for(int i=0; i<Z;i++){
		    locn = neighborLocs[i]; 
		    if (label[locn] != -1 && label[locn] != Q)
			Qvol += LY_labelCC(locn, Q, spin[locn], &surf, &vol, &maxdim);
		}
	    }
	}
	// now take all of the Q clusters and relabel as -1
	for (start = 0; start < nums; ++start)
	    if (label[start] == Q)
		label[start] = -1;
    
	highbitzerotwos();
	depth++;

    }

    // depth is how deep the droplets are nested.
    cout << "\nDepth=" << depth <<"\n\n"; 


}

void Window::LY_analyse_one_avalanche()
{
    int start, locn;
    int surf, vol, maxdim;
    int nums = N;
    int Q;
    int Qvol;

    if (label == NULL) 
	label = new int[nums];

    //////////////////    Connected-Cluster description of the avalanche      //////////////////////////
    int v_CC;
    int a_CC;
    int bv_CC;
    int md_CC;
    Vec_DP EV_CC;

    cout << "\n cluster_count: (up/down) volume/surface/boxvol/maxdim/Asphericity/BoundaryTouching\n";
   
    LY_modified_reportCC(&Q, &Qvol, &v_CC, &a_CC, &bv_CC, &md_CC, &EV_CC);

    // Largest spin-DOWN boundarytouching connected cluster label is now Q .Y.L. 02/28/06
    // but specialize: set all Q spins to have label -1
    for (start = 0; start < nums; ++start)
	if (label[start] == Q)
	    label[start] = -1;


    /////////////////     Domain description of the avalanche       ////////////////////////////////////

    int v_D;
    int a_D;
    //int bv_D;
    //int md_D;
    int depth = 0 ;
    cout << "\n Depth:Domain_count: (up/down) volume/surface \n";

    while (Qvol != nums) 
    {
	LY_reportD_with_notlabel(depth, -1, &v_D, &a_D);  
	//LY_reportD_with_notlabel(depth, -1, &v_D, &a_D, &bv_D, &md_D);  

	for (start = 0; start < nums; ++start) {
	    if (label[start] == -1) {
		GetNeighbors(start);
		for(int i=0; i<Z;i++){
		    locn = neighborLocs[i]; 
		    if (label[locn] != -1 && label[locn] != Q)
			Qvol += LY_labelCC(locn, Q, spin[locn], &surf, &vol, &maxdim);
		}
	    }
	}
	// now take all of the Q clusters and relabel as -1
	for (start = 0; start < nums; ++start)
	    if (label[start] == Q)
		label[start] = -1;
    
	highbitzerotwos();
	depth++;
    }

    // depth is how deep the droplets are nested.
    cout << "\nDepth=" << depth <<"\n\n"; 

    // show result
    cout << "\n\n This Avalanche has:"
 	 << "\n v_CC= " << v_CC << " a_CC= " << a_CC << " bv_CC= " << bv_CC << " md_CC= " << md_CC 
	 << "\n v_D= " << v_D << " a_D= " << a_D << endl;
    cout << "Shape tensor information:" << endl;
    for(int i=0;i<D+8;i++)
	cout <<"EVinfo["<<i<<"]= " << EV_CC[i] << endl;

}

void Window::LY_analyse_one_avalanche(char* filename)
{
    ofstream outfile(filename, ios_base::out);
    outfile << "//v_D a_D v_CC a_CC bv_CC md_CC and  A1 A2 Asph. Prol. tr_P2 trQ_2 tr_P3 trQ_3 ev[0]...ev[D-1]\n";

    int start, locn;
    int surf, vol, maxdim;
    int nums = N;
    int Q;
    int Qvol;

    if (label == NULL) 
	label = new int[nums];

    //////////////////    Connected-Cluster description of the avalanche      //////////////////////////
    int v_CC;
    int a_CC;
    int bv_CC;
    int md_CC;
    Vec_DP EV_CC;

    cout << "\n cluster_count: (up/down) volume/surface/boxvol/maxdim/Asphericity/BoundaryTouching\n";
   
    LY_modified_reportCC(&Q, &Qvol, &v_CC, &a_CC, &bv_CC, &md_CC, &EV_CC);

    // Largest spin-DOWN boundarytouching connected cluster label is now Q .Y.L. 02/28/06
    // but specialize: set all Q spins to have label -1
    for (start = 0; start < nums; ++start)
	if (label[start] == Q)
	    label[start] = -1;


    /////////////////     Domain description of the avalanche       ////////////////////////////////////

    int v_D;
    int a_D;
    //int bv_D;
    //int md_D;
    int depth = 0 ;
    cout << "\n Depth:Domain_count: (up/down) volume/surface \n";

    while (Qvol != nums) 
    {
	LY_reportD_with_notlabel(depth, -1, &v_D, &a_D);  


	for (start = 0; start < nums; ++start) {
	    if (label[start] == -1) {
		GetNeighbors(start);
		for(int i=0; i<Z;i++){
		    locn = neighborLocs[i]; 
		    if (label[locn] != -1 && label[locn] != Q)
			Qvol += LY_labelCC(locn, Q, spin[locn], &surf, &vol, &maxdim);
		}
	    }
	}
	// now take all of the Q clusters and relabel as -1
	for (start = 0; start < nums; ++start)
	    if (label[start] == Q)
		label[start] = -1;
    
	highbitzerotwos();
	depth++;
    }

    // depth is how deep the droplets are nested.
    cout << "\nDepth=" << depth; 

    // show result
    cout << "\n This Avalanche has:"
 	 << "\n v_CC= " << v_CC << " a_CC= " << a_CC << " bv_CC= " << bv_CC << " md_CC= " << md_CC 
	 << "\n v_D= " << v_D << " a_D= " << a_D << endl;
    cout << "Shape tensor information:" << endl;
    for(int i=0;i<D+8;i++)
	cout << EV_CC[i] << " ";
    cout << endl;

    outfile << v_D    << " ";
    outfile << a_D    << " ";
    outfile << v_CC   << " ";
    outfile << a_CC   << " ";
    outfile << bv_CC  << " ";
    outfile << md_CC  << " ";
    for(int i=0;i<D+8; i++)
   	outfile << EV_CC[i] << " ";
    outfile << endl;
    outfile.close();


}

int Window::LY_analyse_many_avalanches(char* mainname)
{
    char fname1[256];
    char fname2[256];
    sprintf(fname1,"%s_CC",mainname);
    sprintf(fname2,"%s_D",mainname);

    int start, locn;
    int nums = N;

    int surf, vol, maxdim;

    if (label == NULL) 
	label = new int[nums];

    int Q;
    int Qvol;

    cout << "\n cluster_count: (up/down) volume/surface/boxvol/maxdim/Asphericity/BoundaryTouching\n";
   
    int UPcluster_count = LY_modified_reportCC(&Q, &Qvol, fname1);

    // Largest spin-DOWN boundarytouching connected cluster label is now Q .Y.L. 02/28/06
    // but specialize: set all Q spins to have label -1
    for (start = 0; start < nums; ++start)
	if (label[start] == Q)
	    label[start] = -1;

    int depth = 0 ;
    cout << "\n Depth:Domain_count: (up/down) volume/surface/boxvol/maxdim \n";

    while (Qvol != nums) 
    {
	LY_reportD_with_notlabel(depth, -1, fname2); 

	for (start = 0; start < nums; ++start) {
	    if (label[start] == -1) {
		GetNeighbors(start);
		for(int i=0; i<Z;i++){
		    locn = neighborLocs[i]; 
		    if (label[locn] != -1 && label[locn] != Q)
			Qvol += LY_labelCC(locn, Q, spin[locn],&surf, &vol, &maxdim);
		}
	    }
	}
	// now take all of the Q clusters and relabel as -1
	for (start = 0; start < nums; ++start)
	    if (label[start] == Q)
		label[start] = -1;
    
	highbitzerotwos();
	depth++;

    }

    // depth is how deep the droplets are nested or the number of levels of nesting.
    cout << "\nDepth=" << depth <<"\n\n"; 

    //return depth;
    return UPcluster_count;
}

void Window::LY_analyse_one_avalanche(bool mark, int seed, char* filename)
{
    ofstream outfile(filename, ios_base::app); // append mode
    if(mark)
	outfile << "/ seed= "<< seed<<" v_D a_D v_CC a_CC bv_CC md_CC and  A1 A2 Asph. Prol. tr_P2 trQ_2 tr_P3 trQ_3 ev[0]...ev[D-1]\n";

    int start;//, locn;
    //int surf, vol, maxdim;
    int nums = N;
    int Q;
    int Qvol;

    if (label == NULL) 
	label = new int[nums];

    //////////////////    Connected-Cluster description of the avalanche      //////////////////////////
    int v_CC;
    int a_CC;
    int bv_CC;
    int md_CC;
    Vec_DP EV_CC;

    LY_modified_reportCC(&Q, &Qvol, &v_CC, &a_CC, &bv_CC, &md_CC, &EV_CC);

    // Largest spin-DOWN boundarytouching connected cluster label is now Q .Y.L. 02/28/06
    // but specialize: set all Q spins to have label -1
    for (start = 0; start < nums; ++start)
	if (label[start] == Q)
	    label[start] = -1;


    /////////////////     Domain description of the avalanche       ////////////////////////////////////
    int v_D;
    int a_D;
    //int bv_D;
    //int md_D;
    int depth = 0 ;


    LY_reportD_with_notlabel(depth, -1, &v_D, &a_D);  
 
    outfile << v_D    << " ";
    outfile << a_D    << " ";
    outfile << v_CC   << " ";
    outfile << a_CC   << " ";
    outfile << bv_CC  << " ";
    outfile << md_CC  << " ";
    for(int i=0;i<D+8; i++)
   	outfile << EV_CC[i] << " ";
    outfile << endl;
    outfile.close();


}

int Window::LY_modified_analyse_one_avalanche(bool &mark, int seed, char* filename)
{
    ofstream outfile(filename, ios_base::app); // append mode
    if(mark)
    {
	outfile << "/ seed= "<< seed
		<< "v_D a_D v_CC a_CC bv_CC md_CC A1 A2 Asp. Pro. tr_P2 trQ_2 tr_P3 trQ_3 ev[0]...ev[D-1] Nhole\n";
	mark=false;
    }
    //int start, locn;
    //int surf, vol, maxdim;
    int nums = N;
    int Q;
    int Qvol;

    if (label == NULL) 
	label = new int[nums];

    //////////////////    Connected-Cluster description of the avalanche      //////////////////////////
    int v_CC;
    int a_CC;
    int bv_CC;
    int md_CC;
    Vec_DP EV_CC;
    vector<int> HoleInfo(3);

    LY_modified_reportCC(&Q, &Qvol, 
			 &v_CC, &a_CC, &bv_CC, &md_CC, 
			 &EV_CC, 
			 &HoleInfo);


    /////////////////     Domain description of the avalanche       ////////////////////////////////////
    int v_D;
    int a_D;
    int Nhole;
    v_D    = v_CC + HoleInfo[0];
    a_D    = a_CC - HoleInfo[1];
    Nhole = HoleInfo[2]; 


    ////////////////      write data to a file                   ///////////////////////////////////////
    outfile << v_D    << " ";
    outfile << a_D    << " ";
    outfile << v_CC   << " ";
    outfile << a_CC   << " ";
    outfile << bv_CC  << " ";
    outfile << md_CC  << " ";
    for(int i=0;i<D+8; i++)
   	outfile << EV_CC[i] << " ";
    outfile << Nhole;
    outfile << endl;
    outfile.close();

    return v_CC;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// the same as the above one, but also record the external field
int Window::LY_modified_analyse_one_avalanche(bool &mark, int seed, char* filename, double H)
{
    ofstream outfile(filename, ios::app); // append mode
    if(mark)
    {
 	outfile <<"/ seed= "<< seed
		<<" v_D a_D v_CC a_CC bv_CC md_CC A1 A2 Asp. Pro. tr_P2 trQ_2 tr_P3 trQ_3 ev[0]...ev[D-1] Nhole H\n";
	mark=false;
    }

    //int start, locn;
    //int surf, vol, maxdim;
    int nums = N;
    int Q;
    int Qvol;

    if (label == NULL) 
	label = new int[nums];

    //////////////////    Connected-Cluster description of the avalanche      //////////////////////////
    int v_CC;
    int a_CC;
    int bv_CC;
    int md_CC;
    Vec_DP EV_CC;
    vector<int> HoleInfo(3);

    LY_modified_reportCC(&Q, &Qvol, 
			 &v_CC, &a_CC, &bv_CC, &md_CC, 
			 &EV_CC, 
			 &HoleInfo);


    /////////////////     Domain description of the avalanche       ////////////////////////////////////
    int v_D;
    int a_D;
    int Nhole;
    v_D    = v_CC + HoleInfo[0];
    a_D    = a_CC - HoleInfo[1];
    Nhole = HoleInfo[2]; 


    ////////////////      write data to a file                   ///////////////////////////////////////
    outfile << v_D    << ' ';
    outfile << a_D    << ' ';
    outfile << v_CC   << ' ';
    outfile << a_CC   << ' ';
    outfile << bv_CC  << ' ';
    outfile << md_CC  << ' ';
    for(int i=0;i<D+8; i++)
   	outfile << EV_CC[i] << ' ';
    outfile << Nhole << ' ';
    outfile << H << endl;
    outfile.close();


    return v_CC;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void Window::LY_analyse_one_small_avalanche(bool &mark, int seed, int size, char* filename, double H)
{
    ofstream outfile(filename, ios_base::app); // append mode

    if(mark)
    {
	outfile <<"/ seed= "<< seed
		<<" v_D a_D v_CC a_CC bv_CC md_CC A1 A2 Asp. Pro. tr_P2 trQ_2 tr_P3 trQ_3 ev[0]...ev[D-1] Nhole H\n";
	mark=false;
    }

    if(size==1)
    {
	outfile << 1     << ' ';//v_D    << ' ';
	outfile << 2*D   << ' ';//a_D    << ' ';
	outfile << 1     << ' ';//v_CC   << ' ';
	outfile << 2*D   << ' ';//a_CC   << ' ';
	outfile << 1     << ' ';//bv_CC  << ' ';
	outfile << 1     << ' ';//md_CC  << ' ';
	for(int i=0;i<D+8; i++)
	    outfile << -1<< ' ';//EV_CC[i] << ' ';
    }// end of if size==1
    else if (size==2)
    {
	outfile << 2     << ' ';//v_D    << ' ';
	outfile << 4*D-2 << ' ';//a_D    << ' ';
	outfile << 2     << ' ';//v_CC   << ' ';
	outfile << 4*D-2 << ' ';//a_CC   << ' ';
	outfile << 2     << ' ';//bv_CC  << ' ';
	outfile << 2     << ' ';//md_CC  << ' ';

	if (D==2)
	{//0 -1 1 -1 0.03125 0.0625 0 0.015625 0.25 0 
	    outfile << 0          << ' ';//EV_CC[0]:Anisotropy1
	    outfile << -1         << ' ';//EV_CC[1]:Anisotropy2
	    outfile << 1          << ' ';//EV_CC[2]:Asphericity
	    outfile << -1         << ' ';//EV_CC[3]:Prolateness
	    outfile << 0.03125    << ' ';//EV_CC[4]:tr_P2
	    outfile << 0.0625     << ' ';//EV_CC[5]:trQ_2
	    outfile << 0          << ' ';//EV_CC[6]:tr_P3
	    outfile << 0.015625   << ' ';//EV_CC[7]:trQ_3
	    outfile << 0.25       << ' ';//EV_CC[8]:lambda[0]
	    outfile << 0          << ' ';//EV_CC[9]:lambda[1]
	}
	else if(D==3)
	{//0 0 1 1 0.0416667 0.0625 0.00347222 0.015625 0.25 0 0 
	    outfile << 0          << ' ';//EV_CC[0]:Anisotropy1
	    outfile << 0          << ' ';//EV_CC[1]:Anisotropy2
	    outfile << 1          << ' ';//EV_CC[2]:Asphericity
	    outfile << 1          << ' ';//EV_CC[3]:Prolateness
	    outfile << 0.0416667  << ' ';//EV_CC[4]:tr_P2
	    outfile << 0.0625     << ' ';//EV_CC[5]:trQ_2
	    outfile << 0.00347222 << ' ';//EV_CC[6]:tr_P3
	    outfile << 0.015625   << ' ';//EV_CC[7]:trQ_3
	    outfile << 0.25       << ' ';//EV_CC[8]:lambda[0]
	    outfile << 0          << ' ';//EV_CC[9]:lambda[1]
	    outfile << 0          << ' ';//EV_CC[10]:lambda[2]
	}
    }// end of if size==2

    outfile << 0 << ' '; //Nhole=0
    outfile << H << endl;
    outfile.close();
    
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////




//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int Window::LY_analyse_onespin_avalanche(bool &mark, int seed, char* filename, double H)
{
    ofstream outfile(filename, ios_base::app); // append mode
    if(mark)
    {
	outfile <<"/ seed= "<< seed
		<<" v_D a_D v_CC a_CC bv_CC md_CC A1 A2 Asp. Pro. tr_P2 trQ_2 tr_P3 trQ_3 ev[0]...ev[D-1] Nhole H\n";
	mark=false;
    }

    outfile << 1     << ' ';//v_D    << ' ';
    outfile << 2*D   << ' ';//a_D    << ' ';
    outfile << 1     << ' ';//v_CC   << ' ';
    outfile << 2*D   << ' ';//a_CC   << ' ';
    outfile << 1     << ' ';//bv_CC  << ' ';
    outfile << 1     << ' ';//md_CC  << ' ';
    for(int i=0;i<D+8; i++)
	outfile << -1 << ' ';//EV_CC[i] << ' ';

    outfile << 0 << ' '; //Nhole=0
    outfile << H << endl;
    outfile.close();
    
    return 1;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////




//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int Window::LY_analyse_many_avalanches(bool mark, int seed, char* mainname, double H)
{
    char fname1[256];
    char fname2[256];
    sprintf(fname1,"%s_CC",mainname);
    sprintf(fname2,"%s_D",mainname);

    int start, locn;
    int nums = N;

    int surf, vol, maxdim;

    if (label == NULL) 
	label = new int[nums];

    int Q;
    int Qvol;

    int UPcluster_count = LY_modified_reportCC(&Q, &Qvol, mark, seed, fname1, H);

    // Largest spin-DOWN boundarytouching connected cluster label is now Q .Y.L. 02/28/06
    // but specialize: set all Q spins to have label -1
    for (start = 0; start < nums; ++start)
	if (label[start] == Q)
	    label[start] = -1;

    int depth = 0 ;
    while (Qvol != nums) 
    {
	LY_reportD_with_notlabel(depth, -1, mark, seed, fname2, H); 

	for (start = 0; start < nums; ++start) {
	    if (label[start] == -1) {
		GetNeighbors(start);
		for(int i=0; i<Z;i++){
		    locn = neighborLocs[i]; 
		    if (label[locn] != -1 && label[locn] != Q)
			Qvol += LY_labelCC(locn, Q, spin[locn], &surf, &vol, &maxdim);
		}
	    }
	}
	// now take all of the Q clusters and relabel as -1
	for (start = 0; start < nums; ++start)
	    if (label[start] == Q)
		label[start] = -1;
    
	highbitzerotwos();
	depth++;

    }

    //return depth;
    return UPcluster_count;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////






