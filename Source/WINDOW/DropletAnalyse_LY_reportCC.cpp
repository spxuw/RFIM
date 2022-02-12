//////////////////////////////////////////////////////////////////////////////////////////////////////
//  The following function is almost the same as the function 
//    void Window::report_connected_clusters(int* bigidx, int* bigvol)
// The difference is that it use the 
// LY_labelCC_getQ(...) instead of label_connected_cluster(...)
// so it can give the shape tensor information. 
//  It will label domains of connected spins, pass back index and volume of largest cluster.
//////////////////////////////////////////////////////////////////////////////////////////////////////
void Window::LY_reportCC(int* bigidx, int* bigvol)
{
    int volume, surface, boxvol, maxdim;
    int nums = N;
  
    int cluster_count = 0; // by Liu Yang 07/21/04
    Vec_DP EVinfo; // added by Yang Liu 02/25/06           
  
    char s;

    *bigidx = -1;
    *bigvol = 0;

    for (int start = 0; start < nums; start++) 
    {
	// skip if already visited 
	if (spin[start] & BLACK) continue;
	
        // init statistics of this droplet. 
	s = spin[start];
	surface = 0;
	volume = LY_labelCC_getQ(start, start, s, &surface, &boxvol, &maxdim, &EVinfo);
	
	if (volume > *bigvol) 
	{
	    *bigvol = volume; 
	    *bigidx = start;
	}
	
	cout << cluster_count++ << ":  (" 
	    //   << (int)(spin[start]!=BLACK) << ")  " 	// Y.L 02/24/06, to check whether this cluster is spin-up or not 
	     << (spin[start] & 1) << ")  "    
	     << volume << "/" << surface << "/" 
	     << boxvol << "/" << maxdim << "/"
	     << EVinfo[D] << "\n";
	
    }
    //cout << "\n";
    highbitzerotwos();     // all 4th bits are now zeroed, i.e., set all GRAY or BLACK spins to WHITE spins
}
//////////////////////////////////////////////////////////////////////////////////////////////////////



//////////////////////////////////////////////////////////////////////////////////////////////////////
// The following function is based on the above function, but it considers
// the case when the largest connected cluster is spin UP. This could be true in the hysteresis model.
// Though we are considering non-spanning avalanches, but some of them could be huge, and could be the 
// largest connected cluster in the system. In this case, if we still label them with -1, then the 
// domain valume we get will not be correct. It also considers the case when the largest spin-DOWN
// connected cluster is not boundary-touching. Then we will also get the wrong domain information!
// Y.L. 02/28/06
// The idea is that: we only label the largest spin-down boundarytouching connected cluster! //02/28/06
//////////////////////////////////////////////////////////////////////////////////////////////////////
void Window::LY_modified_reportCC(int* bigidx, int* bigvol)

{
    int volume, surface, boxvol, maxdim;
    Vec_DP EVinfo; 
    bool BoundaryTouching;

    int nums = N;
  
    int cluster_count = 0; // by Liu Yang 07/21/04

    char s;

    *bigidx = -1;
    *bigvol = 0;


    ////////////////////////////////// 04/06/06  ////////////////////////////////////////////////////////////////////
    int bigidx_b=-1;  //bigest boundary-touching spin-DOWN cluster's index
    int bigvol_b=0;  //bigest boundary-touching spin-DOWN cluster's volume
    int bigidx_i=-1;  //bigest inside (non-boundary-touching) spin-DOWN cluster's index
    int bigvol_i=0;  //bigest inside (non-boundary-touching) spin-DOWN cluster's volume
    ////////////////////////////////// 04/06/06  ////////////////////////////////////////////////////////////////////


    for (int start = 0; start < nums; start++) 
    {
	// skip if already visited 
	if (spin[start] & BLACK) continue;
	
        // init statistics of this droplet. 
	s = spin[start];
	surface = 0;
	volume = LY_labelCC_getQ(start, start, s, &surface, &boxvol, &maxdim, &EVinfo, &BoundaryTouching);


	////////////////////////////////// 04/06/06  ////////////////////////////////////////////////////////////////
	// Note here, we focus on spin-DOWN (and boundarytouching) cluster, which will be labelled with -1
	// First, the label-1 cluster should be spin-DOWN, if there are spin-DOWN clusters; 
        //        but if there is no spin-DOWN cluster at all, this is the all spins spin-UP case, this will be 
        //        dealed with later.
	//
	// Second, if there are spin-DOWN clusters which are boundarytouching, then label-1 cluster should be 
        //         the largest boundarytouching spin-Down clustes. 
        //         If there is no boundarytouching spin-DOWN clusters, then label-1 cluster should be the 
        //         largest non-boundarytouching spin-DOWN clusters. 
	if ( !(spin[start] & 1) )
	{
	    if ( BoundaryTouching && (volume > bigvol_b) )
	    {
		bigvol_b = volume; 
		bigidx_b = start;
	    }

	    if ( !BoundaryTouching && (volume > bigvol_i) )
	    {
		bigvol_i = volume; 
		bigidx_i = start;
	    }
	}
	////////////////////////////////// 04/06/06  //////////////////////////////////////////////////////////////////
	    

	cout << cluster_count++ << ":  (" 
	    // << (int)(spin[start]!=BLACK) << ")  " 	// Y.L 02/24/06, to check whether this cluster is spin-up or not 
	     << (spin[start] & 1) << ")  "    
	     << volume << "/" << surface << "/" 
	     << boxvol << "/" << maxdim << "/"
	     << EVinfo[D] << "/"
	     << BoundaryTouching << "\n";
    }
    //cout << "\n";

    cout << "bigvol_b=" << bigvol_b << "  bigidx_b=" << bigidx_b << endl;
    cout << "bigvol_i=" << bigvol_i << "  bigidx_i=" << bigidx_i << endl;

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //Y.L. 04/06/06 for all the cases:
    if(bigvol_b>0)                                   // there are boundary-touching spin-DOWN clusters
    {

	*bigvol = bigvol_b; 
	*bigidx = bigidx_b;
    }
    else if(bigvol_i>0)                              // there is no boundary-touching spin-DOWN cluster
    {                                                // but there are non-boundary-touching spin-DOWN clusters

	*bigvol = bigvol_i; 
	*bigidx = bigidx_i;
    }
    else if( (cluster_count==1)  && (spin[0] & 1) )   // there is only one spin-UP cluster
    {
	//cout << "cluster_count=" << cluster_count << endl;	
	*bigvol = N; 
	*bigidx = 0;
    }
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
	

    highbitzerotwos();     // all 4th bits are now zeroed, i.e., set all GRAY or BLACK spins to WHITE spins
}
//////////////////////////////////////////////////////////////////////////////////////////////////////





//////////////////////////////////////////////////////////////////////////////////////////////////////
/* The following function is overloaded from the 
    void Window::LY_modified_report_connected_clusters(int* bigidx, int* bigvol)
   Here, we focus on the spin configurations with many avalanches(spin-UP). The following function 
   will write the information about those avalanches to a file. 03/02/06
////////////////////////////////////////////////////////////////////////////////////////////////////*/ 
//void Window::LY_modified_reportCC(int* bigidx, int* bigvol, char* filename) 
int Window::LY_modified_reportCC(int* bigidx, int* bigvol, char* filename) 
{
    ofstream outfile(filename, ios_base::out);
    outfile << "//v_CC a_CC bv_CC md_CC and  A1 A2 Asph. Prol. tr_P2 trQ_2 tr_P3 trQ_3 ev[0]...ev[D-1]\n";

    int volume, surface, boxvol, maxdim;
    Vec_DP EVinfo; 
    bool BoundaryTouching;

    int nums = N;
  
    int cluster_count = 0; // by Liu Yang 07/21/04

    int UPcluster_count = 0; // added by Yang Liu 06/16/06  

    char s;

    *bigidx = -1;
    *bigvol = 0;


    ////////////////////////////////// 04/06/06  ////////////////////////////////////////////////////////////////////
    int bigidx_b=-1;  //bigest boundary-touching spin-DOWN cluster's index
    int bigvol_b=0;  //bigest boundary-touching spin-DOWN cluster's volume
    int bigidx_i=-1;  //bigest inside (non-boundary-touching) spin-DOWN cluster's index
    int bigvol_i=0;  //bigest inside (non-boundary-touching) spin-DOWN cluster's volume
    ////////////////////////////////// 04/06/06  ////////////////////////////////////////////////////////////////////


    for (int start = 0; start < nums; start++) 
    {
	// skip if already visited 
	if (spin[start] & BLACK) continue;
	
        // init statistics of this droplet. 
	s = spin[start];
	surface = 0;
	if(spin[start]&1)
	    volume = LY_labelCC_getQ(start, start, s, &surface, &boxvol, &maxdim, &EVinfo, &BoundaryTouching);
	else
	    volume = LY_labelCC(start, start, s, &BoundaryTouching);


	////////////////////////////////// 04/06/06  ////////////////////////////////////////////////////////////////
	// Note here, we focus on spin-DOWN (and boundarytouching) cluster, which will be labelled with -1
	// First, the label-1 cluster should be spin-DOWN, if there are spin-DOWN clusters; 
        //        but if there is no spin-DOWN cluster at all, this is the all spins spin-UP case, this will be 
        //        dealed with later.
	//
	// Second, if there are spin-DOWN clusters which are boundarytouching, then label-1 cluster should be 
        //         the largest boundarytouching spin-Down clustes. 
        //         If there is no boundarytouching spin-DOWN clusters, then label-1 cluster should be the 
        //         largest non-boundarytouching spin-DOWN clusters. 
	if ( !(spin[start] & 1) )
	{
	    if ( BoundaryTouching && (volume > bigvol_b) )
	    {
		bigvol_b = volume; 
		bigidx_b = start;
	    }

	    if ( !BoundaryTouching && (volume > bigvol_i) )
	    {
		bigvol_i = volume; 
		bigidx_i = start;
	    }
	}
	////////////////////////////////// 04/06/06  //////////////////////////////////////////////////////////////////


	/*cout << cluster_count << ":  (" 
	    // << (int)(spin[start]!=BLACK) << ")  "// Y.L 02/24/06, to check whether this cluster is spin-up or not 
	     << (spin[start] & 1) << ")  "    
	     << volume << "/" << surface << "/" 
	     << boxvol << "/" << maxdim << "/"
	     << EVinfo[D] << "/"
	     << BoundaryTouching << "\n";
	*/
	if(spin[start]&1)
	{
	    outfile << volume          << " ";
	    outfile << surface         << " ";
	    outfile << boxvol          << " ";
	    outfile << maxdim          << " ";
	    
	    for(int i=0;i<D+8; i++)
		outfile << EVinfo[i] << " ";
	    
	    outfile << endl;

	    UPcluster_count++; //Y.L. 06/16/06

	}
	cluster_count++;
    }


    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //Y.L. 04/06/06 for all the cases:
    if(bigvol_b>0)                                   // there are boundary-touching spin-DOWN clusters
    {

	*bigvol = bigvol_b; 
	*bigidx = bigidx_b;
    }
    else if(bigvol_i>0)                              // there is no boundary-touching spin-DOWN cluster
    {                                                // but there are non-boundary-touching spin-DOWN clusters

	*bigvol = bigvol_i; 
	*bigidx = bigidx_i;
    }
    else if( (cluster_count==1)  && (spin[0] & 1) )   // there is only one spin-UP cluster
    {
	//cout << "cluster_count=" << cluster_count << endl;	
	*bigvol = N; 
	*bigidx = 0;
    }
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////




    highbitzerotwos(); // all 4th bits are now zeroed, i.e., set all GRAY or BLACK spins to WHITE spins
    outfile.close();

    return UPcluster_count;

}
//////////////////////////////////////////////////////////////////////////////////////////////////////







//////////////////////////////////////////////////////////////////////////////////////////////////////
/* The following function is overloaded from the 
    void Window::LY_modified_report_connected_clusters(int* bigidx, int* bigvol)
   Here, we focus on the spin configurations with just one avalanche(spin-UP). The following function 
   will pass back information about this single avalanche. 03/01/06
   here, v_CC, a_CC, bv_CC, md_CC, EV_CC carry the 
   volume, suface area, boxvol, maxdim and shape tensor information of the only one spin-UP 
   connected cluster(CC), i.e. the single avalanche.
////////////////////////////////////////////////////////////////////////////////////////////////////*/ 
void Window::LY_modified_reportCC(int* bigidx, int* bigvol, 
			          int* v_CC, int* a_CC, int* bv_CC, int* md_CC, Vec_DP* EV_CC)
{
    int volume, surface, boxvol, maxdim;
    Vec_DP EVinfo; 
    bool BoundaryTouching;

    int nums = N;
  
    int cluster_count = 0; // by Liu Yang 07/21/04

    char s;

    *bigidx = -1;
    *bigvol = 0;


    ////////////////////////////////// 04/06/06  ////////////////////////////////////////////////////////////////////
    int bigidx_b=-1;  //bigest boundary-touching spin-DOWN cluster's index
    int bigvol_b=0;  //bigest boundary-touching spin-DOWN cluster's volume
    int bigidx_i=-1;  //bigest inside (non-boundary-touching) spin-DOWN cluster's index
    int bigvol_i=0;  //bigest inside (non-boundary-touching) spin-DOWN cluster's volume
    ////////////////////////////////// 04/06/06  ////////////////////////////////////////////////////////////////////


    for (int start = 0; start < nums; start++) 
    {
	// skip if already visited 
	if (spin[start] & BLACK) continue;
	
        // init statistics of this droplet. 
	s = spin[start];
	surface = 0;

	// Y.L. 03/13/06 for spin-DOWN cluster, we don't want their structure information.
	if(spin[start]&1)
	    volume = LY_labelCC_getQ(start, start, s, &surface, &boxvol, &maxdim, &EVinfo, &BoundaryTouching);
	else
	    volume = LY_labelCC(start, start, s, &BoundaryTouching);

	

	////////////////////////////////// 04/06/06  ////////////////////////////////////////////////////////////////
	// Note here, we focus on spin-DOWN (and boundarytouching) cluster, which will be labelled with -1
	// First, the label-1 cluster should be spin-DOWN, if there are spin-DOWN clusters; 
        //        but if there is no spin-DOWN cluster at all, this is the all spins spin-UP case, this will be 
        //        dealed with later.
	//
	// Second, if there are spin-DOWN clusters which are boundarytouching, then label-1 cluster should be 
        //         the largest boundarytouching spin-Down clustes. 
        //         If there is no boundarytouching spin-DOWN clusters, then label-1 cluster should be the 
        //         largest non-boundarytouching spin-DOWN clusters. 
	if ( !(spin[start] & 1) )
	{
	    if ( BoundaryTouching && (volume > bigvol_b) )
	    {
		bigvol_b = volume; 
		bigidx_b = start;
	    }

	    if ( !BoundaryTouching && (volume > bigvol_i) )
	    {
		bigvol_i = volume; 
		bigidx_i = start;
	    }
	}
	////////////////////////////////// 04/06/06  //////////////////////////////////////////////////////////////////
	    


	// if spin-UP, this is the only one avalanche!!
	if (spin[start] & 1)
	{
	    *v_CC = volume;
	    *a_CC = surface;
	    *bv_CC= boxvol;
	    *md_CC= maxdim;
	    *EV_CC= EVinfo;
	}

	/*cout << cluster_count++ << ":  (" 
	    // << (int)(spin[start]!=BLACK) << ")  "// Y.L 02/24/06, to check whether this cluster is spin-up or not 
	     << (spin[start] & 1) << ")  "    
	     << volume << "/" << surface << "/" 
	     << boxvol << "/" << maxdim << "/"
	     << EVinfo[D] << "/"
	     << BoundaryTouching << "\n";
	*/

	cluster_count++;

    }
    //cout << "\n";

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //Y.L. 04/06/06 for all the cases:
    if(bigvol_b>0)                                   // there are boundary-touching spin-DOWN clusters
    {

	*bigvol = bigvol_b; 
	*bigidx = bigidx_b;
    }
    else if(bigvol_i>0)                              // there is no boundary-touching spin-DOWN cluster
    {                                                // but there are non-boundary-touching spin-DOWN clusters

	*bigvol = bigvol_i; 
	*bigidx = bigidx_i;
    }
    else if( (cluster_count==1)  && (spin[0] & 1) )   // there is only one spin-UP cluster
    {
	//cout << "cluster_count=" << cluster_count << endl;	
	*bigvol = N; 
	*bigidx = 0;
    }
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
	




    highbitzerotwos(); // all 4th bits are now zeroed, i.e., set all GRAY or BLACK spins to WHITE spins
}
//////////////////////////////////////////////////////////////////////////////////////////////////////



//////////////////////////////////////////////////////////////////////////////////////////////////////
/* The following function is overloaded from the 
    void Window::LY_modified_reportCC(int* bigidx, int* bigvol, 
			          int* v_CC, int* a_CC, int* bv_CC, int* md_CC, Vec_DP* EV_CC)
   Here, we still focus on the spin configurations with just one avalanche(spin-UP).
   The only difference is that this function will pass back the information of the possible holes inside the avalanche.
////////////////////////////////////////////////////////////////////////////////////////////////////*/ 
void Window::LY_modified_reportCC(int* bigidx, int* bigvol, 
				  int* v_CC, int* a_CC, int* bv_CC, int* md_CC, Vec_DP* EV_CC,
				  vector<int>* Holeinfo)
{
    int volume, surface, boxvol, maxdim;
    Vec_DP EVinfo;
    bool BoundaryTouching;
    vector<int> HoleInfo(3,0); // total volume, surface and the number of the possible holes 

    int nums = N;
  
    int cluster_count = 0; // by Liu Yang 07/21/04

    char s;

    *bigidx = -1;
    *bigvol = 0;


    ////////////////////////////////// 04/06/06  ////////////////////////////////////////////////////////////////////
    int bigidx_b=-1;  //bigest boundary-touching spin-DOWN cluster's index
    int bigvol_b=0;  //bigest boundary-touching spin-DOWN cluster's volume
    int bigidx_i=-1;  //bigest inside (non-boundary-touching) spin-DOWN cluster's index
    int bigvol_i=0;  //bigest inside (non-boundary-touching) spin-DOWN cluster's volume
    ////////////////////////////////// 04/06/06  ////////////////////////////////////////////////////////////////////


    for (int start = 0; start < nums; start++) 
    {
	// skip if already visited 
	if (spin[start] & BLACK) continue;
	
        // init statistics of this droplet. 
	s = spin[start];
	surface = 0;

	// Y.L. 03/13/06 for spin-DOWN cluster, we don't want their structure information.
	// 03/14/06, but we spin-DOWN hole, we do need their structure information to deduce the domain description of the
        // avalanche without using any functions explicitly.  
	if(spin[start]&1)
	    volume = LY_labelCC_getQ(start, start, s, 
				     &surface, &boxvol, &maxdim, 
				     &EVinfo, &BoundaryTouching);
	else                                                                     
	    volume = LY_labelCC(start, start, s, 
				&surface, &boxvol, &maxdim, 
				&BoundaryTouching);


	////////////////////////////////// 04/06/06  ////////////////////////////////////////////////////////////////
	// Note here, we focus on spin-DOWN (and boundarytouching) cluster, which will be labelled with -1
	// First, the label-1 cluster should be spin-DOWN, if there are spin-DOWN clusters; 
        //        but if there is no spin-DOWN cluster at all, this is the all spins spin-UP case, this will be 
        //        dealed with later.
	//
	// Second, if there are spin-DOWN clusters which are boundarytouching, then label-1 cluster should be 
        //         the largest boundarytouching spin-Down clustes. 
        //         If there is no boundarytouching spin-DOWN clusters, then label-1 cluster should be the 
        //         largest non-boundarytouching spin-DOWN clusters. 
	if ( !(spin[start] & 1) )
	{
	    if ( BoundaryTouching && (volume > bigvol_b) )
	    {
		bigvol_b = volume; 
		bigidx_b = start;
	    }

	    if ( !BoundaryTouching && (volume > bigvol_i) )
	    {
		bigvol_i = volume; 
		bigidx_i = start;
	    }
	}
	////////////////////////////////// 04/06/06  //////////////////////////////////////////////////////////////////
	    



	// for possible hole
	if ( !BoundaryTouching && !(spin[start] & 1) )
	{
	    HoleInfo[0] += volume;
	    HoleInfo[1] += surface;
	    HoleInfo[2]++;
	    *Holeinfo = HoleInfo;
	}

	// if spin-UP, this is the only one avalanche!!
	if (spin[start] & 1)
	{
	    *v_CC = volume;
	    *a_CC = surface;
	    *bv_CC= boxvol;
	    *md_CC= maxdim;
	    *EV_CC= EVinfo;
	}


	/*cout << cluster_count++ << ":  (" 
	    // << (int)(spin[start]!=BLACK) << ")  "// Y.L 02/24/06, to check whether this cluster is spin-up or not 
	     << (spin[start] & 1) << ")  "    
	     << volume << "/" << surface << "/" 
	     << boxvol << "/" << maxdim << "/"
	     << EVinfo[D] << "/"
	     << BoundaryTouching << "\n";
	*/

	cluster_count++;

    }




    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //Y.L. 04/06/06 for all the cases:
    if(bigvol_b>0)                                   // there are boundary-touching spin-DOWN clusters
    {

	*bigvol = bigvol_b; 
	*bigidx = bigidx_b;
    }
    else if(bigvol_i>0)                              // there is no boundary-touching spin-DOWN cluster
    {                                                // but there are non-boundary-touching spin-DOWN clusters

	*bigvol = bigvol_i; 
	*bigidx = bigidx_i;
    }
    else if( (cluster_count==1)  && (spin[0] & 1) )   // there is only one spin-UP cluster
    {
	//cout << "cluster_count=" << cluster_count << endl;	
	*bigvol = N; 
	*bigidx = 0;
    }
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
	


    highbitzerotwos(); // all 4th bits are now zeroed, i.e., set all GRAY or BLACK spins to WHITE spins
}
//////////////////////////////////////////////////////////////////////////////////////////////////////






//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/* The following function is overloaded from the 
    void Window::LY_modified_report_connected_clusters(int* bigidx, int* bigvol)
   Here, we focus on the spin configurations with many avalanches(spin-UP). The following function 
   will write the information about those avalanches to a file. 03/02/06
*/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//void Window::LY_modified_reportCC(int* bigidx, int* bigvol, bool mark, int seed, char* filename) 

int Window::LY_modified_reportCC(int* bigidx, int* bigvol, bool mark, int seed, char* filename, double H) 
{   // return UPcluster_count, Yang Liu 06/16/06

    ofstream outfile(filename, ios_base::app); // append mode
    if(mark)
	outfile << "/ seed= "<< seed<< " v_CC a_CC bv_CC md_CC and  A1 A2 Asph. Prol. tr_P2 trQ_2 tr_P3 trQ_3 ev[0]...ev[D-1] H\n";

    int volume, surface, boxvol, maxdim;
    Vec_DP EVinfo; 
    bool BoundaryTouching;

    int nums = N;
  
    int cluster_count = 0; // added by Yang Liu 07/21/04

    int UPcluster_count = 0; // added by Yang Liu 06/16/06  

    char s;

    *bigidx = -1;
    *bigvol = 0;


    ////////////////////////////////// 04/06/06  ////////////////////////////////////////////////////////////////////
    int bigidx_b=-1;  //bigest boundary-touching spin-DOWN cluster's index
    int bigvol_b=0;  //bigest boundary-touching spin-DOWN cluster's volume
    int bigidx_i=-1;  //bigest inside (non-boundary-touching) spin-DOWN cluster's index
    int bigvol_i=0;  //bigest inside (non-boundary-touching) spin-DOWN cluster's volume
    ////////////////////////////////// 04/06/06  ////////////////////////////////////////////////////////////////////



    for (int start = 0; start < nums; start++) 
    {
	// skip if already visited 
	if (spin[start] & BLACK) continue;
	
        // init statistics of this droplet. 
	s = spin[start];
	surface = 0;

	// Y.L. 03/13/06 for spin-DOWN cluster, we don't want their structure information.
	if(spin[start]&1)
	    volume = LY_labelCC_getQ(start, start, s, &surface, &boxvol, &maxdim, &EVinfo, &BoundaryTouching);
	else
	    volume = LY_labelCC(start, start, s, &BoundaryTouching);


	////////////////////////////////// 04/06/06  ////////////////////////////////////////////////////////////////
	// Note here, we focus on spin-DOWN (and boundarytouching) cluster, which will be labelled with -1
	// First, the label-1 cluster should be spin-DOWN, if there are spin-DOWN clusters; 
        //        but if there is no spin-DOWN cluster at all, this is the all spins spin-UP case, this will be 
        //        dealed with later.
	//
	// Second, if there are spin-DOWN clusters which are boundarytouching, then label-1 cluster should be 
        //         the largest boundarytouching spin-Down clustes. 
        //         If there is no boundarytouching spin-DOWN clusters, then label-1 cluster should be the 
        //         largest non-boundarytouching spin-DOWN clusters. 
	if ( !(spin[start] & 1) )
	{
	    if ( BoundaryTouching && (volume > bigvol_b) )
	    {
		bigvol_b = volume; 
		bigidx_b = start;
	    }

	    if ( !BoundaryTouching && (volume > bigvol_i) )
	    {
		bigvol_i = volume; 
		bigidx_i = start;
	    }
	}
	////////////////////////////////// 04/06/06  //////////////////////////////////////////////////////////////////




        //Y.L. 03/06/06 :   The size of the single avalanche could be obtained before we run the
        //   DropletAnalyse code. However, for the spin configuration, we are not so lucky.
	//if((spin[start]&1) && (volume >3)) 
	if(spin[start] & 1)  //Y.L. 06/16/06
	{
	    outfile << volume          << " ";
	    outfile << surface         << " ";
	    outfile << boxvol          << " ";
	    outfile << maxdim          << " ";
	    
	    for(int i=0;i<D+8; i++)
		outfile << EVinfo[i] << " ";

	    outfile << H << endl;

	    UPcluster_count++;
	}
	cluster_count++;
    }



    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //Y.L. 04/06/06 for all the cases:
    if(bigvol_b>0)                                   // there are boundary-touching spin-DOWN clusters
    {

	*bigvol = bigvol_b; 
	*bigidx = bigidx_b;
    }
    else if(bigvol_i>0)                              // there is no boundary-touching spin-DOWN cluster
    {                                                // but there are non-boundary-touching spin-DOWN clusters

	*bigvol = bigvol_i; 
	*bigidx = bigidx_i;
    }
    else if( (cluster_count==1)  && (spin[0] & 1) )   // there is only one spin-UP cluster
    {
	//cout << "cluster_count=" << cluster_count << endl;	
	*bigvol = N; 
	*bigidx = 0;
    }
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
	





    highbitzerotwos(); // all 4th bits are now zeroed, i.e., set all GRAY or BLACK spins to WHITE spins
    outfile.close();

    return UPcluster_count;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

