/////////////////////////////////////////////////////////////////////////////////////////////////////
// The following functions, adapted from the kuntz's hysteresis code,  are used for 
// the function : get_shape_tensor()  Yang Liu 02/20/06
/////////////////////////////////////////////////////////////////////////////////////////////////////
inline int Window::GetLoc(const int* coords) const
{
    int loc=0;
    loc+=coords[D-1];
    for(int i=D-2;i>=0;i--)
    {
	loc+=coords[i]*stride[i];
    }
    return loc;
}

inline void Window::GetCoords(int loc, int* newCoords) const
{
    long remainder;
    remainder = loc%stride[0];
    newCoords[0] = loc/stride[0];

    for(int i=1;i<D;i++) {
	newCoords[i] = remainder/stride[i];
	remainder %= stride[i];
    }
}

inline void Window::GetNeighbors(int loc, const int* coords) const
{
    for(int i=0;i<D;i++) {
	for(int j=0;j<D;j++) {
	    if(j==i) {
		if(coords[j]==size[j]-1)//if(coords[j]==L-1)
		{
		    neighborLocs[2*j]=loc-(size[j]-1)*stride[j];
		    neighborLocs[2*j+1]=loc-stride[j];
		    neighborCoords[2*j][j]=coords[j]-(size[j]-1);
		    neighborCoords[2*j+1][j]=coords[j]-1;
		}
		else if(coords[j]==0)
		{
		    neighborLocs[2*j]=loc+stride[j];
		    neighborLocs[2*j+1]=loc+(size[j]-1)*stride[j];
		    neighborCoords[2*j][j]=coords[j]+1;
		    neighborCoords[2*j+1][j]=coords[j]+(size[j]-1);
		}
		else
		{
		    neighborLocs[2*j]=loc+stride[j];
		    neighborLocs[2*j+1]=loc-stride[j];
		    neighborCoords[2*j][j]=coords[j]+1;
		    neighborCoords[2*j+1][j]=coords[j]-1;
		}
	    }
	    else neighborCoords[2*i+1][j]=neighborCoords[2*i][j]=coords[j];
	}
    }
}


// for given loc, this function gets its nns' locs 
inline void Window::GetNeighbors(int loc) const
{
    int* coords = new int[D];
    GetCoords(loc, coords);             // get the coordinates of the starting spin with the function 

    for(int i=0;i<D;i++) {
	for(int j=0;j<D;j++) {
	    if(j==i) {
		if(coords[j]==size[j]-1)//if(coords[j]==L-1)
		{
		    neighborLocs[2*j]=loc-(size[j]-1)*stride[j];
		    neighborLocs[2*j+1]=loc-stride[j];
		}
		else if(coords[j]==0)
		{
		    neighborLocs[2*j]=loc+stride[j];
		    neighborLocs[2*j+1]=loc+(size[j]-1)*stride[j];
		}
		else
		{
		    neighborLocs[2*j]=loc+stride[j];
		    neighborLocs[2*j+1]=loc-stride[j];
		}
	    }
	}
    }

    delete [] coords;

}



inline void Window::GetNeighbors(const Point point) const
{
    int* coords = new int[D];   

    if(D==2){
	coords[1] = point.x();	coords[0] = point.y();
    }
    else if(D==3){
	coords[2] = point.x();	coords[1] = point.y();	coords[0] = point.z();
    }

    int loc=0;
    loc+=coords[D-1];
    for(int i=D-2;i>=0;i--)
    {
	loc+=coords[i]*stride[i];
    }

    for(int i=0;i<D;i++) {
	for(int j=0;j<D;j++) {
	    if(j==i) {
		if(coords[j]==size[j]-1)//if(coords[j]==L-1)
		{
		    neighborLocs[2*j]=loc-(size[j]-1)*stride[j];
		    neighborLocs[2*j+1]=loc-stride[j];
		}
		else if(coords[j]==0)
		{
		    neighborLocs[2*j]=loc+stride[j];
		    neighborLocs[2*j+1]=loc+(size[j]-1)*stride[j];
		}
		else
		{
		    neighborLocs[2*j]=loc+stride[j];
		    neighborLocs[2*j+1]=loc-stride[j];
		}
	    }
	}
    }

    delete [] coords;
}






// The following two functions don't consider the PBC, just get the neighbors brutely.
// Y.L 02/20/06
inline void Window::GetPseudoNeighbors(const int* coords) const
{
    for(int i=0;i<D;i++) {
	for(int j=0;j<D;j++) {
	    if(j==i) 
	    {
		pseudoneighborCoords[2*j][j]=coords[j]+1;
		pseudoneighborCoords[2*j+1][j]=coords[j]-1;
	    }
	    else pseudoneighborCoords[2*i+1][j]=pseudoneighborCoords[2*i][j]=coords[j];
	}
    }
}

inline void Window::GetPseudoNeighbors(const Point point) const
{
    int* coords = new int[D];   

    if(D==2){
	coords[1] = point.x();	coords[0] = point.y();
    }
    else if(D==3){
	coords[2] = point.x();	coords[1] = point.y();	coords[0] = point.z();
    }

    for(int i=0;i<D;i++) {
	for(int j=0;j<D;j++) {
	    if(j==i) 
	    {
		pseudoneighborCoords[2*j][j]=coords[j]+1;
		pseudoneighborCoords[2*j+1][j]=coords[j]-1;
	    }
	    else pseudoneighborCoords[2*i+1][j]=pseudoneighborCoords[2*i][j]=coords[j];
	}
    }

    delete [] coords;
}

// The following funcitons can get the real location for given coordinates,
// Note that PBC are seriously considered // Yang Liu 02/20/06
inline int Window::GetLoc_PBC(const int* coords) const
{
    int* tmpCoords = new int[D];

    for(int i=D-1;i>=0;i--)
    {
	tmpCoords[i] = coords[i];

	// to make sure that tmpCoords[i] is in the range: [0,L)
	if(tmpCoords[i] >= size[i])         
	    tmpCoords[i] -= size[i];
	else if(tmpCoords[i]<0)    
	    tmpCoords[i] += size[i];
    }


    /* for debug 
       cout << "real Coords:(";
       for(int i=D-1;i>=0;i--)
       cout << tmpCoords[i] << ",";
       cout << ")";
       cout << "size:(";
       for(int i=D-1;i>=0;i--)
       cout << size[i] << ",";
       cout << ")";

       for debug */ 

    int loc=0;
    loc += tmpCoords[D-1];
    for(int i=D-2;i>=0;i--)
    {
	loc += tmpCoords[i]*stride[i];
    }

    delete [] tmpCoords;
    
    return loc;
}


inline int Window::GetLoc_PBC(const Point point) const
{
    
    int* Coords = new int[D];

    if(D==2)
    {
	Coords[1] = point.x();
	Coords[0] = point.y();
    }
    else if(D==3)
    {
	Coords[2] = point.x();
	Coords[1] = point.y();
	Coords[0] = point.z();
    }

    for(int i=D-1;i>=0;i--)
    {
	// to make sure that Coords[i] is in the range: [0,L)
	if(Coords[i] >= size[i])         
	    Coords[i] -= size[i];
	else if(Coords[i]<0)    
	    Coords[i] += size[i];
    }

    int loc=0;
    loc += Coords[D-1];
    for(int i=D-2;i>=0;i--)
    {
	loc += Coords[i]*stride[i];
    }

    delete [] Coords;
    
    return loc;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////



//////////////////////////////////////////////////////////////////////////////////////////////////////
//The follwoing function will solve the eigenvalue problem for a given point list(a cluster)   
// Note that: the function will return the eigenvalue informations:
/*
  EVinfo[0] = Anisotropy1;
  EVinfo[1] = (D==3) ? Anisotropy2 : -1; 

  EVinfo[2] = Asphericity;
  EVinfo[3] = (D==3) ? Prolateness : -1; // actually just make senses in 3D.

  EVinfo[4] = tr_P2;
  EVinfo[5] = trQ_2;
  EVinfo[6] = tr_P3;
  EVinfo[7] = trQ_3;

  for(int i=0;i<D;i++) 
  EVinfo[i+8] = lambda[i];
*/


Vec_DP Window::get_eigenvalue(list<Point>  pointList)
{
    ////////     calculate the average of x,y,z,xy,yz,zx,x^2,y^2,z^2: start   ///////
    int totM      = 0;     // total mass
    
    double x_ave  = 0;     // average of x  (center of mass: x component)
    double y_ave  = 0;     // average of y  (center of mass: y component)
    double z_ave  = 0;     // average of z  (center of mass: z component)
    
    double xy_ave = 0;     // average of xy 
    double yz_ave = 0;     // average of yz 
    double xz_ave = 0;     // average of xz

    double xx_ave = 0;     // average of x^2 
    double yy_ave = 0;     // average of y^2 
    double zz_ave = 0;     // average of z^2 

    //cout << "the pointList is:";
    double x,y,z; // for convenience of computation
    for(list<Point>::iterator p = pointList.begin(); p != pointList.end(); p++)
    {
	totM++;
        x = (*p).x();	y = (*p).y();	z = (*p).z();
	x_ave  += x;	y_ave  += y;	z_ave  += z;
	xy_ave += x*y;	yz_ave += y*z;	xz_ave += z*x;
	xx_ave += x*x;	yy_ave += y*y;	zz_ave += z*z;
	//(*p).Draw(); 
    }

    x_ave  /= totM;    y_ave  /= totM;    z_ave  /= totM;
    xy_ave /= totM;    yz_ave /= totM;    xz_ave /= totM;
    xx_ave /= totM;    yy_ave /= totM;    zz_ave /= totM;
    ////////     calculate the average of x,y,z,xy,yz,zx,x^2,y^2,z^2: end   /////////

    ////////     get the shape  tensor:start     ///////////////////////////////////////////
    double Qxx, Qxy, Qxz;
    double      Qyy, Qyz;
    double           Qzz;

    Qxx = xx_ave - x_ave*x_ave;   Qxy = xy_ave - x_ave*y_ave;   Qxz = xz_ave - x_ave*z_ave; 
    Qyy = yy_ave - y_ave*y_ave;   Qyz = yz_ave - y_ave*z_ave; 
    Qzz = zz_ave - z_ave*z_ave;
    ////////     get the shape  tensor:end     /////////////////////////////////////////////


    ////////    diagonalize the shape tensor: start ////////////////////////////////////////
    double Q_d[3*3]=
	{Qxx, Qxy, Qxz,
	 Qxy, Qyy, Qyz,
	 Qxz, Qyz, Qzz};

    int nrot;
    Mat_DP Q(Q_d,3,3);
    Vec_DP lambda(3);
    Mat_DP v(3,3);
    
    //cout << "\n\n shape tensor:\n" << Q << endl;
    NR::jacobi(Q,lambda,v,nrot); 
    NR::eigsrt(lambda,v);
    //cout << "eigenvalues = " << lambda << endl;
    ////////    diagonalize the shape tensor: end //////////////////////////////////////////

    double trQ   = 0;    
    for(int i=0;i<D;i++) 
	trQ += lambda[i];
    double trQ_2 = trQ*trQ;
    double trQ_3 = trQ_2*trQ;
    double lambda_ave = trQ/double(D);
    
    // introduce a traceless matrix P = Q - lambda_ave * I 
    double tr_P2 = 0;    
    double tr_P3 = 0;    
    //double det_P = 0;
    for(int i=0;i<D;i++) 
    {
	double temp = lambda[i]-lambda_ave; 
	tr_P2 +=  temp*temp;
	tr_P3 +=  temp*temp*temp;
	//det_P *= t;
    }

    //////////// Anisotropy defined by F. Family et. al. (PRL 1985)
    double Anisotropy1 = (lambda[D-1]/lambda[0]);
    //////////// Asphericity and Prolateness defined by Aronovitz and Nelson (JDP 1986)
    double Asphericity = D/double(D-1)*tr_P2/trQ_2;

    //////// for D>=3 ////////////////
    double Anisotropy2 = 0;
    double Prolateness = 0;
    if(D>=3) // actually just make senses in 3D.
    {
	Anisotropy2 = (lambda[D-2]/lambda[0]);
	Prolateness = D*D/double((D-1)*(D-2))*tr_P3/trQ_3;
    }

/*    cout <<"Asphericity = " << Asphericity << endl;
      cout <<"Anisotropy1 = " << Anisotropy1 << endl;
      if(D>=3)
      {
      cout <<"Anisotropy2 = " << Anisotropy2 << endl;
      cout <<"Prolateness = " << Prolateness << endl;
      }*/

    /////// store the eigenvalue information into a vector /// 
    Vec_DP EVinfo(-1, D+8); // initialize with value -1
    // Note that Anisotropy2 and Prolatesness actually just make senses in 3D.
    // so in 2D, we jsut let them be -1, as a mark.

    EVinfo[0] = Anisotropy1;
    EVinfo[1] = (D==3) ? Anisotropy2 : -1; 

    EVinfo[2] = Asphericity;
    EVinfo[3] = (D==3) ? Prolateness : -1; // actually just make senses in 3D.

    EVinfo[4] = tr_P2;
    EVinfo[5] = trQ_2;
    EVinfo[6] = tr_P3;
    EVinfo[7] = trQ_3;

    for(int i=0;i<D;i++) 
	EVinfo[i+8] = lambda[i];

    return EVinfo;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////




//////////////////////////////////////////////////////////////////////////////////////////////////////
// overloaded from the above functions,  only difference is that if volume <3, we can use pre-determined 
// eigenvalues without doing any calculation to save time  // Yang Liu 06/25/06
Vec_DP Window::get_eigenvalue(list<Point>  pointList, int size)
{
    if(size==1)
    {
	Vec_DP EVinfo(-1, D+8); // initialize with value -1
	return EVinfo;

    }// end of if size==1
    else if(size==2)
    {
	Vec_DP EVinfo(-1, D+8); // initialize with value -1
	
	if (D==2)
	{//0 -1 1 -1 0.03125 0.0625 0 0.015625 0.25 0 
	    EVinfo[0] =  0;
	    EVinfo[1] = -1; 
	    EVinfo[2] =  1;
	    EVinfo[3] = -1;
	    EVinfo[4] = 0.03125;
	    EVinfo[5] = 0.0625;
	    EVinfo[6] =  0;
	    EVinfo[7] = 0.015625;
	    EVinfo[8] = 0.25;
	    EVinfo[9] = 0;
	}
	else if(D==3)
	{//0 0 1 1 0.0416667 0.0625 0.00347222 0.015625 0.25 0 0 
	    EVinfo[0] =  0;
	    EVinfo[1] =  0; 
	    EVinfo[2] =  1;
	    EVinfo[3] =  1;
	    EVinfo[4] = 0.0416667;
	    EVinfo[5] = 0.0625;
	    EVinfo[6] = 0.00347222;
	    EVinfo[7] = 0.015625;
	    EVinfo[8] = 0.25;
	    EVinfo[9] = 0;
	    EVinfo[10]= 0;
	}
	return EVinfo;

    }// end of if size==2
    else
    {
	////////     calculate the average of x,y,z,xy,yz,zx,x^2,y^2,z^2: start   ///////
	int totM      = 0;     // total mass
    
	double x_ave  = 0;     // average of x  (center of mass: x component)
	double y_ave  = 0;     // average of y  (center of mass: y component)
	double z_ave  = 0;     // average of z  (center of mass: z component)
    
	double xy_ave = 0;     // average of xy 
	double yz_ave = 0;     // average of yz 
	double xz_ave = 0;     // average of xz

	double xx_ave = 0;     // average of x^2 
	double yy_ave = 0;     // average of y^2 
	double zz_ave = 0;     // average of z^2 

	//cout << "the pointList is:";
	double x,y,z; // for convenience of computation
	for(list<Point>::iterator p = pointList.begin(); p != pointList.end(); p++)
	{
	    totM++;
	    x = (*p).x();	y = (*p).y();	z = (*p).z();
	    x_ave  += x;	y_ave  += y;	z_ave  += z;
	    xy_ave += x*y;	yz_ave += y*z;	xz_ave += z*x;
	    xx_ave += x*x;	yy_ave += y*y;	zz_ave += z*z;
	    //(*p).Draw(); 
	}

	x_ave  /= totM;    y_ave  /= totM;    z_ave  /= totM;
	xy_ave /= totM;    yz_ave /= totM;    xz_ave /= totM;
	xx_ave /= totM;    yy_ave /= totM;    zz_ave /= totM;
	////////     calculate the average of x,y,z,xy,yz,zx,x^2,y^2,z^2: end   /////////

	////////     get the shape  tensor:start     ///////////////////////////////////////////
	double Qxx, Qxy, Qxz;
	double      Qyy, Qyz;
	double           Qzz;

	Qxx = xx_ave - x_ave*x_ave;   Qxy = xy_ave - x_ave*y_ave;   Qxz = xz_ave - x_ave*z_ave; 
	Qyy = yy_ave - y_ave*y_ave;   Qyz = yz_ave - y_ave*z_ave; 
	Qzz = zz_ave - z_ave*z_ave;
	////////     get the shape  tensor:end     /////////////////////////////////////////////


	////////    diagonalize the shape tensor: start ////////////////////////////////////////
	double Q_d[3*3]=
	    {Qxx, Qxy, Qxz,
	     Qxy, Qyy, Qyz,
	     Qxz, Qyz, Qzz};

	int nrot;
	Mat_DP Q(Q_d,3,3);
	Vec_DP lambda(3);
	Mat_DP v(3,3);
    
	//cout << "\n\n shape tensor:\n" << Q << endl;
	NR::jacobi(Q,lambda,v,nrot); 
	NR::eigsrt(lambda,v);
	//cout << "eigenvalues = " << lambda << endl;
	////////    diagonalize the shape tensor: end //////////////////////////////////////////

	double trQ   = 0;    
	for(int i=0;i<D;i++) 
	    trQ += lambda[i];
	double trQ_2 = trQ*trQ;
	double trQ_3 = trQ_2*trQ;
	double lambda_ave = trQ/double(D);
    
	// introduce a traceless matrix P = Q - lambda_ave * I 
	double tr_P2 = 0;    
	double tr_P3 = 0;    
	//double det_P = 0;
	for(int i=0;i<D;i++) 
	{
	    double temp = lambda[i]-lambda_ave; 
	    tr_P2 +=  temp*temp;
	    tr_P3 +=  temp*temp*temp;
	    //det_P *= t;
	}

	//////////// Anisotropy defined by F. Family et. al. (PRL 1985)
	double Anisotropy1 = (lambda[D-1]/lambda[0]);
	//////////// Asphericity and Prolateness defined by Aronovitz and Nelson (JDP 1986)
	double Asphericity = D/double(D-1)*tr_P2/trQ_2;

	//////// for D>=3 ////////////////
	double Anisotropy2 = 0;
	double Prolateness = 0;
	if(D>=3) // actually just make senses in 3D.
	{
	    Anisotropy2 = (lambda[D-2]/lambda[0]);
	    Prolateness = D*D/double((D-1)*(D-2))*tr_P3/trQ_3;
	}

	/////// store the eigenvalue information into a vector /// 
	Vec_DP EVinfo(-1, D+8); // initialize with value -1
	// Note that Anisotropy2 and Prolatesness actually just make senses in 3D.
	// so in 2D, we jsut let them be -1, as a mark.

	EVinfo[0] = Anisotropy1;
	EVinfo[1] = (D==3) ? Anisotropy2 : -1; 

	EVinfo[2] = Asphericity;
	EVinfo[3] = (D==3) ? Prolateness : -1; // actually just make senses in 3D.

	EVinfo[4] = tr_P2;
	EVinfo[5] = trQ_2;
	EVinfo[6] = tr_P3;
	EVinfo[7] = trQ_3;

	for(int i=0;i<D;i++) 
	    EVinfo[i+8] = lambda[i];

	return EVinfo;

    }// end of if size > 2

}
//////////////////////////////////////////////////////////////////////////////////////////////////////



//////////////////////////////////////////////////////////////////////////////////////////////////////
// Note that the following function CAN be applied to avalanches with multiple parts connected with PBC.
// Of course, it can deal with a single avalanche if it has just one part (without PBC problem). 
// However, it still can NOT deal with those spin configurations with more than one cluster. 
//  ----Yang Liu 02/20/06
//
// Basic Idea: get the pseudo cluster/avalanche first, if the cluster/avalanche has just one part, then fine.
//             otherwise, to calculate the anisotropy, we have to shift somes part to image cells to make 
//             a pseudo cluster/avalnche with only one part.
Vec_DP Window::get_shape_tensor() 
{
    queue<Point> pointQueue;            // store the coordinates of the pseudo cluster/avalanche
    list<Point>  pointList;             // store the coordinates of the pseudo cluster/avalanche

    char* mark = new char[N];           // store the marks telling whether the spin has been scanned 
    for(int i=0; i<N; i++) mark[i] = 0; // initialize the mark array

    int Loc = startingSpinLoc;          // the starting spin when scanning the cluster/avalanche  
    mark[Loc] = 1;                      // mark the starting spin
    int* Coords = new int[D];           // the coordinates of the spin  
    GetCoords(Loc, Coords);             // get the coordinates of the starting spin 
    Point point(D, Coords);                
    pointList.push_back(point);         // push it into the queue and the list  
    pointQueue.push(point);

    //cout << "starting point=";
    //point.Draw();
    //cout << endl;

    while(!pointQueue.empty())         
    {
	Point point = pointQueue.front();  // get the point from the queue 
 	pointQueue.pop();                  // pop it from the queue

	GetPseudoNeighbors(point);         // get its pseudo neighbors  

	// loop over this spin's pseudo neighbors and record them 
	// if the neighbor is also spin-up and has not been pushed into the queue 		
	for(int i=0;i<Z;i++) 
	{
	    int neighborLoc = GetLoc_PBC(pseudoneighborCoords[i]); // get pseudo neighbors' real locations. 
	    if( spin[neighborLoc]==SPINUP && mark[neighborLoc]==0) // if the corresponding spin is up,  
	    {                                                      // which means it belongs to the avalanche, and it has 
                                                                   // not been pushed into the queue, 
		mark[neighborLoc] = 1;                             // then mark it 
		Point point(D, pseudoneighborCoords[i]);           
		pointList.push_back(point);                        // and push it into the queue and the list. 
		pointQueue.push(point);
	    }
	}// end of for each neighbor loop
    }

    delete [] Coords;
    delete [] mark;

    // Once we get the point list of the cluster/avalanche, the following work is trivial. !!!
    return (get_eigenvalue(pointList));
}
//////////////////////////////////////////////////////////////////////////////////////////////////////
