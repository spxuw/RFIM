void perpsquare(int colorcode, int i, int j, int k, int in, int jn, int kn,  ofstream & gout)
{
    // If two points are not too different, draw a square separating them.
    // Builds up a surface.

    // Color code according to whether 1, -1 or 2
    double r = 0;
    double g = 0;
    double b = 0;
  
    /*if (colorcode == 1) {r = g = b = 1.0;}
      else {if (colorcode < 0)  { r = g = 1.0; b = 0.0; }
      else { b = 1.0; r = g = 0.0; }}*/
    switch(colorcode) // by yang liu 07/19/04
    {
	case 1:  r = g = b = 1.0;     break;
	case -1: r = g = 1.0; b = 0.0; break;
	case 2:  b = 1.0; r = g = 0.0; break;
    }

    double csize = 0.5;
  
    if (abs(i-in) > 1 || abs (j - jn) > 1 || abs (k - kn) > 1) return;
    if (k != kn) {// r = 0.0;
	gout << i - csize << " " << j - csize << " " << k + (kn - k) * csize << " " << r << " " << g << " " << b << " 1\n";
	gout << i + csize << " " << j - csize << " " << k + (kn - k) * csize << " " << r << " " << g << " " << b << " 1\n";
	gout << i + csize << " " << j + csize << " " << k + (kn - k) * csize << " " << r << " " << g << " " << b << " 1\n";
	gout << i - csize << " " << j + csize << " " << k + (kn - k) * csize << " " << r << " " << g << " " << b << " 1\n";
    }
    if (j != jn) {// g = 0.0;
	gout << i - csize << " " << j + (jn - j) * csize << " " << k - csize << " " << r << " " << g << " " << b << " 1\n";
	gout << i + csize << " " << j + (jn - j) * csize << " " << k - csize << " " << r << " " << g << " " << b << " 1\n";
	gout << i + csize << " " << j + (jn - j) * csize << " " << k + csize << " " << r << " " << g << " " << b << " 1\n";
	gout << i - csize << " " << j + (jn - j) * csize << " " << k + csize << " " << r << " " << g << " " << b << " 1\n";
    }
    if (i != in) {// b = 0.0;
	gout << i + (in - i) * csize<< " " << j - csize << " " << k - csize << " " << r << " " << g << " " << b << " 1\n";
	gout << i + (in - i) * csize<< " " << j + csize << " " << k - csize << " " << r << " " << g << " " << b << " 1\n";
	gout << i + (in - i) * csize<< " " << j + csize << " " << k + csize << " " << r << " " << g << " " << b << " 1\n";
	gout << i + (in - i) * csize<< " " << j - csize << " " << k + csize << " " << r << " " << g << " " << b << " 1\n";
    }
}

void perpsquare(int i, int j, int k, int in, int jn, int kn,  ofstream & gout)
{
    // If two points are not too different, draw a square separating them.
    // Builds up a surface.
    double r, g, b;
    r = g = b = 1.0;
    double csize = 0.5;
  
    if (abs(i-in) > 1 || abs (j - jn) > 1 || abs (k - kn) > 1) return;
    if (k != kn) { r = 0.0;
    gout << i - csize << " " << j - csize << " " << k + (kn - k) * csize << " " << r << " " << g << " " << b << " 1\n";
    gout << i + csize << " " << j - csize << " " << k + (kn - k) * csize << " " << r << " " << g << " " << b << " 1\n";
    gout << i + csize << " " << j + csize << " " << k + (kn - k) * csize << " " << r << " " << g << " " << b << " 1\n";
    gout << i - csize << " " << j + csize << " " << k + (kn - k) * csize << " " << r << " " << g << " " << b << " 1\n";
    }
    if (j != jn) { g = 0.0;
    gout << i - csize << " " << j + (jn - j) * csize << " " << k - csize << " " << r << " " << g << " " << b << " 1\n";
    gout << i + csize << " " << j + (jn - j) * csize << " " << k - csize << " " << r << " " << g << " " << b << " 1\n";
    gout << i + csize << " " << j + (jn - j) * csize << " " << k + csize << " " << r << " " << g << " " << b << " 1\n";
    gout << i - csize << " " << j + (jn - j) * csize << " " << k + csize << " " << r << " " << g << " " << b << " 1\n";
    }
    if (i != in) { b = 0.0;
    gout << i + (in - i) * csize<< " " << j - csize << " " << k - csize << " " << r << " " << g << " " << b << " 1\n";
    gout << i + (in - i) * csize<< " " << j + csize << " " << k - csize << " " << r << " " << g << " " << b << " 1\n";
    gout << i + (in - i) * csize<< " " << j + csize << " " << k + csize << " " << r << " " << g << " " << b << " 1\n";
    gout << i + (in - i) * csize<< " " << j - csize << " " << k + csize << " " << r << " " << g << " " << b << " 1\n";
    }
}

void perpsquaredelta(int i, int j, int k,
		     int di, int dj, int dk,  ofstream & gout)
{
    // If two points are not too different, draw a square separating them.
    // Builds up a surface.
    double r, g, b;
    r = g = b = 1.0;
    double csize = 0.5;
  
    //  cout << "psd " << i << "," << j << "," << k << ":" << di << dj << dk << " ";
    //  cout << "\n";
    if (dk != 0) { r = 0.0;
    gout << i - csize << " " << j - csize << " " << k + dk * csize << " " << r << " " << g << " " << b << " 1\n";
    gout << i + csize << " " << j - csize << " " << k + dk * csize << " " << r << " " << g << " " << b << " 1\n";
    gout << i + csize << " " << j + csize << " " << k + dk * csize << " " << r << " " << g << " " << b << " 1\n";
    gout << i - csize << " " << j + csize << " " << k + dk * csize << " " << r << " " << g << " " << b << " 1\n";
    }
    if (dj != 0) { g = 0.0;
    gout << i - csize << " " << j + dj * csize << " " << k - csize << " " << r << " " << g << " " << b << " 1\n";
    gout << i + csize << " " << j + dj * csize << " " << k - csize << " " << r << " " << g << " " << b << " 1\n";
    gout << i + csize << " " << j + dj * csize << " " << k + csize << " " << r << " " << g << " " << b << " 1\n";
    gout << i - csize << " " << j + dj * csize << " " << k + csize << " " << r << " " << g << " " << b << " 1\n";
    }
    if (di != 0) { b = 0.0;
    gout << i + di * csize<< " " << j - csize << " " << k - csize << " " << r << " " << g << " " << b << " 1\n";
    gout << i + di * csize<< " " << j + csize << " " << k - csize << " " << r << " " << g << " " << b << " 1\n";
    gout << i + di * csize<< " " << j + csize << " " << k + csize << " " << r << " " << g << " " << b << " 1\n";
    gout << i + di * csize<< " " << j - csize << " " << k + csize << " " << r << " " << g << " " << b << " 1\n";
    }
}

// draw two squares displaced x, y, or z plane, so can draw up to 3 distinct sets.
void cube(int i, int j, int k, int spin, int plane,  ofstream & gout)
{
    double r, g, b;
    double csize = 0.15;
    if (spin == 0) {r = 1.0; g = b = 0.0; return;}
    else  {g = 1.0; r = b = 0.0;}
    switch (plane) {
	case 0:
	    gout << i - csize << " " << j - csize << " " << k << " " << r << " " << g << " " << b << " 1\n";
	    gout << i + csize << " " << j - csize << " " << k << " " << r << " " << g << " " << b << " 1\n";
	    gout << i + csize << " " << j + csize << " " << k << " " << r << " " << g << " " << b << " 1\n";
	    gout << i - csize << " " << j + csize << " " << k << " " << r << " " << g << " " << b << " 1\n";
	    break;
	case 1:
	    gout << i - csize << " " << j << " " << k - csize << " " << r << " " << g << " " << b << " 1\n";
	    gout << i + csize << " " << j << " " << k - csize << " " << r << " " << g << " " << b << " 1\n";
	    gout << i + csize << " " << j << " " << k + csize << " " << r << " " << g << " " << b << " 1\n";
	    gout << i - csize << " " << j << " " << k + csize << " " << r << " " << g << " " << b << " 1\n";
	    break;
	case 2:
	    gout << i << " " << j - csize << " " << k - csize << " " << r << " " << g << " " << b << " 1\n";
	    gout << i << " " << j + csize << " " << k - csize << " " << r << " " << g << " " << b << " 1\n";
	    gout << i << " " << j + csize << " " << k + csize << " " << r << " " << g << " " << b << " 1\n";
	    gout << i << " " << j - csize << " " << k + csize << " " << r << " " << g << " " << b << " 1\n";
	    break;
    }
}

void Window::geomview(char *filename, int plane)
{
    int spincount = 0; 
    ofstream geom_out(filename);
  
    geom_out << "CQUAD\n";
    for (int k = 0; k < H; ++k) {
	for (int j = 0; j < L; ++j) {
	    for (int i = 0; i < W; ++i) {
		if (spin[spincount] == SPINUP)
		    cube(i, j, k, 1, plane, geom_out);
		else
		    cube(i, j, k, 0, plane, geom_out);
		spincount++;
	    }
	}
    }
}


//////////////////////////////////////////////////////////////////////////////////////////////////////
void Window::postlayer(char *fname)
{
    int spincount = 0;
    ofstream fout(fname);

    fout << "%!Adobe-PS1.0\n";
    fout << "newpath 50 100 moveto 500 0 rlineto 0 500 rlineto -500 0 rlineto 0 -500 rlineto stroke\n";
    fout << "/s {newpath moveto 1 0 rlineto 0 1 rlineto -1 0 rlineto closepath fill} def\n";
    fout << "50 100 translate\n" << 500.0/L << " " << 500.0/L << " scale\n";
    /* k = 0 layer. */
    for (int j = L-1; j >=0; j--) {//Y.L. 03/02/06
//    for (int j = 0; j < L; ++j) {
	for (int i = 0; i < W; ++i) {
	    if (spin[spincount] == SPINUP)
		fout << i << " " << j << " s\n";
	    spincount++;
	}
    }
    fout.close();
}


//////////////////////////////////////////////////////////////////////////////////////////////////////
// assumes periodic in y, z directions
int Window::depthtrace(int i, int j, int k) {
    int idx;		// index of spin
    int in, jn, kn;

    // set this spin to 2
    // follow through spins that have spin 0 (xor gives 1 when spins differ.)
    // count "dead ends": start from two, facing onto 1.

    idx = i + W * (j + L * k);
    //  cout << i << " " << j << " " << k << " " << int(spin[idx]) << "\n";
    if (spin[idx] == 1) return 1; // facing onto terminating site.
    if (spin[idx] != 0) return 0; // this site has already been visited

    int count = 0;
    spin[idx] = 2;

    in = i;
    jn = j;
    kn = k;
    if (i < W - 1) {
	in = i + 1;
	count += depthtrace(in, jn, kn);
    }
    if (i > 0) {
	in = i - 1;
	count += depthtrace(in, jn, kn);
    }
    in = i;

    jn = (j + L - 1) % L;
    count += depthtrace(in, jn, kn);
    jn = (j + 1) % L;
    count += depthtrace(in, jn, kn);
    jn = j;

    kn = (k + H - 1) % H;
    count += depthtrace(in, jn, kn);
    kn = (k + 1) % H;
    count += depthtrace(in, jn, kn);
    kn = k;

    return count;
}
  
// assumes periodic in y, z directions
// maybe too much recursion on large systems (128^3)?
// stack overflow, most likely
// override with ulimit?


//////////////////////////////////////////////////////////////////////////////////////////////////////
int Window::highbitdepthtrace(int i, int j, int k) {
    int idx;		// index of spin
    int in, jn, kn;

    // set this spin to 32
    // follow through spins that have spin 0 (xor gives 1 when spins differ.)
    // count "dead ends": start from two, facing onto 1.

    idx = i + W * (j + L * k);
    //  cout << i << " " << j << " " << k << " " << int(spin[idx]) << "\n";
    if (spin[idx] & 16) return 1; // facing onto terminating site.
    if (spin[idx] & 32) return 0; // this site has already been visited

    int count = 0;
    spin[idx] |= 32;		// mark as visited

    in = i;
    jn = j;
    kn = k;
    if (i < W - 1) {
	in = i + 1;
	count += highbitdepthtrace(in, jn, kn);
    }
    if (i > 0) {
	in = i - 1;
	count += highbitdepthtrace(in, jn, kn);
    }
    in = i;

    jn = (j + L - 1) % L;
    count += highbitdepthtrace(in, jn, kn);
    jn = (j + 1) % L;
    count += highbitdepthtrace(in, jn, kn);
    jn = j;

    kn = (k + H - 1) % H;
    count += highbitdepthtrace(in, jn, kn);
    kn = (k + 1) % H;
    count += highbitdepthtrace(in, jn, kn);
    kn = k;

    return count;
}



int Window::highbitbreadthtrace(int i, int j, int k) {
    int idx, idxn;		// index of spin, nbr
    int in, jn, kn;
    int headQ;
    list <int> Q;
    // assume all vertices are "white" - have 0 in bits 4 and 5

    int count = 0;

    idx = i + W * (j + L * k);
    spin[idx] |= GRAY;		// "gray"

    Q.push_back(idx);
    while (!Q.empty()) 
    {
	headQ = *(Q.begin());
	Q.pop_front();
	
        // get coordinates
	k = kn = headQ / L / W;
	j = jn = (headQ / W) % L;
	i = in = headQ % W;

	// now look at n.n.: if white, color gray and queue
	if (i < W - 1) 
	{
	    in = i + 1;
	    idxn = in + W * (jn + L * kn);
	    if ((spin[idxn] & NOGO)) { //don't include! but increment count
		count++;
	    } 
	    else if ((spin[idxn] & NOTWHITE) == 0)
	    {spin[idxn] |= GRAY; Q.push_back(idxn);}
	}
	if (i > 0) {
	    in = i - 1;
	    idxn = in + W * (jn + L * kn);
	    if ((spin[idxn] & NOGO)) { //don't include! but increment count
		count++;
	    } else
		if ((spin[idxn] & NOTWHITE) == 0){spin[idxn]|=GRAY;Q.push_back(idxn);}
	}
	in = i;

	jn = (j + L - 1) % L;
	idxn = in + W * (jn + L * kn);
	if ((spin[idxn] & NOGO)) { //don't include! but increment count
	    count++;
	} else
	    if ((spin[idxn] & NOTWHITE) == 0){spin[idxn]|=GRAY;Q.push_back(idxn);}
	jn = (j + 1) % L;
	idxn = in + W * (jn + L * kn);
	if ((spin[idxn] & NOGO)) { //don't include! but increment count
	    count++;
	} else
	    if ((spin[idxn] & NOTWHITE) == 0){spin[idxn]|=GRAY;Q.push_back(idxn);}
	jn = j;
    
	kn = (k + H - 1) % H;
	idxn = in + W * (jn + L * kn);
	if ((spin[idxn] & NOGO)) { //don't include! but increment count
	    count++;
	} else
	    if ((spin[idxn] & NOTWHITE) == 0){spin[idxn]|=GRAY;Q.push_back(idxn);}
	kn = (k + 1) % H;
	idxn = in + W * (jn + L * kn);
	if ((spin[idxn] & NOGO)) { //don't include! but increment count
	    count++;
	} else
	    if ((spin[idxn] & NOTWHITE) == 0){spin[idxn]|=GRAY;Q.push_back(idxn);}
	kn = k;
    
	spin[headQ] -= GRAY; // better have this bit set! (should by algorithm)
	spin[headQ] |= BLACK;
    }
    // all 4th bits are now zeroed
    highbitzerotwos();
    return count;
}
  

//////////////////////////////////////////////////////////////////////////////////////////////////////
void Window::highbitzerotwos() {// set all thirtytwos back to zero
    // used for cleaning up after highbitdepthtrace
    int idx;
    int i, j, k;

    for (k = 0; k < H; ++k)
	for (j = 0; j < L; ++j)
	    for (i = 0; i < W; ++i) {
		idx = i + W * (j + L * k);
		if (spin[idx] & GRAY) spin[idx] -= GRAY;
		if (spin[idx] & BLACK) spin[idx] -= BLACK;
	    }
}


//////////////////////////////////////////////////////////////////////////////////////////////////////
void Window::zerotwos() {// set all twos back to zero
    // used for cleaning up after depthtrace
    int idx;
    int i, j, k;

    for (k = 0; k < H; ++k)
	for (j = 0; j < L; ++j)
	    for (i = 0; i < W; ++i) {
		idx = i + W * (j + L * k);
		if (spin[idx] == 2) spin[idx] = 0;
	    }
}


//////////////////////////////////////////////////////////////////////////////////////////////////////
void Window::geomviewsurf(char *filename)
    // geomview of the surface from 1 to 0 spin, 
    // used for viewing comparisons.
{
    int spincount = 0; 
    ofstream geom_out(filename);
  
    geom_out << "CQUAD\n";
    // red (x=0) green (x=W) bounding planes at x = 0, W, for geomview reference
    //  geom_out << "0 0 0 1 0 0 1\n";
    //  geom_out << "0 0 " << H << " 1 0 0 1\n";
    //  geom_out << "0 " << L << " " << H << " 1 0 0 1\n";
    //  geom_out << "0 " << L << " 0 1 0 0 1\n";
    //  geom_out << W << " 0 0 0 1 0 1\n";
    //  geom_out << W << " 0 " << H << " 0 1 0 1\n";
    //  geom_out << W << " " << L << " " << H << " 0 1 0 1\n";
    //  geom_out << W << " " << L << " 0 0 1 0 1\n";
  
    int i, j, k, di, dj, dk, dir;
    for (k = 0; k < H; ++k)
	for (j = 0; j < L; ++j)
	    for (i = 0; i < W; ++i)
		for (dir = 0; dir < 6; ++dir) {
		    spincount = i + W * (j + L * k);
		    if (spin[spincount]) {
			di = dj = dk = 0;
			switch (dir) {
			    case 0: di = 1; break;
			    case 1: dj = 1; break;
			    case 2: dk = 1; break;
			    case 3: dk = -1; break;
			    case 4: dj = -1; break;
			    case 5: di = -1; break;
			}
			int nbrcount = ((i+di+W)%W)
			    + W * ( ((j+dj+L)%L) + L * ((k+dk+H)%H) );
			if (!spin[nbrcount])
			    perpsquaredelta(i, j, k, di, dj, dk, geom_out);
		    }
		}
}


//////////////////////////////////////////////////////////////////////////////////////////////////////
void errsym(char *x) { fprintf(stderr,"Error in symmetric diff: %s\n", x); exit(1);}


//////////////////////////////////////////////////////////////////////////////////////////////////////
// Symmetric difference on Windows.
// Note that + operator is now ^ operator, before 4/14/99, was !(A^B)
Window operator + (const Window & A, const Window & B)
{
    Window x(A);

    if (B.W != x.W) errsym("Bad coords");
    if (B.L != x.L) errsym("Bad coords");
    if (B.H != x.H) errsym("Bad coords");
    /* don't do anything with Wsrc, Lsrc, Hsrc */
    int nspins = A.W * A.L * A.H;

    for (int i = 0; i < nspins; ++i)
	x.spin[i] = (A.spin[i]) ^ (B.spin[i]);

    return x;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////
// Symmetric difference on Windows.
// Note that + operator is now ^ operator, before 4/14/99, was !(A^B)
void Window::myxor(const Window & A, const Window & B)
{
    if (B.W != W) errsym("Bad coords");
    if (B.L != L) errsym("Bad coords");
    if (B.H != H) errsym("Bad coords");
    if (A.W != W) errsym("Bad coords");
    if (A.L != L) errsym("Bad coords");
    if (A.H != H) errsym("Bad coords");
    /* don't do anything with Wsrc, Lsrc, Hsrc */
    int nspins = A.W * A.L * A.H;

    for (int i = 0; i < nspins; ++i)
	spin[i] = (A.spin[i]) ^ (B.spin[i]);
}

//////////////////////////////////////////////////////////////////////////////////////////////////////
// Symmetric difference on Windows.
// Note that + operator is now ^ operator, before 4/14/99, was !(A^B)
// xor this and A, putting result into bit 4 of this
void Window::highbitxor(const Window & A)
{
    if (A.W != W) errsym("Bad coords");
    if (A.L != L) errsym("Bad coords");
    if (A.H != H) errsym("Bad coords");
    /* don't do anything with Wsrc, Lsrc, Hsrc */
    int nspins = A.W * A.L * A.H;

    for (int i = 0; i < nspins; ++i) {
	// set bit 4, if 0th bit xor
	spin[i] |= 16;		// set by default
	if ( !((spin[i] ^ A.spin[i]) & 1) ) { // unset the bit if not distinct
	    spin[i] -= 16;
	}
    }
}




//////////////////////////////////////////////////////////////////////////////////////////////////////
int Window::getnbr(int start, int dir)
{
    int startz, starty, startx, nbrx, nbry, nbrz;
    startz = (start)/(W*L);
    starty = (start - startz*W*L)/W;
    startx = (start) % W;
    nbrx = startx; nbry = starty; nbrz = startz;
    switch (dir) {
	case 0: nbrx = (startx + 1) % W; break;
	case 1: nbry = (starty + 1) % L; break;
	case 2: nbrz = (startz + 1) % H; break;
	case 3: nbrz = (startz + H - 1) % H; break;
	case 4: nbry = (starty + L - 1) % L; break;
	case 5: nbrx = (startx + W - 1) % W; break;
    }
    return nbrx + W * (nbry + L * nbrz);
}


//////////////////////////////////////////////////////////////////////////////////////////////////////
// Note that this code cannot be applied to avalanches with multiple parts connected 
// with periodical boundaries. ----Yang Liu 06/28/05
// It can deal with a single avalanche if there is no pbc problem. 
// It can also deal with several clusters if there is no pbc problem.
void Window::get_inertia_tensor()
{

    list<Point> pointList;

    ////////     calculate the center of mass   /////// 
    int totM = 0;     // total mass
    int xc= 0;     // center of mass
    int yc= 0;    
    int zc= 0;

    cout << "\nthe pointList is:";    
    for (int k = 0; k < H; k++) {
	for (int j = 0; j < L; j++) {
	    for (int i = 0; i < W; i++) {
		if (spin[i + W * (j + L * k)] == SPINUP) {
		    Point point(i,j,k);
		    pointList.push_back(point);
		    totM++;
		    xc += i;
		    yc += j;
		    zc += k;

		    point.Draw(); 
		}
	    }
	}
    }

    xc /= totM;
    yc /= totM;
    zc /= totM;
    ////////     calculate the center of mass   /////// 


    ////////     get the inertia tensor: start         ///////
    double Ixx, Ixy, Ixz;
    double      Iyy, Iyz;
    double           Izz;
    Ixx=Ixy=Ixz=Iyy=Iyz=Izz=0;

    for(list<Point>::iterator p = pointList.begin(); p != pointList.end(); p++)
    {
	(*p).Shift(xc,yc,zc); //    shift the coordinates into the fram of center-of-mass   

	Ixx += (*p).GetIxx();
	Ixy += (*p).GetIxy();
	Ixz += (*p).GetIxz();
	Iyy += (*p).GetIyy();
	Iyz += (*p).GetIyz();
	Izz += (*p).GetIzz();
    }
    ////////    get the inertia tensor:  end       ///////


    ////////    diagonalize //////////////////////////
    double I_d[3*3]=
	{Ixx, Ixy, Ixz,
	 Ixy, Iyy, Iyz,
	 Ixz, Iyz, Izz};

    int nrot;
    Mat_DP I(I_d,3,3);
    Vec_DP d(3);
    Mat_DP v(3,3);

    
    cout << "\n\n inertia tensor:\n";
    cout << I;
    NR::jacobi(I,d,v,nrot);
    NR::eigsrt(d,v);

    cout << "eigenvalues: " << endl;
    for (int j=0;j<3;j++) 
    {
	cout << ' ' << d[j];
	if ((j+1) % 5 == 0) cout << endl;
    }

    cout << endl << "eigenvectors:" << endl;
    for (int j=0;j<3;j++) 
    {
	cout << "number" << (j+1);
	for (int k=0;k<3;k++) 
	{
	    cout << ' ' << v[k][j];
	    if ((k+1) % 5 == 0) cout << endl;
	}
	cout << endl;
    }


}


