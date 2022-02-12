#ifndef DROPLETANALYSE_H
#define DROPLETANALYSE_H
#define SPINUP 1
#define SPINDN 0
#include <stdio.h>
#include <fstream>
#include <iostream>
#include <queue>
#include <vector>
#include <list>
#include "nr.h"                // added by Yang Liu 06/27/05
#include "Point.h"             // added by Yang Liu 07/20/05
#include "State.h"             // added by Yang Liu 08/14/06 
using namespace std;

class Window
{
 public:
    Window(std::istream& x);                           // constructor: init window from an ASCII stream
    Window(const Window& x);                        // constructor: init window from another window 
    Window(int dimension, int length, list<int>& spinList); // constructor: init window from a list
    Window(int dimension, int length, char* spinCluster); // constructor: init window from an array
    Window(const State& x);                         // constructor: init window from a state x

    ~Window() 
	{
	    delete[] spin; 
	    if (label != NULL) 
		delete[] label;
	    
	    delete [] stride;
	    delete [] size;//added by Y.L. 02/26/06
	    delete [] neighborLocs;
	    for(int i=0;i<Z;i++)
	    {
		delete [] neighborCoords[i];
		delete [] pseudoneighborCoords[i];
	    }
	    delete [] neighborCoords;
	    delete [] pseudoneighborCoords;
	};
  

    /////////////////////     some small tool-functions    ///////////////////////////////////////////////////////////////
    void printme();
    void geomview(char* fname, int plane);             // for viewing by geomview
    void geomviewsurf(char* fname);
    void postlayer(char* fname);                       // output a layer for postscript file.
    int* slicecount(int dir);
    int* slicecountcenter(int dir);
    /* count up how many bonds on wall, moving through regions of zeros starting from the indexed boundary wall */
    int wallarea(int startboundary);  
    /* count up how many bonds on wall, moving through regions of zeros in bit 4,starting from the indexed boundary wall*/
    int highbitwallarea(int startboundary); 
    int depthtrace(int i, int j, int k);
    int highbitdepthtrace(int i, int j, int k);
    int highbitbreadthtrace(int i, int j, int k);
    friend Window operator+(const Window& A, const Window& B);
    void myxor(const Window& A, const Window& B);
    void highbitxor(const Window& A);
    int getnbr(int start, int dir);
    void zerotwos();
    void highbitzerotwos();
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    //////////////////////////////  AAM's original droplet-analysis code  ///////////////////////////////////////////////
    int label_connected_cluster(int clustindex, int newlabel, char s, int* surface, int* boxvol, int* maxdim);
    void report_connected_clusters(int* bigidx, int* bigvol);
    int domain_with_notlabel(int clustindex, int notlabel, int* surface, int* boxvol, int* maxdim);
    void report_domains_with_notlabel(int notlabel);
    void print_bubble_surfaces();
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


    //////// the following tool-funcitons are for the calculation of the shape tensor        /////////////// ////////////
    // adopted from the Hysteresis Code. Y.L.
    void GetNeighbors(int currentLoc) const; //Y.L. 02/27/05
    void GetNeighbors(const Point point) const;// Y.L. 02/27/05
    void GetNeighbors(int currentLoc, const int* currentCoords) const;
    void GetPseudoNeighbors(const int* currentCoords) const;// added by Yang Liu 07/20/05
    void GetPseudoNeighbors(const Point point) const;// added by Yang Liu 07/22/05
    void GetCoords(int loc, int* newCoords) const;
    int  GetLoc(const int* coords) const;
    int  GetLoc_PBC(const int* coords) const; // added by Yang Liu 07/20/05
    int  GetLoc_PBC(const Point point) const; // added by Yang Liu 02/26/06
    void get_inertia_tensor();
    Vec_DP get_shape_tensor();
    Vec_DP get_eigenvalue(list<Point>  pointList);
    Vec_DP get_eigenvalue(list<Point>  pointList, int size); // added by Yang Liu 06/25/06
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



    ///////////  LY's droplet-analysis code: more accurate in calculation of maxdim and domain volume        //////////
    //1, for connected-cluster description
    //  see DropletAnalyse_LY_labelCC.cpp
    int  LY_labelCC(int clustindex, int newlabel, char s, int* surface, int* boxvol, int* maxdim);
    int  LY_labelCC(int clustindex, int newlabel, char s, int* surface, int* boxvol, int* maxdim, bool* Touch);
    int  LY_labelCC_getQ(int loc, int newlabel, char s, int* surface, int* boxvol, int* maxdim, Vec_DP* EV);
    int  LY_labelCC_getQ(int loc, int newlabel, char s, int* surface, int* boxvol, int* maxdim, Vec_DP* EV, bool* Touch);

    int  LY_labelCC(int loc, int newlabel, char s, bool* BoundaryTouching);  //for spin-DOWN connected cluster
    int  LY_labelCC_getQ(int loc, int newlabel, char s, 
			 int* surface, int* boxvol, int* maxdim, 
			 Vec_DP* EV, bool* Touch,
			 list<Point>* newList, vector<int>* newSize);


    //  see DropletAnalyse_LY_reportCC.cpp
    void LY_reportCC(int* bigidx, int* bigvol);
    void LY_modified_reportCC(int* bigidx,int* bigvol);

    //void LY_modified_reportCC(int* bigidx, int* bigvol, char* filename) ;
    int LY_modified_reportCC(int* bigidx, int* bigvol, char* filename) ; //Yang Liu 06/16/06

    void LY_modified_reportCC(int* bigidx,int* bigvol,
			      int* v,int* a,int* bv,int* md,Vec_DP* EV); // for single aval.
    void LY_modified_reportCC(int* bigidx,int* bigvol,
			      int* v,int* a,int* bv,int* md, Vec_DP* EV,
			      list<Point>* newList, vector<int>* newSize); // for single aval.
    void LY_modified_reportCC(int* bigidx,int* bigvol,
			      int* v,int* a,int* bv,int* md, Vec_DP* EV,
			      vector<int>* Holeinfo); // for single aval.

    //void LY_modified_reportCC(int* bigidx, int* bigvol, bool mark, int seed, char* filename);    // for spin config.
    int LY_modified_reportCC(int* bigidx, int* bigvol, bool mark, int seed, char* filename, double H);    
    //  return the number of UP clusters, Yang Liu 06/16/06


    //2, for domain description
    //  see DropletAnalyse_LY_labelD.cpp
    int  LY_domain_with_notlabel(int loc, int notlabel, int* surface, int* boxvol, int* maxdim);
    int  LY_domain_with_notlabel(int loc, int notlabel, int* surface);
    int  LY_domain_with_notlabel(int loc, int notlabel); //for spin-DOWN domain

    //  see DropletAnalyse_LY_reportD.cpp
    void LY_reportD_with_notlabel(int notlabel);  
    void LY_reportD_with_notlabel(int depth, int notlabel,int* v_D, int* a_D, int* bv_D, int* md_D);
    void LY_reportD_with_notlabel(int depth, int notlabel,int* v_D, int* a_D);
    void LY_reportD_with_notlabel(int depth, int notlabel, char* filename);  
    void LY_reportD_with_notlabel(int depth, int notlabel, bool mark, int seed, char* filename, double H); 


    //3, the higher level functions, see DropletAnalyse_LY.cpp
    void LY_print_bubble_surfaces();
    void LY_analyse_one_avalanche();
    void LY_analyse_one_avalanche(char* filename);   
    int  LY_analyse_many_avalanches(char* filename); 
    void LY_analyse_one_avalanche(bool mark, int seed, char* filename);   // use in the Hysteresis code

    int  LY_analyse_many_avalanches(bool mark, int seed, char* mainname, double H); // use in the Hysteresis code
    // return the number of UP clusters, not the depth,  Yang Liu 06/16/06

    //void LY_modified_analyse_one_avalanche(bool mark, int seed, char* filename);   // use in the Hysteresis code
    int LY_modified_analyse_one_avalanche(bool &mark, int seed, char* filename);   // use in the Hysteresis code
    int LY_modified_analyse_one_avalanche(bool &mark, int seed, char* filename,double H);   // use in the AIGS code
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    void LY_analyse_one_small_avalanche(bool &mark, int seed, int size, char* fname, double H); //for Hysteresis Code 03/03/06
    // those avalanches with size <=3, the structure information is unique and can be precalculated

    //for AIGS Code 06/26/06, when the number of flipping spin is just one, of course there is just one avalanche.
    //But if the number >=2, there could be multiple avalanches if the resolution is too low.
    int LY_analyse_onespin_avalanche(bool &mark, int seed, char* fname, double H); 



private:
    char* spin;                                 // the actual data -- need only a bit each, but easier here    
    
    int D;               // dimension, added by Yang Liu 07/21/05
    int W;                                                           // dimensions of this window
    int L;
    int H;          
    int Z;	         // Z=2D
    int N;               // total number of spins in the system (lattice), added by Yang Liu 07/21/05   

    int* label;
    int startingSpinLoc; // starting location when scanning the avalanche, added by Yang Liu 07/21/05

    int* stride;
    int* size; // added by Y.L. 02/25/06
    mutable int* neighborLocs;
    mutable int** neighborCoords;
    mutable int* pseudoneighborLocs;
    mutable int** pseudoneighborCoords; 
    /////////////////////////////////////////////////////////////////////////////////////////////

};



#endif /* !DROPLETANALYSE_H */










