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
#include "nr.h"                
#include "Point.h"             
#include "State.h"             
using namespace std;

class Window
{
 public:
    Window(std::istream& x);                           
    Window(const Window& x);                        
    Window(int dimension, int length, list<int>& spinList); 
    Window(int dimension, int length, char* spinCluster); 
    Window(const State& x);                         

    ~Window() 
	{
	    delete[] spin; 
	    if (label != NULL) 
		delete[] label;
	    
	    delete [] stride;
	    delete [] size;
	    delete [] neighborLocs;
	    for(int i=0;i<Z;i++)
	    {
		delete [] neighborCoords[i];
		delete [] pseudoneighborCoords[i];
	    }
	    delete [] neighborCoords;
	    delete [] pseudoneighborCoords;
	};
  

    
    void printme();
    void geomview(char* fname, int plane);             
    void geomviewsurf(char* fname);
    void postlayer(char* fname);                       
    int* slicecount(int dir);
    int* slicecountcenter(int dir);
   
    int wallarea(int startboundary);  
   
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
    

    
    int label_connected_cluster(int clustindex, int newlabel, char s, int* surface, int* boxvol, int* maxdim);
    void report_connected_clusters(int* bigidx, int* bigvol);
    int domain_with_notlabel(int clustindex, int notlabel, int* surface, int* boxvol, int* maxdim);
    void report_domains_with_notlabel(int notlabel);
    void print_bubble_surfaces();
    


    
    
    void GetNeighbors(int currentLoc) const; 
    void GetNeighbors(const Point point) const;
    void GetNeighbors(int currentLoc, const int* currentCoords) const;
    void GetPseudoNeighbors(const int* currentCoords) const;
    void GetPseudoNeighbors(const Point point) const;
    void GetCoords(int loc, int* newCoords) const;
    int  GetLoc(const int* coords) const;
    int  GetLoc_PBC(const int* coords) const; 
    int  GetLoc_PBC(const Point point) const; 
    void get_inertia_tensor();
    Vec_DP get_shape_tensor();
    Vec_DP get_eigenvalue(list<Point>  pointList);
    Vec_DP get_eigenvalue(list<Point>  pointList, int size); 
    



    
    
    
    int  LY_labelCC(int clustindex, int newlabel, char s, int* surface, int* boxvol, int* maxdim);
    int  LY_labelCC(int clustindex, int newlabel, char s, int* surface, int* boxvol, int* maxdim, bool* Touch);
    int  LY_labelCC_getQ(int loc, int newlabel, char s, int* surface, int* boxvol, int* maxdim, Vec_DP* EV);
    int  LY_labelCC_getQ(int loc, int newlabel, char s, int* surface, int* boxvol, int* maxdim, Vec_DP* EV, bool* Touch);

    int  LY_labelCC(int loc, int newlabel, char s, bool* BoundaryTouching);  
    int  LY_labelCC_getQ(int loc, int newlabel, char s, 
			 int* surface, int* boxvol, int* maxdim, 
			 Vec_DP* EV, bool* Touch,
			 list<Point>* newList, vector<int>* newSize);


    
    void LY_reportCC(int* bigidx, int* bigvol);
    void LY_modified_reportCC(int* bigidx,int* bigvol);

    
    int LY_modified_reportCC(int* bigidx, int* bigvol, char* filename) ; 

    void LY_modified_reportCC(int* bigidx,int* bigvol,
			      int* v,int* a,int* bv,int* md,Vec_DP* EV); 
    void LY_modified_reportCC(int* bigidx,int* bigvol,
			      int* v,int* a,int* bv,int* md, Vec_DP* EV,
			      list<Point>* newList, vector<int>* newSize); 
    void LY_modified_reportCC(int* bigidx,int* bigvol,
			      int* v,int* a,int* bv,int* md, Vec_DP* EV,
			      vector<int>* Holeinfo); 

    
    int LY_modified_reportCC(int* bigidx, int* bigvol, bool mark, int seed, char* filename, double H);    
    


    
    
    int  LY_domain_with_notlabel(int loc, int notlabel, int* surface, int* boxvol, int* maxdim);
    int  LY_domain_with_notlabel(int loc, int notlabel, int* surface);
    int  LY_domain_with_notlabel(int loc, int notlabel); 

    
    void LY_reportD_with_notlabel(int notlabel);  
    void LY_reportD_with_notlabel(int depth, int notlabel,int* v_D, int* a_D, int* bv_D, int* md_D);
    void LY_reportD_with_notlabel(int depth, int notlabel,int* v_D, int* a_D);
    void LY_reportD_with_notlabel(int depth, int notlabel, char* filename);  
    void LY_reportD_with_notlabel(int depth, int notlabel, bool mark, int seed, char* filename, double H); 


    
    void LY_print_bubble_surfaces();
    void LY_analyse_one_avalanche();
    void LY_analyse_one_avalanche(char* filename);   
    int  LY_analyse_many_avalanches(char* filename); 
    void LY_analyse_one_avalanche(bool mark, int seed, char* filename);   

    int  LY_analyse_many_avalanches(bool mark, int seed, char* mainname, double H); 
    

    
    int LY_modified_analyse_one_avalanche(bool &mark, int seed, char* filename);   
    int LY_modified_analyse_one_avalanche(bool &mark, int seed, char* filename,double H);   
    

    void LY_analyse_one_small_avalanche(bool &mark, int seed, int size, char* fname, double H); 
    

    
    
    int LY_analyse_onespin_avalanche(bool &mark, int seed, char* fname, double H); 



private:
    char* spin;                                 
    
    int D;               
    int W;                                                           
    int L;
    int H;          
    int Z;	         
    int N;               

    int* label;
    int startingSpinLoc; 

    int* stride;
    int* size; 
    mutable int* neighborLocs;
    mutable int** neighborCoords;
    mutable int* pseudoneighborLocs;
    mutable int** pseudoneighborCoords; 
    

};



#endif










