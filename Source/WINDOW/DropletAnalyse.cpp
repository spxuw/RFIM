#define GRAY 32
#define BLACK 64
#define NOGO 16
#define NOTWHITE 96
//////////////////////////////////////////////////////////////////////////////////////////////////////
#include "nr.h"                                                                       // Y.L. 06/27/05
#include "DropletAnalyse.h"
//////////////////////////////////////////////////////////////////////////////////////////////////////
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <list>
#include "DropletAnalyse_initpart.cpp"          //initialize part                     // Y.L. 02/27/06
#include "DropletAnalyse_GetShapeTensor.cpp"   // part about calculation of shape tensor// Y.L. 02/27/06
#include "DropletAnalyse_prolixpart.cpp"     //unused/unimportant functions are here!  // Y.L. 02/25/06

#include "DropletAnalyse_AAM.cpp"  //AAM's orginal code (bug about maxdim removed)        Y.L.02/28/06

// Here comes LY's update on AAM's original code                       02/28/06
#include "DropletAnalyse_LY_labelCC.cpp"   
#include "DropletAnalyse_LY_reportCC.cpp"   
#include "DropletAnalyse_LY_labelD.cpp"   
#include "DropletAnalyse_LY_reportD.cpp"   
#include "DropletAnalyse_LY.cpp"   


using namespace std;
 

