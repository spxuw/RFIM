/*######################################################################################################################
  #  Avalanches in the Ground States of Random Field Ising Model (AIGS-RFIM)                                             #
  #  Implementation of Vives's algorithm.                                                                                #
  #  Yang Liu                                                                                                            #
  #  06/26/06                                                                                                            #
  #######################################################################################################################*/
#include <iostream>
#include <list>
#include <queue>
#include "types.h"   
#include "hi_pr.h" 
#include "State.h"   
#include "LineSegment.h"
#include "CrossingField.h"               
#include "DropletAnalyse.h"    
#include "MakeDirectory.h"

using namespace std;
//////////////////////////////////      GLOBAL CONSTANT          ///////////////////////////////////////////////////////
const hType  RESOLUTION=(hType)1e17;// used to distinguish two arbitrary random fields in the lattice
                                    // resolution = 1/discreteness

const double EPSILON   = 1e-19;     // the precision of energy difference
const double FASTPOINT = 0.15;      // the critical density above which using the eariler solution (ES) will speed up 
                                    // the calculation of the M-H curve

int dist;                       // the random field distribution
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



///////////////////////////    For Hypercubic Lattice with periodic boundary condition and random fields  //////////////
int main (int argc, char** argv)
{
  if(argc!=7 && argc!=8)
    {
      cout <<"\nNumber of arguments is incorrect.";
      cout <<"\nCommand format:  ./AIGS mode D L R dist seed (H)";
      cout <<"\nmode = 0 : Calculate the Ground States of Random Field Ising Model (CIGS-RFIM) at given H.";
      cout <<"\nmode = 1 : Calculate the Avalanches in the Ground State evolution of the Random Field Ising Model (AIGS-RFIM).\n";	exit(0);
    }

  float t0   = timer();                   // get the beginning time    

  int mode   = atoi(argv[1]);             // mode: (0) calculate the GS with field H; (1) calculate the entire M-H curve
  int D      = atoi(argv[2]);             // dimensionality
  int L      = atoi(argv[3]);             // linear size of the hypercubic lattice
  double R   = atof(argv[4]);             // disorder parameter(standard deviation of the Gaussian)
  dist   = atoi(argv[5]);             // distribution: 0 : bimodal 1: gaussian 2: rectangular 
  int seed   = atoi(argv[6]);             // random seed

  int N=1;
  for(int i=0;i<D;i++) 
    N*=L;
  double c = 2*D;

  char parameters[256];
  sprintf(parameters,"D%d-L%d-R%.3lf-dist%d-seed%d", D, L, R, dist, seed);


  if(mode==0)
    {

      double H = atof(argv[7]);  
      cout<<"-------------------------------------------------------------------------------------------------------\n";
      cout<<"Spin Clusters in the Ground States of Random Field Ising Model (CIGS-RFIM):\n";
      cout << " D="<< D <<" L="<< L <<" R="<< R  << " dist=" << dist << " seed=" << seed << " H=" << H  
	   <<" RESOLUTION="<< (double)RESOLUTION << endl;
      MakeDirectory();


      //////////////  Find the Ground State at the field Hext=0. Begin //////////////////////////////
      //State C(D, L, R, seed, DOWN);          // ground state(GS) C0 with all spins DOWN, M=-N
      //C.SetRandomField(D, L, R, seed);       // set the quenched-disorder(Gaussian random field) for all States
      State C(N, R, seed);
      //C.Set_RandomField_RandomBond_Lattice(D, L, R, seed, 0); //PBC=0
      C.Set_RandomField_RandomBond_Lattice(D, L, R, seed, 1); //PBC=1


      C.SetEffectiveField(H);                // set the effective local field with external field H
      Calculate_GS_M(C, false, NONE);        // Using Pre-Relabel Algorithm to get the Ground State    

      float t = timer() - t0;                
      cout << "\n To get the GS takes " << t << " s \n";

      C.CalE0();                     // calculate the internal energy of State C
      C.CalE(H);                     // calculate the total energy of State C(Hext=H)
      double Egs = C.GetE();          // the ground state energy
      double Mgs = C.Getm();          // the ground state magnetization
      cout << "Egs= " << Egs << "  Mgs=" << Mgs << endl;

      cout << " Now do the Data Analysis ......\n"; 
      t0 = timer();

      char fname[256]; 
      sprintf(fname,"./data/%s-H%e",parameters, H);
      C.WriteGML(fname);
      C.WriteGML(UP,fname);


      t = timer() - t0;                
      cout << " Data Analysis takes " << t << " s \n";

    }// end of if mode==0

     
  else if (mode==1)
    {

      cout<<"-------------------------------------------------------------------------------------------------------\n";
      cout<<"Avalanches in the Ground States of Random Field Ising Model (AIGS-RFIM):\n";
      cout << " D="<< D <<" c=" << 2*D << " L="<< L <<" R="<< R  << " dist=" << dist << " seed=" << seed 
	   <<" RESOLUTION="<< (double)RESOLUTION << " fastpoint= " << FASTPOINT << endl;

      MakeDirectory(N,c,R,seed);

      ////////////////////////////////////////////////////////////////////////////////////////////////////////////
      //State C0(D, L, R, seed, DOWN);          // ground state(GS) C0 with all spins DOWN, M=-N
      //C0.SetRandomField(D, L, R, seed);       // set the quenched-disorder(Gaussian random field) for all States
 
      State C0(N, R, seed, DOWN);
      //C0.Set_RandomField_RandomBond_Lattice(D, L, R, seed, 0); // PBC=0
      C0.Set_RandomField_RandomBond_Lattice(D, L, R, seed, 1); // PBC=1
      C0.SetE0(DOWN);                         // C0 has internal energy = -N*D + sum(hi) with all spins DOWN

      //State C1(D, L, R, seed, UP);            // state C1 with all spins UP, M=+N
      State C1(N, R, seed, UP);
      C1.SetE0(UP);                           // C1 has internal energy = -N*D - sum(hi) with all spins UP

      C0.Save_State(0);                       // save the spin configuration into a file with name: `0`
      C1.Save_State(1);                       // save the spin configuration into a file with name: `1`

      LineSegment LS0(C0,0);                  // save C0's macro information as linesegment LS0 with index `0`
      LineSegment LS1(C1,1);                  // save C1's macro information as linesegment LS1 with index `1`

      double hmin = C0.Get_hmin();            
      double hmax = C0.Get_hmax();            
      double havg = C0.Get_havg();            
      double Hx0  = -havg;                    // =Get_GrossField(C0, C1) only valid for the initial C0 and C1
      LS0.Set_Hboundaries(-hmax, Hx0);        // initial H-boundaries for LS0, should be (-\infty, Hx0)
      LS1.Set_Hboundaries(Hx0, -hmin);        // initial H-boundaries for LS1, should be (Hx0, +\infty)

      list<LineSegment> List;
      int sizeofList = 0;
      List.push_back(LS0);                    // push into the list of LineSegment (GS's macro information)
      List.push_back(LS1);
      sizeofList = 2;

      queue<CrossingField> Queue;
      CrossingField q0(Hx0,0,1);              // State C0 and C1 have the crossing field Hx0.
      Queue.push(q0);                         // push into the queue of crossing fiels(CF)



      /////////////////////////////////////////////////////////////////////////////////////////////////////////
      // As long as the queue Q is not empty, find the successive crossing fields.
      list<LineSegment>::iterator p, p_LSk, p_LSl; 
      int i, k, l;                            // GS index
      int deltaM;                             // magnetization difference between GS Ck and Cl
      double b,Hx,B;                          // crossing fields 
      bool uES;                               // a flag whether we can use earlier solution(ES)
      char FZdir;                             // if we can use ES, then what is the frozen direction 

      State C(N,R,seed);            // allocate memory for state C

      while(!Queue.empty()) 
	{
	  CrossingField q = Queue.front();    // retrieve the next crossing field from the queue
	  Hx = q.GetHx();                                                             
	  k  = q.Getk();                      // get one GS index 
	  l  = q.Getl();                      // get another GS index
	  Queue.pop();                        // remove this crossing field from the queue
	  cout << "\n Hx= "; cout.width(10); cout << Hx << ",# of GS:" << sizeofList << ",  "; //test

	  for(p=List.begin(),i=0; p!=List.end(); p++,i++)
	    {
	      if(i==k)  p_LSk = p;            // get the pointer to the kth GS  
	      if(i==l)  p_LSl = p;            // get the pointer to the lth GS   
	    }

	  deltaM = (*p_LSk).GetM()-(*p_LSl).GetM();

	  // In Middletion's note, he suggests that one just consider big avalanches to speed up the code.
	  // But this will cause trouble, since small avalanches will be considered as big one.
	  // This size distribution must be changed. Only the size distribution of avalanches with size > Smin 
	  // can be correctly calculated.
	  // if( fabs(deltaM) > 2* Smin ) // if we just consider big avalanches with size > Smin

	  if ((deltaM-2) !=0 && (deltaM+2) != 0) // if at least one spins flips between Ck and Cl, 
	    {
	      //////////////  Find the Ground State at the crossing field Hx. Begin //////////////////////////////
	      if((*p_LSk).Get_nup() > FASTPOINT)              
		{	    
		  //cout << " with frozen UP spins !"; //test
		  uES = true; FZdir = UP;
		  C.Initialize_fromfile(k);           // copy state C with Ck(from a file) with UP spins frozen
		  //C.Show_State();//test
		}
	      else if((*p_LSl).Get_ndn() > FASTPOINT)
		{
		  //cout << " with frozen DOWN spins !"; //test
		  uES = true; FZdir = DOWN;
		  C.Initialize_fromfile(l);           // copy state C with Cl(from a file) with DOWN spins frozen
		  //C.Show_State();//test
		}
	      else
		{
		  uES = false; FZdir = NONE;
		  C.Initialize_fromfield(Hx);         // state C with initial all spins DOWN if Hx<0 (UP if Hx>0)
		  //C.Show_State();  //test
		}


	      C.SetEffectiveField(Hx);       // set the effective local field with external field Hx
	      Calculate_GS_M(C, uES, FZdir); // Using Pre-Relabel Algorithm to get the Ground State    
	      C.CalE0();                     // calculate the internal energy of State C
	      C.CalE(Hx);                    // calculate the total energy of State C(Hext=Hx)
	      //////////////  Find the Ground State at the crossing field Hx. End  //////////////////////////////

	      //////////////  Compare the total energy of state G and Ck(Cl) at Hx. Begin////////////////////////
	      (*p_LSk).CalE(Hx);             // calculate the total energy of State Ck(Hx), which should be 
	      // equal to Cl(Hx) due to  the definition of Hx.

	      if( (*p_LSk).GetE()-C.GetE() >= EPSILON &&  // Ek > E will cause sth/0=inf problem
		  (*p_LSk).GetM() != C.GetM() &&          // I think the strong condition: Mk != M != Ml 
		  (*p_LSl).GetM() != C.GetM() )           // will avoid the sth/0=inf problem
		{                                                                         
		  C.Save_State(sizeofList);               // save the spin configuration into a file: `sizeofList`	
		  LineSegment LS(C,sizeofList);           // get the GS C's macro information, saved as linesegment LS
		  // with index `sizeofList`

		  B = Get_CrossField((*p_LSk), LS);       // calculate the CF for GS Ck and C                     
		  b = Get_CrossField(LS, (*p_LSl));       // calculate the CF for GS C and Cl                      
                                                                       
		  (*p_LSk).Set_Hupper(B);                 // update H-boundaries for GS Ck 
		  (*p_LSl).Set_Hlower(b);                 // update H-boundaries for GS Cl                             

		  CrossingField q1(B,k,sizeofList);       // B is the CF for GS Ck and C(index=current size of list) 
		  CrossingField q2(b,sizeofList,l);       // b is the CF for GS C(index=current size of list) and Cl
      
		  Queue.push(q1);                         // push into the queue                            
		  Queue.push(q2);                                                    
                                                                       
		  LS.Set_Hboundaries(B, b);             // set H-boundaries for GS C                            
		  List.push_back(LS);                     // push into the list of GS                            
		  sizeofList++;                                                      

		}// end of `if E(C) is smaller`
	      C.Reset();
	      //////////////  Compare the total energy of state C and Ck(Cl) at Hx. End ////////////////////////

	    }// end of `if at least one spins flips between Ck and Cl'

	}// end of the while loop
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////


      float t = timer() - t0;                
      cout << "\n The full ramp of the M-H curve takes " << t << " s \n";



      ///////////////////////////////////////////////////////////////////////////////////////////
      char fname1[256]; 
      sprintf(fname1,"./data/MH-%s",parameters);
      //char fname2[256]; 
      //sprintf(fname2,"./data/Surface-Q-Avala-D%d-L%d-R%.3lf-dist%d",D,L,R,dist);
      ofstream outf1(fname1, ios::out);
      ///////////////////////////////////////////////////////////////////////////////////////////



      cout << " Now do the Data Analysis ......\n"; 
      t0 = timer();

      List.sort(LineSegment()); // sort the list according to the Hlower value of each LineSegment


      double Hlower, Hupper, m; // LineSegment's information
      list<int> avalanche;      // list of locs of those spins participating in the avalanche
      //int size;
      int wholesize;      // size of the avalanche (should be equal to the whole size of those flipping spins)
      // i.e. there is one and only one avalanche during each M jump
      //bool mark=true;           // if mark=true, then write a comment line in the data file to record the value of 
      // the seed and the named of all the recorded quantities  

      // The state X is used to store the avalanche information (color) during the whole M(H) process
      State X(N,R,seed);                         // allocate memory for state X
      int avalancheindex = 0;                      // avalanche index

      p = List.begin();
      for(int j=1; j<=sizeofList; j++)  
	{
	  LineSegment LSa = *p;
	  State A(N,R,seed);                         // allocate memory for state A
	  A.Initialize_fromfile(LSa.Get_index());      // copy state A from a file with name given by LSa's index 

	  Hlower = LSa.Get_Hlower();
	  Hupper = LSa.Get_Hupper();
	  m      = LSa.Getm();

	  //if(mark)
	  //outf1 << "/ seed="<<seed<<"    H               m\n";
	  outf1.width(15); outf1 << Hlower << ' ';
	  outf1.width(15); outf1 << m       << endl;
	  outf1.width(15); outf1 << Hupper << ' ';
	  outf1.width(15); outf1 << m       << endl;

	  //cout << Hlower << ' ' << m << endl; //test

	  if(j<sizeofList)
	    {
	      LineSegment LSb = *(++p);       // ATTENTION PLEASE: ++p or p++ somewhere, don't do them both!! 
	      State B(N,R,seed);            // allocate memory for state B
	      B.Initialize_fromfile(LSb.Get_index());  // copy state B from a file with name given by LSb's index 

	      Cal_Avalanche(A,B,avalanche);
	      Window W(D, L, avalanche);

	      wholesize = (B.GetM() - A.GetM())/2;

	      //draw this avalanche in the state X with index: avalancheindex 
	      //		X.DrawOneAvalanche(avalancheindex,avalanche);
	      avalancheindex++;

	      /*
	      //if we are not sure whether RESOLUTION is high enough, we can check this by checking whethere there 
	      // are multiple avalanches during each M jump.
	      if(wholesize==1)
	      size  = W.LY_analyse_onespin_avalanche(mark, seed, fname2, Hupper);
	      else
	      size  = W.LY_modified_analyse_one_avalanche(mark,seed, fname2, Hupper); 

	      if(size < wholesize)
	      cout << "At Hupper="<< Hupper 
	      <<"; # of avalanche > 1 ==> Resolution is too low! " << endl; 
	      */

	      /*
	      // if we choose pretty high RESOLUTION, then we are pretty sure that multiple-avalanches in one M jump 
	      // are impossible! 
	      if(wholesize < 3)
		{
		  size = wholesize; 
		  //cout << "wholesize=" << wholesize << endl; //test
		  W.LY_analyse_one_small_avalanche(mark,seed, size, fname2, Hupper);
		} 
	      else
		{
		  //cout << "wholesize=" << wholesize << endl; //test
		  size  = W.LY_modified_analyse_one_avalanche(mark,seed, fname2, Hupper); 
		}
	      */
	    }

	}

      // draw all the avalanches (with different colors) for postscript file.
      //	char fname[256]; 
      //	sprintf(fname,"./data/Avalanches-%s",parameters);
      //	X.DrawAllAvalanches(fname);

      ///////////////////////////////////////////////////////////////////////////////////////////
      outf1.close();
      t = timer() - t0;                
      cout << " Data Analysis takes " << t << " s \n";

      cout << " \n Delete temporary ground-state-spin-configuration files.\n";
      DeleTempfiles(N,c,R,seed);

      char cmd[256];
      sprintf(cmd,"./showps ./ %s &",parameters);
      system(cmd);

    } // end of if mode==1


  cout << "\n Done!  &:) \n";
  cout << "--------------------------------------------------------------------------------------------------------\n";
  exit(0);
}
////////////////////////////////       THE END      /////////////////////////////////////////////////////////////////////








