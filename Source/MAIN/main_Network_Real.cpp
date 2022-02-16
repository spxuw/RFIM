/*########################################################################################################################
  #  Avalanches in the Ground States of Random Field Ising Model (AIGS-RFIM)                                             #                                                                                                     #
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
const hType  RESOLUTION=(hType)1e15;// used to distinguish two arbitrary random fields in the lattice
                                    // resolution = 1/discreteness

const double EPSILON   = 1e-19;     // the precision of energy difference
const double FASTPOINT = 0.15;      // the critical density above which using the eariler solution (ES) will speed up 
                                    // the calculation of the M-H curve

int dist;                       // the random field distribution
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////




/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int main (int argc, char** argv)
{
  if(argc!=4 && argc!=5 && argc!=6 && argc!=7 && argc!=9)
    {
      cout <<"\nNumber of arguments is incorrect.";
      cout <<"\nmode = 0 :  Calculate the Ground State of Random-Bond Random-Field Ising Model (RBHMS) on a REAL network at given external field H.";
      cout <<"\nCommand format:  ./HMS 0  path file H [benchmarkgenefile]\n";

      cout <<"\nmode = 1 :  Calculate the Ground State of Random-Bond Random-Field Ising Model (RBHMS) on a REAL network at given external field H.";
      cout <<"\nCommand format:  ./HMS 0  path file H benchmarkpath benchmarkgenefile \n";
      
      cout <<"\nmode = 11 : Calculate the equilibrium M(H) curve, i.e., the Ground State evolution of RBHMS on a REAL network. The H boundaries [H1, Hb] are automatically chosen.";	
      cout <<"\nCommand format:  ./HMS 11 path file \n";

      cout <<"\nmode = 12 : Calculate the equilibrium M(H) curve, i.e., the Ground State evolution of RBHMS on a REAL network. The H boundaries [H1, Hb] are given!";	
      cout <<"\nCommand format:  ./HMS 12 path file Ha Hb \n";

      cout <<"\nmode = 21 : Calculate the equilibrium M(H) curve, i.e., the Ground State evolution of RBHMS on a REAL network for a given number (nH) of field values. The H boundaries [H1, Hb] are automatically chosen.\n";	
      cout <<"\nCommand format:  ./HMS 21 path file nH [benchmarkgenefile]\n";

      cout <<"\nmode = 211 : Calculate the equilibrium M(H) curve, i.e., the Ground State evolution of RBHMS on a REAL network for a given number (nH) of field values. The H boundaries [H1, Hb] are automatically chosen. The node weights are randomized\n";	
      cout <<"\nCommand format:  ./HMS 211 path file nH seed\n";

      cout <<"\nmode = 212 : Calculate the equilibrium M(H) curve, i.e., the Ground State evolution of RBHMS on a REAL network for a given number (nH) of field values. The H boundaries [H1, Hb] are automatically chosen. The edges are completely rewired (ER-randomization)\n";	
      cout <<"\nCommand format:  ./HMS 211 path file nH seed\n";

      cout <<"\nmode = 22 : Calculate the equilibrium M(H) curve, i.e., he Ground State evolution of RBHMS on a REAL network for a given number (nH) of field values. The H boundaries [H1, Hb] are given.\n";	
      cout <<"\nCommand format:  ./HMS 22 path file nH Ha Hb [benchmarkgenefile]\n";

      cout <<"\nmode = 23 : Calculate the equilibrium M(H) curve, i.e., he Ground State evolution of RBHMS on a REAL network for a given number (nH) of field values. The H boundaries [H1, Hb] are automatically chosen.\n";	
      cout <<"\nCommand format:  ./HMS 23 path file nH \n";

      cout <<"\nmode = 24 : Calculate the equilibrium M(H) curve, i.e., he Ground State evolution of RBHMS on a REAL network for a given number of lambda values with given lambda boundaries [lambda_1, lambda_2]. The gene sets file is given.\n";	
      cout <<"\nCommand format:  ./HMS 24 path file n_lambda lambda_1 lambda_2 benchmarkpath benchmarkgenefile \n";

      cout <<"\nmode = 31 : Calculate the Metastable State of Random-Bond Random-Field Ising Model (RBHMS) on a REAL network at given external field H.";
      cout <<"\nCommand format:  ./HMS 31 path file H [benchmarkgenefile]\n";
      
      cout <<"\nmode = 32 : Calculate the Hysteresis Loop, i.e., the Metastable State evolution of RBHMS on a REAL network. The H boundaries [H1, Hb] are automatically chosen.";	
      cout <<"\nCommand format:  ./HMS 32 path file \n";

      cout <<"\nmode = 41 : Calculate the highest scoring module using Ideker's simulated annealing algorithm for a REAL network.";	
      cout <<"\nCommand format:  ./HMS 41 path file itype\n";
      cout <<"itype=1: random initializatin; 2: pvalue-biased initialization";

      cout <<"\nmode = 42 : Calculate the highest scoring module using Jia's dmGWAS algorithm for a REAL network.";	
      cout <<"\nCommand format:  ./HMS 42 path file r\n";
      cout <<"r: score improvement parameter: Score_new > Score_old*(1+r) will be accepted.";

      exit(0);
    }

  float t0   = timer();                   // get the beginning time    

  int    mode = atoi(argv[1]);             // mode: (0) calculate the GS with field H; (1) calculate the entire M-H curve
  string path = argv[2];
  string file = argv[3];
  string networkfile = path + file;

  char parameters[256];
  sprintf(parameters,"%s", file.c_str());

  
  int N;
  int E; 
  vector<string> namelist;
  Get_N_E_NameList(networkfile, N, E, namelist);  


  if(mode==0)
    {
      double H = atof(argv[4]);
      cout<<"-------------------------------------------------------------------------------------------------------\n";
      cout<<"Spin Clusters in the Ground States of Random-Bond Random-Field Ising Model on a real network (CIGS-RBRFIM):\n";
      cout << " N="<< N << "; H=" << H  <<" RESOLUTION="<< (double)RESOLUTION << endl;
    
      MakeDirectory();
	
      State C(N);  
      C.Read_Fields_Bonds(networkfile);
      C.Get_BackgroundScore(networkfile);
  
      C.SetEffectiveField(H);                // set the effective local field with external field H
      Calculate_GS_M(C, false, NONE);        // Using Pre-Relabel Algorithm to get the Ground State    

      float t = timer() - t0;                
      cout << "\n To get the GS takes " << t << " s \n";
	
      C.CalE0();                     // calculate the internal energy of State C
      C.CalE(H);                     // calculate the total energy of State C(Hext=H)
      double E0 = C.GetE0();
      double E  = C.GetE();          // the ground state energy
      double m  = C.Getm();          // the ground state magnetization
      int Nup = C.Get_Nup();
      double nup = C.Get_nup();

    
      module Module;
      int LCC = C.CalModule(Module);
      double lcc = LCC/(double)N;


      cout << "H= " << H << ", E0= " << E0 << ", E= " << E << ", m=" << m 
	   << "; Nup= "                << Nup 
	   << "; |LCC|= "              << LCC 
	   << "; nup= "                << nup 
	   << "; lcc= "              << lcc
	   <<  ", Score= " << Module.score << ' ' << Module.size << " : " << Green;
      for(Nbl_itr p = Module.nodelist.begin(); p!= Module.nodelist.end(); p++) {
	cout << ' ' << namelist[*p];
      }
      cout << Normal << endl;

      char fname[256];       
      sprintf(fname,"./data/network-%s-H%5.1f",parameters, H);
      C.WriteGML(UP,fname);

    }// end of if mode==0



  else if(mode==1)
    {
      double H = atof(argv[4]);  
      cout<<"-------------------------------------------------------------------------------------------------------\n";
      cout<<"Spin Clusters in the Ground States of Random-Bond Random-Field Ising Model on a real network (CIGS-RBRFIM):\n";
      cout << " N="<< N << "; H=" << H  <<" RESOLUTION="<< (double)RESOLUTION << endl;
    
      char benchmarkfile[256]; 
      sprintf(benchmarkfile,"%s/%s", argv[5],argv[6]);
 
      MakeDirectory();
	
      State C(N);  
      C.Read_Fields_Bonds(networkfile);
      C.DFS();
   

      C.SetEffectiveField(H);                // set the effective local field with external field H
      Calculate_GS_M(C, false, NONE);        // Using Pre-Relabel Algorithm to get the Ground State    

      float t = timer() - t0;                
      cout << "\n To get the GS takes " << t << " s \n";
	
      C.CalE0();                     // calculate the internal energy of State C
      C.CalE(H);                     // calculate the total energy of State C(Hext=H)
      double E0 = C.GetE0();
      double E  = C.GetE();          // the ground state energy
      double m  = C.Getm();          // the ground state magnetization
      int Nup = C.Get_Nup();
      double nup = C.Get_nup();

    
      module Module;
      int LCC = C.CalLCC(Module);
      double lcc = LCC/(double)N;


      cout << "H= " << H << ", E0= " << E0 << ", E= " << E << ", m=" << m 
	   << "; Nup= "                << Nup 
	   << "; |LCC|= "              << LCC 
	   << "; nup= "                << nup 
	   << "; lcc= "              << lcc << endl;

      char fname[256];        
      sprintf(fname,"./data/network-%s-H%e",parameters, H);
      C.WriteGML(fname);
      C.WriteGML(UP,fname);
    }// end of if mode==1



     
  else if (mode==11 || mode==12)
    {
      
      cout<<"-------------------------------------------------------------------------------------------------------\n";
      cout<<"Avalanches in the Ground States of Random-Bond Random-Field Ising Model on a random network (AIGS-RBRFIM):\n";
      cout << " N="<< N << "; RESOLUTION="<< (double)RESOLUTION << endl;

      MakeDirectory(file);

      State S(N);                  
      S.Read_Fields_Bonds(networkfile);
      S.Get_BackgroundScore(networkfile);
  
      double hmin = S.Get_hmin();            
      double hmax = S.Get_hmax();            
      double havg = S.Get_havg();            
      double H0   = -hmax;  
      double H1   = -hmin;  
      double Hx0  = -havg;                    // =Get_GrossField(C0, C1) only valid for the initial C0 and C1

  
      State C0(N);                       // ground state(GS) C0 with H=H0
      State C1(N);                       // ground state(GS) C1 with H=H1

      ////////////////////////////////////////////////////////////////////////////////////////////////////////////
      if(mode==11) {
	C0.Initialize_withdir(DOWN);  
	C0.SetE0(DOWN);                         // C0 has internal energy = -N*D + sum(hi) with all spins DOWN

	C1.Initialize_withdir(UP);  
	C1.SetE0(UP);                           // C1 has internal energy = -N*D - sum(hi) with all spins UP
      }
      else if(mode==12) {
	if(argc!=6) {cout <<"\nCommand format:  ./AIGS 12 path file H0 H1\n"; exit(0); }
	H0 = atof(argv[4]);  
	H1 = atof(argv[5]);  
	C0.SetEffectiveField(H0);                // set the effective local field with external field H
	Calculate_GS_M(C0, false, NONE);        // Using Pre-Relabel Algorithm to get the Ground State    
	C0.CalE0();   
                  
	C1.SetEffectiveField(H1);                // set the effective local field with external field H
	Calculate_GS_M(C1, false, NONE);        // Using Pre-Relabel Algorithm to get the Ground State    
	C1.CalE0();

	Hx0 = (C1.GetE0()-C0.GetE0())/(double)(C1.GetM()-C0.GetM());

      }
	
      C0.Save_State(file, 0);                 // save the spin configuration into a file with name: `0`
      C1.Save_State(file, 1);                 // save the spin configuration into a file with name: `1`
      
      LineSegment LS0(C0,0);                  // save C0's macro information as linesegment LS0 with index `0`
      LineSegment LS1(C1,1);                  // save C1's macro information as linesegment LS1 with index `1`
      
      LS0.Set_Hboundaries(H0, Hx0);        // initial H-boundaries for LS0, should be (-\infty, Hx0)
      LS1.Set_Hboundaries(Hx0, H1);        // initial H-boundaries for LS1, should be (Hx0, +\infty)
      
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

      State C(N);            // allocate memory for state C

      
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

	  if ((deltaM-2) !=0 && (deltaM+2) != 0) // if at least one spins flips between Ck and Cl, 
	    {
	      //////////////  Find the Ground State at the crossing field Hx. Begin //////////////////////////////
	      if((*p_LSk).Get_nup() > FASTPOINT)              
		{	    
		  uES = true; FZdir = UP;
		  C.Initialize_fromfile(file, k);           // copy state C with Ck(from a file) with UP spins frozen
		}
	      else if((*p_LSl).Get_ndn() > FASTPOINT)
		{
		  uES = true; FZdir = DOWN;
		  C.Initialize_fromfile(file, l);           // copy state C with Cl(from a file) with DOWN spins frozen
		}
	      else
		{
		  uES = false; FZdir = NONE;
		  C.Initialize_fromfield(Hx);         // state C with initial all spins DOWN if Hx<0 (UP if Hx>0)
		}

	      C.SetEffectiveField(Hx);       // set the effective local field with external field Hx
	      Calculate_GS_M(C, uES, FZdir); // Using Pre-Relabel Algorithm to get the Ground State 


   
	      C.CalE0();                     // calculate the internal energy of State C
	      C.CalE(Hx);                    // calculate the total energy of State C(Hext=Hx)
	      //////////////  Find the Ground State at the crossing field Hx. End  //////////////////////////////

	      //////////////  Compare the total energy of state G and Ck(Cl) at Hx. Begin////////////////////////
	      (*p_LSk).CalE(Hx);             // calculate the total energy of State Ck(Hx), which should be 

	      if( (*p_LSk).GetE()-C.GetE() >= EPSILON &&  // Ek > E will cause sth/0=inf problem
		  (*p_LSk).GetM() != C.GetM() &&          // I think the strong condition: Mk != M != Ml 
		  (*p_LSl).GetM() != C.GetM() )           // will avoid the sth/0=inf problem
		{                                                                         
		  C.Save_State(file, sizeofList);               // save the spin configuration into a file: `sizeofList`	
		  LineSegment LS(C,sizeofList);           // get the GS C's macro information, saved as linesegment LS

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
      ofstream outf1(fname1, ios::out);
      ///////////////////////////////////////////////////////////////////////////////////////////

      cout << " Now do the Data Analysis ......\n"; 
      t0 = timer();

      List.sort(LineSegment()); // sort the list according to the Hlower value of each LineSegment


      double Hlower, Hupper, m; // LineSegment's information
      list<int> avalanche;      // list of locs of those spins participating in the avalanche
      int wholesize;      // size of the avalanche (should be equal to the whole size of those flipping spins)
      // i.e. there is one and only one avalanche during each M jump
      bool mark=true;           // if mark=true, then write a comment line in the data file to record the value of 
      // the seed and the named of all the recorded quantities  

      // The state X is used to store the avalanche information (color) during the whole M(H) process
      State X(N);                         // allocate memory for state X
      int avalancheindex = 0;                      // avalanche index

      p = List.begin();
      for(int j=1; j<=sizeofList; j++)  
	{
	  LineSegment LSa = *p;
	  State A(N);                         // allocate memory for state A
	  A.Initialize_fromfile(file, LSa.Get_index());      // copy state A from a file with name given by LSa's index 

	  Hlower = LSa.Get_Hlower();
	  Hupper = LSa.Get_Hupper();
	  m      = LSa.Getm();

	  if(mark)
	    outf1 << "# H               m\n";
	  outf1.width(15); outf1 << Hlower << ' ';
	  outf1.width(15); outf1 << m       << endl;
	  outf1.width(15); outf1 << Hupper << ' ';
	  outf1.width(15); outf1 << m       << endl;
	  
	  if(j<sizeofList)
	    {
	      LineSegment LSb = *(++p);       // ATTENTION PLEASE: ++p or p++ somewhere, don't do them both!! 
	      State B(N);            // allocate memory for state B
	      B.Initialize_fromfile(file, LSb.Get_index());  // copy state B from a file with name given by LSb's index 

	      Cal_Avalanche(A,B,avalanche);

	      wholesize = (B.GetM() - A.GetM())/2;

	      //draw this avalanche in the state X with index: avalancheindex 
	      avalancheindex++;
	    }

	}

      ///////////////////////////////////////////////////////////////////////////////////////////
      outf1.close();
      t = timer() - t0;                
      cout << " Data Analysis takes " << t << " s \n";

      cout << " \n Delete temporary ground-state-spin-configuration files.\n";
      DeleTempfiles(file);

      
      char cmd[256];
      sprintf(cmd,"./showps ./ %s &",parameters);
      system(cmd);

      
    } // end of if mode==1



  else if (mode==21 || mode==22 || mode==211|| mode==212) {
  
    cout<<"-------------------------------------------------------------------------------------------------------\n";
    cout<<"Calculate the Ground States of Random-Bond Random-Field Ising Model on a random network at a series of field values :\n";
    cout << " N="<< N << "; RESOLUTION="<< (double)RESOLUTION << endl;

    MakeDirectory();

   char fname1[256]; 


    State C0(N);  
    C0.Read_Fields_Bonds(networkfile);

    if (mode==211) { // perform randomization of node weights (random fields)
      cout << "Randomize the local fields.\n";
      C0.Randomize_Fields(1);
    }
    else if (mode==212) {
      cout << "Rewire the edges completely.\n";
      C0.Randomize_Bonds(1);
    }


    if(argc==6 || argc==8)
      sprintf(fname1,"./data/MH-%s-%s",parameters,argv[argc-1]);
    else
      sprintf(fname1,"./data/MH-%s",parameters);
    
    


    ofstream fout(fname1, ios::out);
    fout  << "# H               m\n";
    cout  << "# H               m\n";

    double nH = 10;
    double Ha=-1;
    double Hb=1;

    if(mode==21 || mode==211 || mode ==212) {
      double hmin = C0.Get_hmin();            
      double hmax = C0.Get_hmax();            

      nH = atoi(argv[4]);  
      Ha = -hmax * 1.1;
      Hb = -hmin * 1.1;

    }
    else if (mode==22) {
      nH = atof(argv[4]);  
      Ha = atof(argv[5]);
      Hb = atof(argv[6]);
    }
    //char fname2[256]; 

    double dH = 0.02;

    double m0=-1;
    double dm_max=-1;
    double Hc = Ha;

    for(double H=Ha; H <= Hb; H+=dH) { 

      State C(N);
      C.SetEffectiveField(H);                // set the effective local field with external field H
      Calculate_GS_M(C, false, NONE);        // Using Pre-Relabel Algorithm to get the Ground State    


      double m = C.Getm();          // the ground state magnetization
      cout << "H= " << H << "; m= "  << m << "; n_up= " << (1+m)/2.0 << endl;
      fout          << H << ' '      << m << endl;

      double dm= m;
      if((1+m)/2.0<0.05) { 
	dm_max = dm;
	Hc = H;
	m0 = m;
      }
    if ((1+m)/2.0>0.05){
      break;
    }
    } // end of increasing H 

    fout.close();

    float t = timer() - t0;                
    cout << "\n The rough ramp of the M-H curve takes " << t << " s \n";

    char cmd[256];
    sprintf(cmd,"./showps ./ %s &",parameters);
    system(cmd);


    cout << "\n The largest avlanche occurs at H=Hc= " << Hc << " with size= " << dm_max*0.5 << " \n";
    char benchmarkfile[256]; 
    sprintf(benchmarkfile,"%s/%s", argv[5],argv[6]);

    C0.Get_BackgroundScore(networkfile);
    C0.ReadGeneSetsFile(benchmarkfile);
    C0.DFS();
    C0.SetEffectiveField(Hc);   
    Calculate_GS_M(C0, false, NONE); 
    
    C0.CalE0();                     // calculate the internal energy of State C
    C0.CalE(Hc);                     // calculate the total energy of State C(Hext=H)
    double E0 = C0.GetE0();
    double E  = C0.GetE();          // the ground state energy
    double m  = C0.Getm();          // the ground state magnetization
    int Nup = C0.Get_Nup();
    double nup = C0.Get_nup();

    module Module;
    int LCC = C0.CalModule(Module);
    double lcc = LCC/(double)N;

   char fname[256];       
   sprintf(fname,"./data/module-%s-H%e.txt",parameters, Hc);
   ofstream fout2(fname,ios_base::out);

    cout << "H= " << Hc << ", E0= " << E0 << ", E= " << E << ", m=" << m 
	 << "; Nup= "                << Nup 
	 << "; |LCC|= "              << LCC 
	 << "; nup= "                << nup 
	 << "; lcc= "              << lcc
	 <<  ", Score= " << Module.score << ' ' << Module.size << " : " << Green;
    for(Nbl_itr p = Module.nodelist.begin(); p!= Module.nodelist.end(); p++) {
      cout << ' ' << namelist[*p];
      fout2 << namelist[*p] << endl;
    }
    cout << Normal << endl;
    fout2.close();
 

    sprintf(fname,"./data/network-%s-H%e",parameters, Hc);
    C0.WriteGML(UP,fname);
 	
  }// end of mode==2


  
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
 else if (mode==23) {
  
    cout<<"-------------------------------------------------------------------------------------------------------\n";
    cout<<"Calculate the Ground States of Random-Bond Random-Field Ising Model on a random network at a series of field values :\n";
    cout << " N="<< N << "; RESOLUTION="<< (double)RESOLUTION << endl;
 
    MakeDirectory();

    State C0(N);  
    C0.Read_Fields_Bonds(networkfile);
    C0.DFS();
 
    char fname1[256]; 

    sprintf(fname1,"./data/MH-%s",parameters);

    ofstream fout(fname1, ios::out);

    int nH = 10;
    double Ha=-1;
    double Hb=1;

    double hmin = C0.Get_hmin();            
    double hmax = C0.Get_hmax();            
    
    nH = atoi(argv[4]);  
    Ha = -hmax * 1.1;
    Hb = -hmin * 1.1;
    
    double dH = (Hb-Ha)/(double)nH;
    cout << "Ha = " << Ha << endl;
    cout << "Hb = " << Hb << endl;
    cout << "dH = " << dH << endl;

    fout  << "# H               m\n";
    cout  << "# H               m\n";

  double m0=-1;
    double dm_max=-1;
    double Hc = Ha;

    for(double H=Ha; H <= Hb; H+=dH) { 

      State C(N);
      C.SetEffectiveField(H);                // set the effective local field with external field H
      Calculate_GS_M(C, false, NONE);        // Using Pre-Relabel Algorithm to get the Ground State    
      double m = C.Getm();          // the ground state magnetization
      int Nup = C.Get_Nup();
      double nup = C.Get_nup();

      module Module;
      int LCC = C.CalLCC(Module);
      double lcc = LCC/(double)N;


      cout << Normal;
      cout << "H= "                    << H 
	   << "; m= "                  << m 
	   << "; Nup= "                << Nup 
	   << "; |LCC|= "              << LCC 
	   << "; nup= "                << nup 
	   << "; lcc= "              << lcc

	   << endl;

      fout <<  H  << ' '                 // 1
	   <<  m  << ' '                 // 2
	   << Nup << ' '                 // 3
	   << LCC << ' '                 // 4
	   << nup << ' '                 // 5
	   << lcc << ' '                 // 6

	   << endl;

      double dm= m-m0;
      if(dm>dm_max) { 
	dm_max = dm;
	Hc = H-dH;
	m0 = m;
      }
        
 
    } // end of increasing H 
    fout.close();

    float t = timer() - t0;                
    cout << "\n The rough ramp of the M-H curve takes " << t << " s \n";


    char cmd[256];
    //    sprintf(cmd,"./showps.MH         ./data-%s %s %s &",file.c_str(),parameters,argv[6]);
    sprintf(cmd,"./showps  ./  %s &",parameters);
    system(cmd);
	
    
    cout << "\n The largest avlanche occurs at H=Hc= " << Hc << " with size= " << dm_max*0.5 << " \n";
    C0.SetEffectiveField(Hc);   
    Calculate_GS_M(C0, false, NONE); 
    
    C0.CalE0();                     // calculate the internal energy of State C
    C0.CalE(Hc);                     // calculate the total energy of State C(Hext=H)
    double E0 = C0.GetE0();
    double E  = C0.GetE();          // the ground state energy
    double m  = C0.Getm();          // the ground state magnetization
    int Nup = C0.Get_Nup();
    double nup = C0.Get_nup();

    module Module;
    int LCC = C0.CalLCC(Module);
    double lcc = LCC/(double)N;
    
  
    cout << "H= " << Hc << ", E0= " << E0 << ", E= " << E << ", m=" << m 
	 << "; Nup= "                << Nup 
	 << "; |LCC|= "              << LCC 
	 << "; nup= "                << nup 
	 << "; lcc= "              << lcc

	 << endl;
 
    char fname[256];        

    sprintf(fname,"./data/network-%s-H%e",parameters, Hc);
    C0.WriteGML(fname);
    C0.WriteGML(UP,fname);

  }// end of mode==23


  ///////////////////////////////////////////////////////////////////////////////////////////////////////
 /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// same as 23, but without calculating the module scores.
 else if (mode==231) {
  
    cout<<"-------------------------------------------------------------------------------------------------------\n";
    cout<<"Calculate the Ground States of Random-Bond Random-Field Ising Model on a random network at a series of field values :\n";
    cout << " N="<< N << "; RESOLUTION="<< (double)RESOLUTION << endl;
    char benchmarkfile[256]; 
    sprintf(benchmarkfile,"%s/%s", argv[5],argv[6]);
 
    MakeDirectory();

    State C0(N);  
    C0.Read_Fields_Bonds(networkfile);
    C0.Get_BackgroundScore(networkfile);
    //int Nbgs = C0.ReadGeneSetsFile(benchmarkfile);
    C0.DFS();
 
    char fname1[256]; 
    //sprintf(fname1,"./data-%s/MH-%s-%s",file.c_str(),parameters,argv[6]);
    sprintf(fname1,"./data/MH-%s-%s",parameters,argv[6]);

    ofstream fout(fname1, ios::out);
    fout  << "# H               m\n";
    cout  << "# H               m\n";

    int nH = 10;
    double Ha=-1;
    double Hb=1;

    double hmin = C0.Get_hmin();            
    double hmax = C0.Get_hmax();            
    
    nH = atoi(argv[4]);  
    Ha = -hmax * 1.1;
    Hb = -hmin * 1.1;
    
 
    double dH = (Hb-Ha)/(double)nH;

  double m0=-1;
    double dm_max=-1;
    double Hc = Ha;

    for(double H=Ha; H <= Hb; H+=dH) { 

      State C(N);
      C.SetEffectiveField(H);                // set the effective local field with external field H
      Calculate_GS_M(C, false, NONE);        // Using Pre-Relabel Algorithm to get the Ground State    
      double m = C.Getm();          // the ground state magnetization
      int Nup = C.Get_Nup();
      double nup = C.Get_nup();

      cout << Normal;
      cout << "H= "                    << H 
	   << "; m= "                  << m 
	   << "; Nup= "                << Nup 
	   << "; nup= "                << nup 
	   << endl; 

      fout <<  H  << ' '                 // 1
	   <<  m  << ' '                 // 2
	   << Nup << ' '                 // 3
	   << nup << ' '                 // 5
	 << endl;
 
      double dm= m-m0;
      if(dm>dm_max) { 
	dm_max = dm;
	Hc = H-dH;
	m0 = m;
      }
        
 
    } // end of increasing H 
    fout.close();

    float t = timer() - t0;                
    cout << "\n The rough ramp of the M-H curve takes " << t << " s \n";
   
    cout << "\n The largest avlanche occurs at H=Hc= " << Hc << " with size= " << dm_max*0.5 << " \n";
    C0.SetEffectiveField(Hc);   
    Calculate_GS_M(C0, false, NONE); 
    
    C0.CalE0();                     // calculate the internal energy of State C
    C0.CalE(Hc);                     // calculate the total energy of State C(Hext=H)
    double E0 = C0.GetE0();
    double E  = C0.GetE();          // the ground state energy
    double m  = C0.Getm();          // the ground state magnetization
    int Nup = C0.Get_Nup();
    double nup = C0.Get_nup();

    module Module;
    int LCC = C0.CalModule(Module);
    double lcc = LCC/(double)N;

   char fname[256];       
 sprintf(fname,"./data/module-%s-H%e.txt",parameters, Hc);
   ofstream fout2(fname,ios_base::out);

    cout << "H= " << Hc << ", E0= " << E0 << ", E= " << E << ", m=" << m 
	 << "; Nup= "                << Nup 
	 << "; |LCC|= "              << LCC 
	 << "; nup= "                << nup 
	 << "; lcc= "              << lcc
	   << "; HighestModuleScore= " << Module.score  
	   << " with size= "           <<  Module.size 
	 << endl;
  
   cout << "Highestscoring module: {" << Green; 
    for(Nbl_itr p = Module.nodelist.begin(); p!= Module.nodelist.end(); p++) {
      cout << ' ' << namelist[*p];
      fout2 << namelist[*p] << endl;
    }
    cout << Normal << "}\n";
    fout2.close();
 
    sprintf(fname,"./data/network-%s-H%e",parameters, Hc);
    C0.WriteGML_HightlightBenchmarkGenes_NearestNeighbors(UP, 0, fname);


  }// end of mode==231


  ////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////
 else if (mode==24) {
  
    cout<<"-------------------------------------------------------------------------------------------------------\n";
    cout<<"Calculate the Ground States of Random-Bond Random-Field Ising Model on a random network at a series of lambda values :\n";
    cout << " N="<< N << "; RESOLUTION="<< (double)RESOLUTION << endl;
    char benchmarkfile[256]; 
    sprintf(benchmarkfile,"%s/%s", argv[7],argv[8]);
 
    MakeDirectory();

    State C0(N);  
    C0.Read_Fields_Bonds(networkfile);
    C0.Get_BackgroundScore(networkfile);
    int Nbgs = C0.ReadGeneSetsFile(benchmarkfile);
    C0.DFS();
 
    char fname1[256]; 
    //sprintf(fname1,"./data-%s/MH-%s-%s",file.c_str(),parameters,argv[6]);
    sprintf(fname1,"./data/MH-%s",parameters);

    ofstream fout(fname1, ios::out);
    fout  << "# lambda               m\n";
    cout  << "# lambda               m\n";

    int n_lambda = 10;
    double lambda_a=0;
    double lambda_b=10;
  
    n_lambda = atoi(argv[4]);  
    lambda_a = atof(argv[5]);
    lambda_b = atof(argv[6]);
    
 
    double d_lambda = (lambda_b-lambda_a)/(double)n_lambda;
    for(double lambda=lambda_a; lambda <= lambda_b; lambda+=d_lambda) { 

      State C(N);
      C.SetEffectiveField_typeII(lambda);                // set the effective local field with external field H
      Calculate_GS_M(C, false, NONE);        // Using Pre-Relabel Algorithm to get the Ground State    
      double m = C.Getm();          // the ground state magnetization
      int Nup = C.Get_Nup();
      double nup = C.Get_nup();

      module Module;
      int LCC = C.CalModule(Module);
      double lcc = LCC/(double)N;
 
      vector<double> P;
      C.FisherExactTest(P);

      cout << Normal;
      cout << "lambda= "                    << lambda
	   << "; m= "                  << m 
	   << "; Nup= "                << Nup 
	   << "; |LCC|= "              << LCC 
	   << "; nup= "                << nup 
	   << "; lcc= "              << lcc
	   << "; HighestModuleScore= " << Module.score  
	   << " with size= "           <<  Module.size 
	   << ". Enrichment analysis: p-values= (";  
      for(unsigned int a=0;a<P.size();a++) 
	cout << P[a] << ',';
      cout << "); "; 
      cout << "Highestscoring module: {" << Green; 
      for(Nbl_itr p = Module.nodelist.begin(); p!= Module.nodelist.end(); p++) 
	cout << namelist[*p] << ' ';
      cout << Normal << "}\n"; 


      fout <<  lambda << ' '                 // 1
	   <<  m  << ' '                 // 2
	   << Nup << ' '                 // 3
	   << LCC << ' '                 // 4
	   << nup << ' '                 // 5
	   << lcc << ' '                 // 6
	   << Module.score  << ' '       // 7 
	   << Module.size << ' ';        // 8
      for(unsigned int a=0; a<P.size(); a++)           
	fout << P[a] << ' ';
      for(Nbl_itr p = Module.nodelist.begin(); p!= Module.nodelist.end(); p++) 
	fout << namelist[*p] << ' ';
      fout << endl;
     
 
    } // end of increasing H 
    fout.close();

    float t = timer() - t0;                
    cout << "\n The rough ramp of the M-lambda curve takes " << t << " s \n";

    // output the benchmark gene set's minimum p-value and the corresponding H value; the minimum p-value (for active LCC only) and the corresponding H value
    char fname2[256]; 
    sprintf(fname2,"./data/GeneSets-%s-%s", parameters,argv[6]);
    C0.WriteGeneSetsFile(fname2);


    char cmd[256];
    sprintf(cmd,"./showps  ./  %s &",parameters);
    system(cmd);



    sprintf(cmd,"./showps.enrichment ./data %s %s %d &",parameters,argv[6],Nbgs);
    system(cmd);
	
  }// end of mode==24




  else if(mode==31)
    {
      double H = atof(argv[4]);  
      cout<<"-------------------------------------------------------------------------------------------------------\n";
      cout<<"Spin Clusters in the Metastable States of Random-Bond Random-Field Ising Model on a real network (CIMS-RBRFIM):\n";
      cout << " N="<< N << "; H=" << H  << endl;
    
      MakeDirectory();
	
      State C(N);  
      C.Read_Fields_Bonds(networkfile);
      C.Get_BackgroundScore(networkfile);
  
      C.CalMetastableState(H);

      float t = timer() - t0;                
      cout << "\n To get the MS takes " << t << " s \n";
	
      C.CalE0();                     // calculate the internal energy of State C
      C.CalE(H);                     // calculate the total energy of State C(Hext=H)
      double E0 = C.GetE0();
      double E  = C.GetE();          // the total energy
      double m  = C.Getm();          // the magnetization

   
      module Module;
      C.CalModule(Module);

      cout << "H= " << H << ", E0= " << E0 << ", E= " << E << ", m=" << m << ", Score= " << Module.score << ' ' << Module.size << " : " << Green;
      for(Nbl_itr p = Module.nodelist.begin(); p!= Module.nodelist.end(); p++) {
	cout << ' ' << namelist[*p];
      }
      cout << Normal << endl;

      char fname[256]; 
      sprintf(fname,"./data/%s-H%e",parameters, H);
      C.WriteGML(fname);
      C.WriteGML(UP,fname);
      C.WriteGML(UP,3, fname);

    }// end of if mode==31



  else if (mode==32) {
  
    cout<<"-------------------------------------------------------------------------------------------------------\n";
    cout<<"Calculate the Hysteresis Loop of Random-Bond Random-Field Ising Model on a random network :\n";
    cout << " N="<< N << endl;

    MakeDirectory(file);

    State C(N);  
    C.Read_Fields_Bonds(networkfile);
    C.Get_BackgroundScore(networkfile);

    char fname1[256]; 
    sprintf(fname1,"./data/MH-%s",parameters);

    C.HysteresisLoop(fname1);

    float t = timer() - t0;                
    cout << "\n The whole M-H curve takes " << t << " s \n";

    char cmd[256];
    sprintf(cmd,"./showps ./ %s &",parameters);
    system(cmd);
  }// end of mode==32




 else if (mode==41) {
  
    cout<<"-------------------------------------------------------------------------------------------------------\n";
    cout <<"Calculate the highestscoring module using Ideker's simulated annealing algorithm for a REAL network:\n";	
    cout << " N="<< N << endl;

    MakeDirectory(file);

    State C(N);  
    C.Read_Fields_Bonds(networkfile);
    C.Get_BackgroundScore(networkfile);

    char fname1[256]; 
    sprintf(fname1,"./data/jActive-%s",parameters);

    int itype = atoi(argv[4]);
    double Ti = 1.0;
    double Tf = 1e-5;
    int seed=1;
    C.HMS_Ideker(itype, Ti, Tf, seed, fname1);  

    float t = timer() - t0;                
    cout << "\n The whole jActive takes " << t << " s \n";

  }// end of mode==41





 else if (mode==42) {
  
    cout<<"-------------------------------------------------------------------------------------------------------\n";
    cout <<"Calculate the highestscoring module using Jia's dmGWAS algorithm for a REAL network:\n";	
    cout << " N="<< N << endl;

    MakeDirectory(file);

    State C(N);  
    C.Read_Fields_Bonds(networkfile);
    C.Get_BackgroundScore(networkfile);

    char fname1[256]; 
    sprintf(fname1,"./data/dmGWAS-%s",parameters);

    double r = atof(argv[4]);
    C.HMS_Jia(r, fname1);

    float t = timer() - t0;                
    cout << "\n The whole dmGWAS takes " << t << " s \n";

  }// end of mode==4


  cout << "\n Done!  &:) \n";
  cout << "--------------------------------------------------------------------------------------------------------\n";
  exit(0);
}
////////////////////////////////       THE END      /////////////////////////////////////////////////////////////////////





