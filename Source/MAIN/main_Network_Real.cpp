#include <iostream>
#include <list>
#include <queue>
#include "types.h"   
#include "hi_pr.h" 
#include "State.h"   
#include "LineSegment.h"
#include "CrossingField.h"               
#include "MakeDirectory.h"
    
using namespace std;

const hType  RESOLUTION=(hType)1e15;
                                    

const double EPSILON   = 1e-19;     
const double FASTPOINT = 0.15;      
                                    

int dist;                       






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

  float t0   = timer();                   

  int    mode = atoi(argv[1]);             
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
  
      C.SetEffectiveField(H);                
      Calculate_GS_M(C, false, NONE);        

      float t = timer() - t0;                
      cout << "\n To get the GS takes " << t << " s \n";
	
      C.CalE0();                     
      C.CalE(H);                     
      double E0 = C.GetE0();
      double E  = C.GetE();          
      double m  = C.Getm();          
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

    }



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
   

      C.SetEffectiveField(H);                
      Calculate_GS_M(C, false, NONE);        

      float t = timer() - t0;                
      cout << "\n To get the GS takes " << t << " s \n";
	
      C.CalE0();                     
      C.CalE(H);                     
      double E0 = C.GetE0();
      double E  = C.GetE();          
      double m  = C.Getm();          
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
    }



     
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
      double Hx0  = -havg;                    

  
      State C0(N);                       
      State C1(N);                       

      
      if(mode==11) {
	C0.Initialize_withdir(DOWN);  
	C0.SetE0(DOWN);                         

	C1.Initialize_withdir(UP);  
	C1.SetE0(UP);                           
      }
      else if(mode==12) {
	if(argc!=6) {cout <<"\nCommand format:  ./AIGS 12 path file H0 H1\n"; exit(0); }
	H0 = atof(argv[4]);  
	H1 = atof(argv[5]);  
	C0.SetEffectiveField(H0);                
	Calculate_GS_M(C0, false, NONE);        
	C0.CalE0();   
                  
	C1.SetEffectiveField(H1);                
	Calculate_GS_M(C1, false, NONE);        
	C1.CalE0();

	Hx0 = (C1.GetE0()-C0.GetE0())/(double)(C1.GetM()-C0.GetM());

      }
	
      C0.Save_State(file, 0);                 
      C1.Save_State(file, 1);                 
      
      LineSegment LS0(C0,0);                  
      LineSegment LS1(C1,1);                  
      
      LS0.Set_Hboundaries(H0, Hx0);        
      LS1.Set_Hboundaries(Hx0, H1);        
      
      list<LineSegment> List;
      int sizeofList = 0;
      List.push_back(LS0);                    
      List.push_back(LS1);
      sizeofList = 2;

      queue<CrossingField> Queue;
      CrossingField q0(Hx0,0,1);              
      Queue.push(q0);                         
      
      
      
      list<LineSegment>::iterator p, p_LSk, p_LSl; 
      int i, k, l;                            
      int deltaM;                             
      double b,Hx,B;                          
      bool uES;                               
      char FZdir;                             

      State C(N);            

      
      while(!Queue.empty()) 
	{
	  CrossingField q = Queue.front();    
	  Hx = q.GetHx();                                                             
	  k  = q.Getk();                      
	  l  = q.Getl();                      
	  Queue.pop();                        
	  cout << "\n Hx= "; cout.width(10); cout << Hx << ",# of GS:" << sizeofList << ",  "; 

	  for(p=List.begin(),i=0; p!=List.end(); p++,i++)
	    {
	      if(i==k)  p_LSk = p;            
	      if(i==l)  p_LSl = p;            
	    }

	  deltaM = (*p_LSk).GetM()-(*p_LSl).GetM();

	  if ((deltaM-2) !=0 && (deltaM+2) != 0) 
	    {
	      
	      if((*p_LSk).Get_nup() > FASTPOINT)              
		{	    
		  uES = true; FZdir = UP;
		  C.Initialize_fromfile(file, k);           
		}
	      else if((*p_LSl).Get_ndn() > FASTPOINT)
		{
		  uES = true; FZdir = DOWN;
		  C.Initialize_fromfile(file, l);           
		}
	      else
		{
		  uES = false; FZdir = NONE;
		  C.Initialize_fromfield(Hx);         
		}

	      C.SetEffectiveField(Hx);       
	      Calculate_GS_M(C, uES, FZdir); 


   
	      C.CalE0();                     
	      C.CalE(Hx);                    
	      

	      
	      (*p_LSk).CalE(Hx);             

	      if( (*p_LSk).GetE()-C.GetE() >= EPSILON &&  
		  (*p_LSk).GetM() != C.GetM() &&          
		  (*p_LSl).GetM() != C.GetM() )           
		{                                                                         
		  C.Save_State(file, sizeofList);               
		  LineSegment LS(C,sizeofList);           

		  B = Get_CrossField((*p_LSk), LS);       
		  b = Get_CrossField(LS, (*p_LSl));       
                                                                       
		  (*p_LSk).Set_Hupper(B);                 
		  (*p_LSl).Set_Hlower(b);                 

		  CrossingField q1(B,k,sizeofList);       
		  CrossingField q2(b,sizeofList,l);       
      
		  Queue.push(q1);                         
		  Queue.push(q2);                                                    
                                                                       
		  LS.Set_Hboundaries(B, b);             
		  List.push_back(LS);                     
		  sizeofList++;                                                      

		}
	      C.Reset();
	      

	    }

	}
      


      float t = timer() - t0;                
      cout << "\n The full ramp of the M-H curve takes " << t << " s \n";

      
      char fname1[256]; 
      sprintf(fname1,"./data/MH-%s",parameters);
      ofstream outf1(fname1, ios::out);
      

      cout << " Now do the Data Analysis ......\n"; 
      t0 = timer();

      List.sort(LineSegment()); 


      double Hlower, Hupper, m; 
      list<int> avalanche;      
      int wholesize;      
      
      bool mark=true;           
      

      
      State X(N);                         
      int avalancheindex = 0;                      

      p = List.begin();
      for(int j=1; j<=sizeofList; j++)  
	{
	  LineSegment LSa = *p;
	  State A(N);                         
	  A.Initialize_fromfile(file, LSa.Get_index());      

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
	      LineSegment LSb = *(++p);       
	      State B(N);            
	      B.Initialize_fromfile(file, LSb.Get_index());  

	      Cal_Avalanche(A,B,avalanche);

	      wholesize = (B.GetM() - A.GetM())/2;

	      
	      avalancheindex++;
	    }

	}

      
      outf1.close();
      t = timer() - t0;                
      cout << " Data Analysis takes " << t << " s \n";

      cout << " \n Delete temporary ground-state-spin-configuration files.\n";
      DeleTempfiles(file);

      
      char cmd[256];
      sprintf(cmd,"./showps ./ %s &",parameters);
      system(cmd);

      
    } 



  else if (mode==21 || mode==22 || mode==211|| mode==212) {
  
    cout<<"-------------------------------------------------------------------------------------------------------\n";
    cout<<"Calculate the Ground States of Random-Bond Random-Field Ising Model on a random network at a series of field values :\n";
    cout << " N="<< N << "; RESOLUTION="<< (double)RESOLUTION << endl;

    MakeDirectory();

   char fname1[256]; 


    State C0(N);  
    C0.Read_Fields_Bonds(networkfile);

    if (mode==211) { 
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
    

    double dH = 0.02;

    double m0=-1;
    double dm_max=-1;
    double Hc = Ha;

    for(double H=Ha; H <= Hb; H+=dH) { 

      State C(N);
      C.SetEffectiveField(H);                
      Calculate_GS_M(C, false, NONE);        


      double m = C.Getm();          
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
    } 

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
    
    C0.CalE0();                     
    C0.CalE(Hc);                     
    double E0 = C0.GetE0();
    double E  = C0.GetE();          
    double m  = C0.Getm();          
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
 	
  }


  
  
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
      C.SetEffectiveField(H);                
      Calculate_GS_M(C, false, NONE);        
      double m = C.Getm();          
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

      fout <<  H  << ' '                 
	   <<  m  << ' '                 
	   << Nup << ' '                 
	   << LCC << ' '                 
	   << nup << ' '                 
	   << lcc << ' '                 

	   << endl;

      double dm= m-m0;
      if(dm>dm_max) { 
	dm_max = dm;
	Hc = H-dH;
	m0 = m;
      }
        
 
    } 
    fout.close();

    float t = timer() - t0;                
    cout << "\n The rough ramp of the M-H curve takes " << t << " s \n";


    char cmd[256];
    
    sprintf(cmd,"./showps  ./  %s &",parameters);
    system(cmd);
	
    
    cout << "\n The largest avlanche occurs at H=Hc= " << Hc << " with size= " << dm_max*0.5 << " \n";
    C0.SetEffectiveField(Hc);   
    Calculate_GS_M(C0, false, NONE); 
    
    C0.CalE0();                     
    C0.CalE(Hc);                     
    double E0 = C0.GetE0();
    double E  = C0.GetE();          
    double m  = C0.Getm();          
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

  }


  
 
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
    
    C0.DFS();
 
    char fname1[256]; 
    
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
      C.SetEffectiveField(H);                
      Calculate_GS_M(C, false, NONE);        
      double m = C.Getm();          
      int Nup = C.Get_Nup();
      double nup = C.Get_nup();

      cout << Normal;
      cout << "H= "                    << H 
	   << "; m= "                  << m 
	   << "; Nup= "                << Nup 
	   << "; nup= "                << nup 
	   << endl; 

      fout <<  H  << ' '                 
	   <<  m  << ' '                 
	   << Nup << ' '                 
	   << nup << ' '                 
	 << endl;
 
      double dm= m-m0;
      if(dm>dm_max) { 
	dm_max = dm;
	Hc = H-dH;
	m0 = m;
      }
        
 
    } 
    fout.close();

    float t = timer() - t0;                
    cout << "\n The rough ramp of the M-H curve takes " << t << " s \n";
   
    cout << "\n The largest avlanche occurs at H=Hc= " << Hc << " with size= " << dm_max*0.5 << " \n";
    C0.SetEffectiveField(Hc);   
    Calculate_GS_M(C0, false, NONE); 
    
    C0.CalE0();                     
    C0.CalE(Hc);                     
    double E0 = C0.GetE0();
    double E  = C0.GetE();          
    double m  = C0.Getm();          
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


  }


  

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
      C.SetEffectiveField_typeII(lambda);                
      Calculate_GS_M(C, false, NONE);        
      double m = C.Getm();          
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


      fout <<  lambda << ' '                 
	   <<  m  << ' '                 
	   << Nup << ' '                 
	   << LCC << ' '                 
	   << nup << ' '                 
	   << lcc << ' '                 
	   << Module.score  << ' '       
	   << Module.size << ' ';        
      for(unsigned int a=0; a<P.size(); a++)           
	fout << P[a] << ' ';
      for(Nbl_itr p = Module.nodelist.begin(); p!= Module.nodelist.end(); p++) 
	fout << namelist[*p] << ' ';
      fout << endl;
     
 
    } 
    fout.close();

    float t = timer() - t0;                
    cout << "\n The rough ramp of the M-lambda curve takes " << t << " s \n";

    
    char fname2[256]; 
    sprintf(fname2,"./data/GeneSets-%s-%s", parameters,argv[6]);
    C0.WriteGeneSetsFile(fname2);


    char cmd[256];
    sprintf(cmd,"./showps  ./  %s &",parameters);
    system(cmd);



    sprintf(cmd,"./showps.enrichment ./data %s %s %d &",parameters,argv[6],Nbgs);
    system(cmd);
	
  }




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
	
      C.CalE0();                     
      C.CalE(H);                     
      double E0 = C.GetE0();
      double E  = C.GetE();          
      double m  = C.Getm();          

   
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

    }



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
  }




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

  }





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

  }


  cout << "\n Done!  &:) \n";
  cout << "--------------------------------------------------------------------------------------------------------\n";
  exit(0);
}






