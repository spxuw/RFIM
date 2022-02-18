#include "State.h"

using namespace std;



void State::CalMetastableState(double H0)
{
  
  M = -N;
  for(int i=0; i<N; i++) {
    spin[i] = DOWN;
    nUp[i] = 0;
  }

  int curLoc = FindNext(UP); 
  
  
  H = H0; 
  

  
  if(curLoc<0) return;
  

  
  spinFlipQueue.push(curLoc);
  

  
  spinFlipQueue.push(-1);

  
  time = 0;

  
  while(!spinFlipQueue.empty()) {
    
    curLoc = spinFlipQueue.front();
    spinFlipQueue.pop();
		
    
    
    
    if(curLoc==-1) {
      time++;
      if(!spinFlipQueue.empty())
	spinFlipQueue.push(-1);
      continue;
    }

    
    
    if(spin[curLoc]!=UP)   
      FlipSpin(curLoc, UP);
    else 
      continue;

    
    
    for(Nbl::iterator p = A[curLoc].begin(); p!= A[curLoc].end(); p++) {
      int nextLoc = (*p);
      nUp[nextLoc]++;
      if( spin[nextLoc]!=UP && ((2.0*nUp[nextLoc]-K[nextLoc])+h[nextLoc]+H)>0.0 ) {
	spinFlipQueue.push(nextLoc);
      }
    }
  }

  AvalancheCount++; 

}




void State::HysteresisLoop(char* fname)
{
  char filename [256];
  sprintf(filename,"%s", fname);
  ofstream MofH_out(filename, ios_base::out);
  cout << "Outputting MofH to " << filename << endl;

  sprintf(filename,"%s.DS", fname);
  ofstream Size_out(filename, ios_base::out);
  cout << "Outputting Avalanche Size to " << filename << endl;


  
  double M_last = -N;
  for(int i=0; i<N; i++) {
    spin[i] = DOWN;
    nUp[i] = 0;
  }
  M = -N;
  m = -1.0;
  AvalancheCount = 0 ;
  FindNext(UP); 



  cout << ": Major Loop,  H ramps Up ......" << endl; 
  while ( M < N)  {   
    double m = Getm(); 
    
    module Module;
    CalModule(Module);

    
    cout << Normal << H << ' '  << m << ' ' << Blue << Module.score  << ' ' << Red << Module.size << " : " << Green;
    MofH_out <<           H << ' '  << m << ' ' <<         Module.score  << ' '        << Module.size << " : ";
    for(Nbl_itr p = Module.nodelist.begin(); p!= Module.nodelist.end(); p++) {
      cout << ' ' << NodeName[*p];
      MofH_out << ' ' << NodeName[*p];
    }
    cout << Normal << endl;
    MofH_out << endl;
    

    FlipNext(UP); 
    
    Size_out << H << ' ' << fabs(M-M_last)/2 << endl;;
    M_last = M;
    MofH_out << H  << ' ' << m << endl; 
  }
  MofH_out << H << ' ' << Getm() << endl;

    
  cout << ": Major Loop,  H ramps Down ......" << endl; 
  while(GetM() > -N)  {
    double m = Getm(); 
    MofH_out << H << ' ' << m << endl;
    
    FlipNext(DOWN); 
    
    Size_out << H << ' ' << fabs(M-M_last)/2 << endl;;
    M_last = M;
    MofH_out << H << ' ' << m << endl; 
  }
  MofH_out << H << ' ' << Getm() << endl;


  cout << "Done! There are totally " << GetAvalancheCount() << " avalanches.\n";
}







void State::FlipNext(char direction)
{
  
  
  if(direction==0)
    direction = -1;   
 

  
  int curLoc = FindNext(direction);
	
  
  if(curLoc<0) return;
  

  
  spinFlipQueue.push(curLoc);
  

  
  spinFlipQueue.push(-1);

  
  time = 0;

  
  while(!spinFlipQueue.empty()) {
    
    curLoc = spinFlipQueue.front();
    spinFlipQueue.pop();
		
    
    
    
    if(curLoc==-1) {
      time++;
      if(!spinFlipQueue.empty())
	spinFlipQueue.push(-1);
      continue;
    }

    
    
    if(spin[curLoc]!=direction)   
      FlipSpin(curLoc, direction);
    else 
      continue;

    
    
    for(Nbl::iterator p = A[curLoc].begin(); p!= A[curLoc].end(); p++) {
      int nextLoc = (*p);
      nUp[nextLoc] += direction;
      if( spin[nextLoc]!=direction && 
	  ((2.0*nUp[nextLoc]-K[nextLoc])+h[nextLoc]+H)*direction>0.0 ) {
	spinFlipQueue.push(nextLoc);
	
      }
    }
  }

  AvalancheCount++;  
}










int State::FindNext(char direction)
{
  
  
  double maxLocalFieldMinusH = -1e300*direction;   

  double localFieldMinusH;
  int nextLoc = -1;

  
  
  for(int i=0;i<N;i++) {
    if(spin[i]!=direction) {
      localFieldMinusH = h[i] + (2.0*nUp[i]-K[i]);
      if(localFieldMinusH*direction>=maxLocalFieldMinusH*direction) {
	maxLocalFieldMinusH=localFieldMinusH;
	nextLoc = i;
      }
    }
  }
  
  if(nextLoc>=0) H = -maxLocalFieldMinusH;
  return nextLoc;
}











