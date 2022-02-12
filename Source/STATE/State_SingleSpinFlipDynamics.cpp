#include "State.h"

using namespace std;

/////////////////////////////////////////////////////////////////////////////////////////////////////////
// Calculate the metastable state at external field H=H0
void State::CalMetastableState(double H0)
{
  // initialization 
  M = -N;
  for(int i=0; i<N; i++) {
    spin[i] = DOWN;
    nUp[i] = 0;
  }

  int curLoc = FindNext(UP); // this will set H = -maxLocalFieldMinusH
  
  ////////////////////////////////////////////////////////////
  H = H0; // we need to reset the H value
  ////////////////////////////////////////////////////////////

  // If there were no spins left to flip, return
  if(curLoc<0) return;
  //else cout << "curLoc= "<< curLoc << " H= " << H << "   m= " << M/(double)N << endl; 

  // Push the first spin onto the queue and the Loc List
  spinFlipQueue.push(curLoc);
  //spinLocList.push_back(curLoc);     

  // Push a marker indicating the end of the shell onto the queue
  spinFlipQueue.push(-1);

  // Set the current time to zero
  time = 0;

  // As long as there are spins left to flip, propagate the avalanche
  while(!spinFlipQueue.empty()) {
    // retrieve the next spin from the queue
    curLoc = spinFlipQueue.front();
    spinFlipQueue.pop();
		
    // if we got the marker, then the shell is over, so we increment
    // the time, and put a new marker back on the end of the queue
    // to indicate the end of the next shell
    if(curLoc==-1) {
      time++;
      if(!spinFlipQueue.empty())
	spinFlipQueue.push(-1);
      continue;
    }

    // if we haven't flipped the spin yet, flip it.  If we have flipped
    // it, continue with the next spin.
    if(spin[curLoc]!=UP)   
      FlipSpin(curLoc, UP);
    else 
      continue;

    // loop over the spin's neighbors and put them on the queue if they
    // are ready to flip
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



/////////////////////////////////////////////////////////////////////////////////////////////////////////
void State::HysteresisLoop(char* fname)
{
  char filename [256];
  sprintf(filename,"%s", fname);
  ofstream MofH_out(filename, ios_base::out);
  cout << "Outputting MofH to " << filename << endl;

  sprintf(filename,"%s.DS", fname);
  ofstream Size_out(filename, ios_base::out);
  cout << "Outputting Avalanche Size to " << filename << endl;


  // initialization 
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
  while ( M < N)  {   // While there are spins left to flip
    double m = Getm(); 
    
    module Module;
    CalModule(Module);

    //MofH_out << H << "     " << m << endl; 
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
/////////////////////////////////////////////////////////////////////////////////////////////////////




/////////////////////////////////////////////////////////////////////////////////////////////////////////
// FlipNext: Flip the next avalanche in the specified direction.
void State::FlipNext(char direction)
{
  // Note that in this code, direction should take values +1 (UP) or -1 (DOWN)
  // But in State.h, we define UP=+1, and DOWN=0.  
  if(direction==0)
    direction = -1;   
 

  // Find the first spin of the avalanche
  int curLoc = FindNext(direction);
	
  // If there were no spins left to flip, return
  if(curLoc<0) return;
  //else cout << "curLoc= "<< curLoc << " H= " << H << "   m= " << M/(double)N << endl; 

  // Push the first spin onto the queue and the Loc List
  spinFlipQueue.push(curLoc);
  //spinLocList.push_back(curLoc);     

  // Push a marker indicating the end of the shell onto the queue
  spinFlipQueue.push(-1);

  // Set the current time to zero
  time = 0;

  // As long as there are spins left to flip, propagate the avalanche
  while(!spinFlipQueue.empty()) {
    // retrieve the next spin from the queue
    curLoc = spinFlipQueue.front();
    spinFlipQueue.pop();
		
    // if we got the marker, then the shell is over, so we increment
    // the time, and put a new marker back on the end of the queue
    // to indicate the end of the next shell
    if(curLoc==-1) {
      time++;
      if(!spinFlipQueue.empty())
	spinFlipQueue.push(-1);
      continue;
    }

    // if we haven't flipped the spin yet, flip it.  If we have flipped
    // it, continue with the next spin.
    if(spin[curLoc]!=direction)   
      FlipSpin(curLoc, direction);
    else 
      continue;

    // loop over the spin's neighbors and put them on the queue if they
    // are ready to flip
    for(Nbl::iterator p = A[curLoc].begin(); p!= A[curLoc].end(); p++) {
      int nextLoc = (*p);
      nUp[nextLoc] += direction;
      if( spin[nextLoc]!=direction && 
	  ((2.0*nUp[nextLoc]-K[nextLoc])+h[nextLoc]+H)*direction>0.0 ) {
	spinFlipQueue.push(nextLoc);
	//spinLocList.push_back(nextLoc); 
      }
    }
  }

  AvalancheCount++;  
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////




/////////////////////////////////////////////////////////////////////////////////////////////////////////
// Find the first spin of the next avalanche
// Here, I am using the BruteForce algorithm.
// I am thinking that due to the fact that inhomogeneous degree (# of nearest neighbors) of each spin
// the SortedList algorithm cannot be applied. Is this true??
int State::FindNext(char direction)
{
  // initialize maxLocalFieldMinusH to a large negative number so that it will
  // be less than any of the actual local fields.
  double maxLocalFieldMinusH = -1e300*direction;   // negative enough

  double localFieldMinusH;
  int nextLoc = -1;

  // loop over all the spins in the lattice, and find the unflipped spin with
  // the largest local field.
  for(int i=0;i<N;i++) {
    if(spin[i]!=direction) {
      localFieldMinusH = h[i] + (2.0*nUp[i]-K[i]);
      if(localFieldMinusH*direction>=maxLocalFieldMinusH*direction) {
	maxLocalFieldMinusH=localFieldMinusH;
	nextLoc = i;
      }
    }
  }
  // change H to flip the spin, and return the location of the spin
  if(nextLoc>=0) H = -maxLocalFieldMinusH;
  return nextLoc;
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////










