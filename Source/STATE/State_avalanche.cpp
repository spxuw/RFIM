#include "State.h"

using namespace std;


void Cal_Avalanche(State &A, State &B, list<int> &avalanche)
{
    avalanche.clear();

    int N = A.GetN();
    for (int i=0; i<N; i++)
    {
	if(A.GetSpin(i) != B.GetSpin(i))
	    avalanche.push_back(i);
    }
}






void Show_Avalanche(State& A, State& B)
{
    list<int> avalanche;
    int N = A.GetN();
    vector<char> lattice(N);

    for (int i=0; i<N; i++)
    {
	if(A.GetSpin(i) != B.GetSpin(i))
	{
	    avalanche.push_back(i);
	    lattice[i] = 1;
	}
	else
	    lattice[i] = 0;
    }

    cout << endl;

    int L = A.GetL();
    for (int j = 0; j < L; ++j) 
    {
	for (int i = 0; i < L; ++i) 
	{
	    if (lattice[i + L * j] == 1) 
		cout << '1';
	    else
		cout << '0';
	}
	cout << "\n";
    }
    
}






void State::Show_State()
{
    cout << "\n\n";
    if(D>=2)
	for (int j = 0; j < L; ++j) 
	{
	    for (int i = 0; i < L; ++i) 
	    {
		if (spin[i + L * j] == UP) 
		    cout << '1';
		else
		    cout << '0';
	    }
	    cout << "\n";
	}
    else
    {
	for (int i = 0; i < L; ++i) 
	{
	    if (spin[i] == UP) 
		cout << '1';
	    else
		cout << '0';
	}
	cout << "\n";
    }

    cout << "\n\n";

}







void State::Save_State(int k)
{
    char fname[256]; 
    if(!network)
      sprintf(fname,"./data/GS-D%d-L%d-R%.3lf-seed%d/%d",D,L,R,seed,k);
    else 
      sprintf(fname,"./data/GS-N%d-c%e-R%.3lf-seed%d/%d",N,c,R,seed,k);


    ofstream outfile(fname, ios::out|ios::binary);
    if (!outfile)  
    {
	cout << "cannot open file:  " << fname << "for writing\n";
	exit(2);
    }
    
    
    outfile.write(reinterpret_cast<char *>(&M),sizeof(int));

    
    
    outfile.write(spin, N*sizeof(char));
	
    
    outfile.close();
}







void State::Save_State(string file, int k)
{
    char fname[256]; 
    sprintf(fname,"./data/GS-%s/%d",file.c_str(),k);

    ofstream outfile(fname, ios::out|ios::binary);
    if (!outfile)  
    {
	cout << "cannot open file:  " << fname << "for writing\n";
	exit(2);
    }
    
    
    outfile.write(reinterpret_cast<char *>(&M),sizeof(int));

    
    
    outfile.write(spin, N*sizeof(char));
	
    
    outfile.close();
}


	    

void State::getRGBcolor(int avalancheindex, double& r,double& g,double& b)
{
    
    

    int Color = (avalancheindex%216); 

    int remainder = Color%36;
    int R = Color/36;

    int G = remainder/6;
    remainder %= 6;

    int B = remainder;

    r = R/5.0;
    g = G/5.0;
    b = B/5.0;
}





void State::DrawOneAvalanche(int index, list<int> & avalanche)
{
    list<int>::iterator itr; 
    for(itr= avalanche.begin(); itr != avalanche.end(); itr++)
    {
	int curLoc  = (*itr); 
	avalancheindex[curLoc] = index;
    }
}





void State::DrawAllAvalanches(char* fname)
{
    char epsfile[256];
    sprintf(epsfile,"%s.eps",fname);

    ofstream fout(epsfile);

    fout << "%!PS-Adobe-3.0 EPSF-3.0\n";
    fout << "%%BoundingBox: -5 -5 " << L + 25 << ' ' << L + 25 << " \n";
    fout << "/square {newpath moveto 1 0 rlineto 0 1 rlineto -1 0 rlineto closepath fill} def\n";
    fout << "10 10 translate\n";

   
    int spincount = 0;
    double r,g,b;

    
    for (int j = L-1; j>=0; j--) {
	for (int i = 0; i < L; i++) {
	    
	    getRGBcolor(avalancheindex[spincount],r,g,b);
	    fout << i << " " << j << " square " << r << ' ' << g << ' ' << b << " setrgbcolor\n";
	    spincount++;
	}
    }
    
    fout.close();
}


