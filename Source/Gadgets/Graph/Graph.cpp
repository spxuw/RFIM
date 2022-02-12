#include "Graph.h"

Graph::Graph(){}
Graph::~Graph(){}

///////////////////////////////////////////////////////////////////////////////
Graph::Graph(MATRIX& LinkWeight, double t)
{
  //cout << "\nConstruct the graph at threshold t= " << t << endl;
  N = LinkWeight.size();
  E = 0;
  A.resize(N);
  K.resize(N);
  
  AM.resize(N);
  for(int i=0; i<N; i++) 
    AM[i].resize(N,0);
  
  for(int i=0; i<N; i++) {
    for(int j=i+1; j<N; j++) {
      if(LinkWeight[i][j]>=t) {
	AddLink(i,j);
	//cout << LinkWeight[i][j] << ' ' << t << endl; //test
      }
    }
  }
}
///////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
// construct Erdos-Renyi random network from given N, p, and seed
// actually this is the so-called Gilbert model, or Gnp model
void Graph::Erdos_Renyi(int n, double p, int seed)
{
  N = n; 
  K.resize(N);
  A.resize(N);
  rand.seed(seed);	

  AM.resize(N);
  for(int i=0; i<N; i++) 
    AM[i].resize(N,0);

  for (int i=0; i<N; i++) 
    K[i] = 0;

  // link nodes randomly according to connection probability p
  E = 0;
  for (int i=0; i<N; i++) {
    for (int j=i+1; j<N; j++) {
      if(rand.ran1() < p) {
	AddLink(i,j);
      }
    }
  }
  cout << "E= " << E << endl; //test
  cout << "L.size() = " << L.size() << endl;  //test

}
///////////////////////////////////////////////////////////////////////////////



///////////////////////////////////////////////////////////////////////////////
// Make a copy of G
void Graph::Backup()
{
  Ncopy = N;
  Ecopy = E;
  
  Kcopy.clear(); Kcopy.resize(N);
  Acopy.clear(); Acopy.resize(N);
  copy(K.begin(), K.end(), Kcopy.begin());
  copy(A.begin(), A.end(), Acopy.begin());
}
///////////////////////////////////////////////////////////////////////////////


///////////////////////////////////////////////////////////////////////////////
// Make a copy of G
void Graph::Restore()
{
  N = Ncopy;

  K.clear(); K.resize(N);
  A.clear(); A.resize(N);
  L.clear();
  E = 0;
  for(int i=0; i<N; i++) {
    for(Nbl::iterator p=Acopy[i].begin(); p!=Acopy[i].end(); p++) {
      int j=*p;
      if(!find(A[i],j)) 
	AddLink(i,j);
    }
  }
    
}
///////////////////////////////////////////////////////////////////////////////


///////////////////////////////////////////////////////////////////////////////
// AWARE: don't double count links
void Graph::AddLink(int i, int j)
{
  Link link(i, j, E+1); 
  L.push_back(link);

  // build a mapping between the edge (i,j) and the edge index e
  stringstream sst; 
  sst << i << ">" << j;  
  MAP[sst.str()] = E ;
  //cout << E << ',' << link << ',' << MAP[sst.str()] << endl; //debug
  E++;
  
 
  A[i].push_back(j); 
  A[j].push_back(i); 
  
  K[i]++;
  K[j]++;

  AM[i][j]=AM[j][i]=1;
}
///////////////////////////////////////////////////////////////////////////////




///////////////////////////////////////////////////////////////////////////////
void Graph::RemoveLink(int i, int j)
{
  if(find(A[i],j)) {
    E--;
  
    A[i].remove(j); 
    A[j].remove(i); 
    
    K[i]--;
    K[j]--;

    AM[i][j]=AM[j][i]=0;
  }
}
///////////////////////////////////////////////////////////////////////////////


///////////////////////////////////////////////////////////////////////////////
void Graph::AddBackLink(int i, int j)
{
  if(!find(A[i],j)) {
    E++;
  
    A[i].push_back(j); 
    A[j].push_back(i); 
  
    K[i]++;
    K[j]++;

    AM[i][j]=AM[j][i]=1;
  }
}
///////////////////////////////////////////////////////////////////////////////


///////////////////////////////////////////////////////////////////////////////
void Graph::RemoveNode(int i)
{

  list<int> Ai;
  //copy(A[i].begin(), A[i].end(), Ai.begin()); // this doesnot work!!!
  for(Nbl_itr p = A[i].begin(); p!= A[i].end(); p++) 
    Ai.push_back(*p);

  for(Nbl_itr p = Ai.begin(); p!= Ai.end(); p++) {
    int j = (*p);
    E--;
  
    A[i].remove(j); 
    A[j].remove(i); 
    
    K[i]--;
    K[j]--;

    AM[i][j] = 0;
    AM[j][i] = 0;
  }
}
///////////////////////////////////////////////////////////////////////////////



/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Remove a set of nodes
void Graph::Remove(set<int>& V) 
{
  for(set<int>::iterator p=V.begin(); p!=V.end(); p++) {
    RemoveNode(*p);
    //cout << "K[" << *p << "]=" << K[*p] << "; ";
  }

  //cout << K[12] << ',' << *(A[12].begin()) << ',' << K[14] << endl; 

}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////





/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
bool Graph::Plot(char* fname)  
{
    FILE *fp;
    int i, j;
    char filename[255];
    char cmd[255];

    sprintf(filename, "%s.dot", fname);
    fp = fopen(filename, "wb");

    fputs("graph vis {\n", fp);

    //fputs("\tgraph [bgcolor=bisque, size=\"5,5\", ratio=fill, center=1];\n", fp);
    fputs("\tgraph [size=\"5,5\", ratio=fill, center=1];\n", fp);
    //fputs("\tgraph [size=\"15,15\", ratio=fill, center=1];\n", fp);
    
    if(N>50)
      fputs("\tnode [label=\"\", shape=point,style=\"filled\"];\n", fp);
    //fputs("\tnode [label=\"\",shape=circle,style=\"filled\"];\n", fp);
    else
      fputs("\tnode [shape=circle,style=\"filled\"];\n", fp);
    //fputs("\tedge [color=blue];\n", fp);
    fputs("\tedge [color=black];\n", fp);



    for(i=0; i<N; i++) {
      fprintf(fp, "\t%d [color=\"gray\",label=\"%d\"];\n", i, i);
    }
     
    for(int l=0; l<E; l++) {
      i = L[l].GetHead(); // starting node of the l-th link 
      j = L[l].GetTail(); // ending node of the l-th link 
      if(K[i]>0 && K[j]>0)
	fprintf(fp, "\t%d -- %d [color=blue];\n", i, j);
    }
  // end: loop of links
    
    fputs("}\n", fp);
    fclose(fp);


    //'dot' is designed to plot directed graph (digraph).
    //undirected graph or just called graph can be plotted with 'neato' or 'fdp'
    //other radio graph can be plotted with 'twopi'
    //circle graph can be plotted with 'circo'
    //All the commands share the same command line arguments.

    
    sprintf(cmd, "neato -Tps %s.dot > %s.neato.eps", fname, fname);
    system(cmd);
    
    /*
    sprintf(cmd, "twopi -Tps %s.dot > %s.twopi.ps", fname, fname);
    system(cmd);
    
    sprintf(cmd, "fdp -Tps %s.dot > %s.fdp.ps", fname, fname);
    system(cmd);
    */

    //sprintf(cmd, "circo -Tps %s.dot > %s.circo.ps", fname, fname);
    //system(cmd);
    
  
    
    sprintf(cmd, "gv %s.neato.eps & ", fname);
    system(cmd);

    /*
    sprintf(cmd, "gv %s.fdp.ps & ", fname);
    system(cmd);

    sprintf(cmd, "gv %s.twopi.ps & ", fname);
    system(cmd);
    */

    //sprintf(cmd, "gv %s.circo.ps & ", fname);
    //system(cmd);
    

}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// the following DFS is designed to find the largest connected component in a directed or undirected graph
// for directed graph, the DFS is performed by searching the A_list, rather than Aout_list, of each node.
// So it cannot give the topological order of nodes.
void Graph::DFS()
{
  vector<bool> visited(N, false);
  AllComponents.clear();
  
  int maxsize  = 0; // size of the largest connected component
  int lccindex = 0; // index of the largest connected component
  int ncc  = -1;    // count the number of connected component
  for(int i=0; i<N; i++) {
    if(!visited[i]) {
      ncc++; // find another connected component
      Nbl Component;
      DFS(i, visited, Component);
      AllComponents.push_back(Component);

      int size = Component.size();
      if(size>maxsize) {
	maxsize = size;
	lccindex = ncc;
      }
    }
  }
  //cout << "\nThe largest connected component (label " << lccindex << ") has size = " << maxsize; //test

  Nlc = maxsize; // return the size of the largest connected component
  LC.clear();
  LC.resize(N);           // this vector labeling the node whether it belongs to the largest connected component

  Ncc = AllComponents.size();
  //cout << "There are total " << Ncc << " connected components.\n";

  ccindex.clear();
  ccindex.resize(N);
  int cc = 0;	
  Nac = 0;
  for(vector<Nbl>::iterator CC = AllComponents.begin(); CC!=AllComponents.end(); CC++) {
    for(Nbl_itr p = (*CC).begin(); p!= (*CC).end(); p++) {
      ccindex[*p]= cc;
    }
    int size = (*CC).size();
    Nac += size;
    cc++;
  }
  Nac = (Nac-Nlc)/(Ncc-1);

  cc = 0;	
  int sumKlc = 0;    // sum of k for the largest connected component
  for(vector<Nbl>::iterator CC = AllComponents.begin(); CC!=AllComponents.end(); CC++) {
    if(cc==lccindex) { // if this is the largest component
	
      int count=0;
      for(Nbl_itr p = (*CC).begin(); p!= (*CC).end(); p++, count++) {
	//cout << (*p) << ' ' << K[(*p)] << endl; //debug
	//cout << (*p) << ' '; //debug

	LClist.push_back(*p);    // push it into the LClist
	LC[(*p)] = 1;            // label it 
	sumKlc += K[(*p)];
      }
      break;
    }
    else
      cc++;
  }
    
  // calculate the number of edges in the largest connected component
  Elc = sumKlc/2;
  mlc = Elc/(double)Nlc;
  //cout << "\n Nlc = " << Nlc << " sumKlc= " << sumKlc << "  Elc=sumKlc/2= " << Elc << " mlc = Elc/Nlc= " << mlc << endl;

  //cout << "\n Nlc = " << Nlc << "  Elc= " << Elc << endl;


}


 
void Graph::DFS(int u, vector<bool>& visited, Nbl&  Component)
{
  visited[u] = true;
  Component.push_back(u);
    
  //cout << "u = " << u ; //test
  // loop over node u's adjacent nodes
  for(Nbl_itr p = A[u].begin(); p!= A[u].end(); p++) {
    //for(Nbl_itr p = pseudoA[u].begin(); p!= pseudoA[u].end(); p++) { // test for directed graph
    int v = (*p);
    //cout << "v= " << v << endl ; //test
    //printf(" v= %d \n", v); //test
    if(!visited[v])
      DFS(v, visited, Component);
  }
}
////////////////////////////////////////////////////////////////////////////////////////////////



////////////////////////////////////////////////////////////////////////////////////////////////
// for a given source node s, find its distance to the given target
int Graph::BFS(int source, int target)
{
    stack<int>  S;             // return nodes in order of non-increasing distance from the source
    vector<int> d(N,-1);

    //initialize the source vertex
    d[source] = 0;
 
    queue<int> Q; // a FIFO queue   
    Q.push(source); //insert source into the queue
    while(!Q.empty())  {
      int u = Q.front();   // put the value of the queue's head to u 
      Q.pop();             // delete the head of the queue 
      S.push(u);
      // loop over node u's adjacent nodes
      for(Nbl_itr p = A[u].begin(); p!= A[u].end(); p++) {
        int v = (*p);

        // if this node has never been visited before.
        if(d[v] < 0) { 
          d[v]     = d[u]+1;    // get its distance to the source
          Q.push(v);            // put it into the queue

          if(v==target) {
            return d[v];
          }
        }
      }// end of loop over node u's adjacent nodes
    }// end of while loop

    if(d[target]==-1) 
      d[target] = 2*N;

    return d[target]; // if d[target] is still -1, this just means those two nodes are not connected 
}
////////////////////////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////////////////////////
// for a given source node s, find its distances to all other nodes in the graph
void Graph::BFS(int source, VECTOR& d)
{
    stack<int>  S;             // return nodes in order of non-increasing distance from the source
    
    //initialize
    for(int i=0; i<N; i++) {
      d[i] = -1;
    }

    //initialize the source vertex
    d[source] = 0;
 
    queue<int> Q; // a FIFO queue   
    Q.push(source); //insert source into the queue
    while(!Q.empty())  {
      int u = Q.front();   // put the value of the queue's head to u 
      Q.pop();             // delete the head of the queue 
      S.push(u);
      // loop over node u's adjacent nodes
      for(Nbl_itr p = A[u].begin(); p!= A[u].end(); p++) {
        int v = (*p);
        // if this node has never been visited before.
        if(d[v] < 0) { 
          d[v]     = d[u]+1;    // get its distance to the source
          Q.push(v);            // put it into the queue
        }
      }// end of loop over node u's adjacent nodes
    }// end of while loop


    while(!S.empty())  {
      int v = S.top();  
      S.pop();          
      if (d[v]>diameter) 
        diameter = d[v];
    }// end of while loop

    for(int i=0; i<N; i++) {
      if(d[i]==-1)
        d[i] = 2*N; // this just means those two nodes are not connected 
    }
}
////////////////////////////////////////////////////////////////////////////////////////////////



////////////////////////////////////////////////////////////////////////////////////////////////
// calculate the neighborhood matrix, each matrix element Mij indicates the value of the shortest path between nodes i and j.  See Eq.(3) of R.F.S. Andrade et al. / Physics Letters A 372 (2008) 5265–5269
void Graph::CalNeighborhoodMatrix()
{
  M.clear();
  M.resize(N);

  diameter = -1;
  for(int i=0; i<N; i++) {
    M[i].resize(N);
    BFS(i, M[i]);
    
    for(int j=0; j<N; j++) {
      if(M[i][j]==2*N) // so we assume if node i and j are not reachable, then their distance is 0. Is this good?
        M[i][j] = 0; 
    }
  }

  //cout << "diameter = " <<  diameter << endl;
  // Note that for a network with zero edges, the diameter would be 0.
}
////////////////////////////////////////////////////////////////////////////////////////////////



////////////////////////////////////////////////////////////////////////////////////////////////
// calculate the distance between two neighborhood matrices. See Eq.(4) of R.F.S. Andrade et al. / Physics Letters A 372 (2008) 5265–5269
// Note that the neighborhood matrices have to be calculated already. 
double Graph::GetNetworkDistance(Graph& G0)
{
  double d = 0;
  for(int i=0; i<N; i++) {
    for(int j=i+1; j<N; j++) {
      double t = M[i][j]/(double)diameter - G0.M[i][j]/(double)G0.diameter;
      d += t*t;
    }
  }
  d *= 2.0/(double)(N*(N-1));
  return d;
}
////////////////////////////////////////////////////////////////////////////////////////////////



////////////////////////////////////////////////////////////////////////////////////////////////
// calculate the distance between two adjacency matrices. 
double Graph::GetMatrixDistance(Graph& G0)
{
  double d = 0;
  for(int i=0; i<N; i++) {
    for(int j=i+1; j<N; j++) {
      double t = AM[i][j]-G0.AM[i][j];
      d += t*t;
    }
  }
  d *= 2.0/(double)(N*(N-1));
  return d;
}
////////////////////////////////////////////////////////////////////////////////////////////////





