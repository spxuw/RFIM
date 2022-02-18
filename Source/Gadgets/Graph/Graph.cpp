#include "Graph.h"

Graph::Graph(){}
Graph::~Graph(){}

Graph::Graph(MATRIX& LinkWeight, double t)
{
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
      }
    }
  }
}

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

  E = 0;
  for (int i=0; i<N; i++) {
    for (int j=i+1; j<N; j++) {
      if(rand.ran1() < p) {
	AddLink(i,j);
      }
    }
  }
  cout << "E= " << E << endl; 
  cout << "L.size() = " << L.size() << endl;  

}


void Graph::Backup()
{
  Ncopy = N;
  Ecopy = E;
  
  Kcopy.clear(); Kcopy.resize(N);
  Acopy.clear(); Acopy.resize(N);
  copy(K.begin(), K.end(), Kcopy.begin());
  copy(A.begin(), A.end(), Acopy.begin());
}

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

void Graph::AddLink(int i, int j)
{
  Link link(i, j, E+1); 
  L.push_back(link);

  stringstream sst; 
  sst << i << ">" << j;  
  MAP[sst.str()] = E ;
  E++;
  
 
  A[i].push_back(j); 
  A[j].push_back(i); 
  
  K[i]++;
  K[j]++;

  AM[i][j]=AM[j][i]=1;
}



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


void Graph::RemoveNode(int i)
{

  list<int> Ai;
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



void Graph::Remove(set<int>& V) 
{
  for(set<int>::iterator p=V.begin(); p!=V.end(); p++) {
    RemoveNode(*p);
  }


}





bool Graph::Plot(char* fname)  
{
    FILE *fp;
    int i, j;
    char filename[255];
    char cmd[255];

    sprintf(filename, "%s.dot", fname);
    fp = fopen(filename, "wb");

    fputs("graph vis {\n", fp);

    fputs("\tgraph [size=\"5,5\", ratio=fill, center=1];\n", fp);
    
    if(N>50)
      fputs("\tnode [label=\"\", shape=point,style=\"filled\"];\n", fp);
    else
      fputs("\tnode [shape=circle,style=\"filled\"];\n", fp);
    fputs("\tedge [color=black];\n", fp);



    for(i=0; i<N; i++) {
      fprintf(fp, "\t%d [color=\"gray\",label=\"%d\"];\n", i, i);
    }
     
    for(int l=0; l<E; l++) {
      i = L[l].GetHead();
      j = L[l].GetTail();
      if(K[i]>0 && K[j]>0)
	fprintf(fp, "\t%d -- %d [color=blue];\n", i, j);
    }
    
    fputs("}\n", fp);
    fclose(fp);

    
    sprintf(cmd, "neato -Tps %s.dot > %s.neato.eps", fname, fname);
    system(cmd);
  
    
    sprintf(cmd, "gv %s.neato.eps & ", fname);
    system(cmd);

    

}


void Graph::DFS()
{
  vector<bool> visited(N, false);
  AllComponents.clear();
  
  int maxsize  = 0; 
  int lccindex = 0;
  int ncc  = -1; 
  for(int i=0; i<N; i++) {
    if(!visited[i]) {
      ncc++; 
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

  Nlc = maxsize; 
  LC.clear();
  LC.resize(N);           
  Ncc = AllComponents.size();

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
  int sumKlc = 0;   
  for(vector<Nbl>::iterator CC = AllComponents.begin(); CC!=AllComponents.end(); CC++) {
    if(cc==lccindex) { 
	
      int count=0;
      for(Nbl_itr p = (*CC).begin(); p!= (*CC).end(); p++, count++) {

	LClist.push_back(*p);   
	LC[(*p)] = 1;            
	sumKlc += K[(*p)];
      }
      break;
    }
    else
      cc++;
  }
    
  Elc = sumKlc/2;
  mlc = Elc/(double)Nlc;

}


 
void Graph::DFS(int u, vector<bool>& visited, Nbl&  Component)
{
  visited[u] = true;
  Component.push_back(u);
    
  for(Nbl_itr p = A[u].begin(); p!= A[u].end(); p++) {
    int v = (*p);
    if(!visited[v])
      DFS(v, visited, Component);
  }
}

int Graph::BFS(int source, int target)
{
    stack<int>  S;             
    vector<int> d(N,-1);

    d[source] = 0;
 
    queue<int> Q;   
    Q.push(source); 
    while(!Q.empty())  {
      int u = Q.front();   
      Q.pop();             
      S.push(u);
      
      for(Nbl_itr p = A[u].begin(); p!= A[u].end(); p++) {
        int v = (*p);

        
        if(d[v] < 0) { 
          d[v]     = d[u]+1;    
          Q.push(v);            

          if(v==target) {
            return d[v];
          }
        }
      }
    }

    if(d[target]==-1) 
      d[target] = 2*N;

    return d[target]; 
}





void Graph::BFS(int source, VECTOR& d)
{
    stack<int>  S;             
    
    
    for(int i=0; i<N; i++) {
      d[i] = -1;
    }

    
    d[source] = 0;
 
    queue<int> Q; 
    Q.push(source); 
    while(!Q.empty())  {
      int u = Q.front();   
      Q.pop();             
      S.push(u);
      
      for(Nbl_itr p = A[u].begin(); p!= A[u].end(); p++) {
        int v = (*p);
        
        if(d[v] < 0) { 
          d[v]     = d[u]+1;    
          Q.push(v);            
        }
      }
    }


    while(!S.empty())  {
      int v = S.top();  
      S.pop();          
      if (d[v]>diameter) 
        diameter = d[v];
    }

    for(int i=0; i<N; i++) {
      if(d[i]==-1)
        d[i] = 2*N; 
    }
}






void Graph::CalNeighborhoodMatrix()
{
  M.clear();
  M.resize(N);

  diameter = -1;
  for(int i=0; i<N; i++) {
    M[i].resize(N);
    BFS(i, M[i]);
    
    for(int j=0; j<N; j++) {
      if(M[i][j]==2*N) 
        M[i][j] = 0; 
    }
  }

  
  
}







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






