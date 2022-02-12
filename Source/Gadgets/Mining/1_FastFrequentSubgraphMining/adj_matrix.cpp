/********************************************************************
 * adj_matrix.cpp
 * Implementation of the class AdjMatrix
 *
 ********************************************************************/
#pragma warning (disable:4786 4018 4267)
#include <iostream>
#include <assert.h>
#include <vector>
#include <string>
#include <fstream>
#include <algorithm>
using namespace std;

#include "adj_matrix.h"
#include "common.h"
#include "Ullman.h"
#if defined (VFLIB)
	#include "argraph.h"
	#include "argedit.h"
	#include "attribute.h"
	#include "match.h"
	#include "ull_sub_state.h"
	#include "vf2_sub_state.h"
	#include "vf2_state.h"
	#include "vf2_sub_state.h"
#endif

bool AdjMatrix::denseEnough(int d, bool final){
     int s = size();
     if( s == 1 || s== 2) return true;
     int e = esize(), i;
     bool flag = true;
	 
#ifdef DENSITY
     int l = lastbit();
     int ad = s-1-l;
     if( final) ad = 0;
     if( s*(s-1)/2 -e -ad >  d) flag = false;
#endif

	 
#ifdef MINIMAL_LABEL
	 int limit = MINIMAL_LABEL + EMPTY_EDGE; 
	 int row = s-1;
	 if( final ) row = s;
	 for(  i=0; i<row && flag; i++){
		GLTYPE c  = getLabel(row-1, i);  //check the second to the last row
		if( c > limit ) 	break;   //must be at least the minimal bound
		//cout << "c is " << c - EMPTY_EDGE<< " ";
	 }
	 //cout << "i is " << i << " " << s << " " << MINIMAL_LABEL << endl;
	 if( i == row) flag = false; //we could not find a label with sufficiently low value;
#endif

	/*
     if( final) {
		 cout << "s e f l d " << s << " " << e <<  " " <<
		     flag << " " << l << " " << d <<endl;
		this->print();
	 }
	 */
     return flag;
}

// fill in vector v with neighbors of node index
vector<GNSIZE> * AdjMatrix::getNeighbors(int index)
{
	return SAdjList[index];
}

void AdjMatrix::getNeighbors(int index, vector<GNSIZE> & v)
{

	copy(SAdjList[index]->begin(), SAdjList[index]->end(), v.begin());

}


bool AdjMatrix::isPath()
{
        int ends = 0;
	for (int i = 1; i < size(); i++) {
            int degree = 0;
            for (int j = 0; j < size(); j++) { 
                if ( getLabel(i,j) != EMPTY_EDGE) {
                     degree++;
                     if (degree > 3)
                         return false;
                }
            }
            if (degree == 2) ends++;
        }

        return ends == 2;
}

// fill in vector v with nodes with max label
void AdjMatrix::getMaxNode(vector<int> &v)
{
	GLTYPE max = getLabel(0,0);
	int i;
	for (i = 1; i < size(); i++) {
		if (max < getLabel(i,i))
			max = getLabel(i,i);
	}

	for (i = 0; i < size(); i++) {
		if (getLabel(i,i) == max)
			v.push_back(i);
	}
}

void AdjMatrix::getMaxNode2(vector<KEYTYPE> &v)
{
	GLTYPE max = getLabel(0,0);
	int i;
	for (i = 1; i < size(); i++) {
		if (max < getLabel(i,i))
			max = getLabel(i,i);
	}

	for (i = 0; i < size(); i++) {
		if (getLabel(i,i) == max)
			v.push_back(aMAKEKEY(aNULLKEY,i));
	}
}


// return true if the matrix is canonical form

int AdjMatrix::isCanonicalForm(){
	
	//int t1 = isCanonicalForm1();
	int t2 = isCanonicalForm2();

	return t2;
	/*if( t1 != t2 ){ 
		print();
		cout << " t1 t2 is " << t1 << "\t" << t2 << endl;
		error("mis matching in canonical form\n"); return 0;
	}
	else return t1;
	*/
	
} 

int AdjMatrix::isCanonicalForm2() 
{
	vector<KEYTYPE>   searchTree;
	vector<GNSIZE>    nodeTree;

	if( size() == 2) return 1;  // for symmetrix adj matrix
	vector<GNSIZE> ancesters(size(), 0);
	int i, j, k, m, n, u, v, w, uplimit1, uplimit2;

	// level 1(single node)
	getMaxNode2(searchTree);

	//cout << "search 2\n" << "search tree size" << searchTree.size() << endl;
	if (getLabel(0,0) != getLabel( aSECONDKEY(searchTree[0]), aSECONDKEY(searchTree[0])) )
		return -1;    //false -> change to -1

	//initialize the nodeTree
	for( m =0; m< searchTree.size(); m++){
		nodeTree.push_back(0);
	}

	//set up l
	int bi, bl;
	last2bitil(bi, bl);    //to set the second to last bit's index and level

	//level wise search
	uplimit1 = 0;
	int code = 1;
	for (i = 1; i < size(); i++) {
		
		// step 1. propse next level candidates to vector next
		uplimit2 = searchTree.size();
		vector<KEYTYPE> refValues;

		//populate the refValues;
		for(m=0; m<=i; m++){
			if( getLabel(i,m)!= EMPTY_EDGE )
				refValues.push_back(aMAKEKEY(i-m, getLabel(i,m)));
		}

		for (j = uplimit1; j < uplimit2; j++) {
			u = searchTree[j];
					
#ifdef DEBUG1
			cout << "level " << i << ", #" <<j<<": " << hex <<  u  << dec << endl;
#endif
			
			//get the ancesters of u
			ancesters[i-1] = aSECONDKEY(u);
			int u1 = u;
			for(m=1; m<i; m++){
				int fkey = aFIRSTKEY(u1);
				u1   = searchTree[fkey];
				ancesters[i-1-m] = aSECONDKEY(u1);
			}


//#ifdef TIMERFLAG
//			Timers.start_timer(TIMER2);
//#endif
			int outdegree = 0;
			for( n = nodeTree[j]; n < i && !outdegree; n++){
				
				//get neighbors
                vector<GNSIZE> * neighbor_vec = SAdjList[ancesters[n]];
	
				// step 2. scan for neighbors that havn't be searched before
				for( k=0; k< neighbor_vec->size(); k++){
					v = (*neighbor_vec)[k];
                    ancesters[i] = v;

					//string *buf = matrix[v] ;
					bool flag = true;

					//populate the cost
					int tcode = 1, mii =0;
					//cout << "code is " << code << "\n";
					if( find(ancesters.begin(), ancesters.begin()+i, v) 
                                               != (ancesters.begin()+i))
						flag = false;
					else{
					   for( m=0; m<=i && tcode==1 ; m++){
						  w = ancesters[m];
					      int refkey = getLabel(v, w);
						  if(  refkey != EMPTY_EDGE){
	
						     if( refValues[mii] < aMAKEKEY(i-m,  refkey )){
						        if( ((i-aFIRSTKEY(refValues[mii])) > bi && i == bl) || i > bl) tcode = 0;
						        else tcode = -1;  //no need to further search
									//cout << "value " << hex << refValues[mii] << " " << aMAKEKEY(i-m, (*buf)[w]) << dec << endl;
							 }
						     else if( refValues[mii] > aMAKEKEY(i-m, refkey) )tcode = 2;
						     mii++;
					       }
					   }
					}
					//cout << "code is " << code  << " " << tcode <<" \n";
					//cout << "flag is " << flag  << " " << v << "\n";

					if( flag ) {
						
						if( tcode < 0 ) return tcode;
						else if( tcode ==0 ) code =0;

						outdegree++;

						if( tcode == 1 ){
							searchTree.push_back(aMAKEKEY(j, v));
							nodeTree.push_back(n);
						}
					} //if( flag) know the result for sure
				}
			}
//#ifdef TIMERFLAG
//			Timers.stop_timer(TIMER2);
//#endif
		} //for (j = uplimit1; j < uplimit2; j++)  
		uplimit1 = uplimit2;
	}

	//assert(false); // never reach here
	return code;
}

void AdjMatrix::print(int pi, string f1, string f2)
{
        //output the pattern
        ofstream of1(f1.c_str(), ios::app);
        int i, j;

        for ( i = 0; i < size(); i++){
          GLTYPE c = getLabel(i, i);
          of1 << "node " << dec << pi << " " << i << " " << (c-EMPTY_EDGE) << endl;
        }
        of1.close();

        //output the pattern
        ofstream of2(f2.c_str(), ios::app);
        for ( i = 1; i < size(); i++){
          for( j=0; j< i; j++){
             GLTYPE c = getLabel(i, j);
             if( c != EMPTY_EDGE)
                 of2 << "edge " << dec << pi << " " << i << " " << j << " " <<
				 (c-EMPTY_EDGE) << endl;
          }
        }
        of2.close();
}


void AdjMatrix::print(int pi, string f1, string f2, string f3, map<int,int> & g_reindex)
{
    ofstream of3(f3.c_str(), ios::app);
    //of3 << "this" << endl;
	int density = size()*(size()-1)/2-eSize;
	of3 << "sup: " << totalSup << " bgf: " << bgf << " s: " << size() << " d: " << density;
 
#if (SAVERESULTS == 2)
	GRAPHS1 ::iterator ip;
	of3 << " tids: " ;
    for( ip = TIDs.begin(); ip != TIDs.end(); ip++){
		if (g_reindex.count(*ip) == 0 ) error("now graph mapping is found");
		
        of3 << dec << (g_reindex[*ip]+1) << ":" << 1  << " " ;
	}
#endif

	//cout << " header " << f3  << endl;
	of3 << endl;
    of3.close();
    print(pi, f1, f2);
}

#ifdef SVAERESULTS
void AdjMatrix::setGraphs(vector<occur*> & occs1)
{
	occs = occs1;
    for(int i=0; i< occs.size(); i++){
       int gi = occs[i]->getGraph();
       if( graphs.count(gi) )
         graphs[gi]++;
       else
         graphs[gi] = 1;

#ifdef SAVERESULTS
	   if( SAVERESULTS == 1 ) delete occs[i];
#endif
	}
	occs.clear();
}

void AdjMatrix :: setTIDs(GRAPHS1 & g,  int tsize)
{
	totalSup = g.size();
	int temp=0, j=0;
	for(int i=0; i< tsize; i++){
		if( j< totalSup && i == g[j] ){
			j++; temp++;
		}
		if( (i+1) % 32 == 0){ TIDs.push_back(temp); temp = 0; }
		else	temp = temp << 1;
	}
	TIDs.push_back(temp);
}
#endif

void AdjMatrix::print()
{
	//output the pattern
	int i, j;
	cout << "adj matrix: " << endl;
	for ( i = 0; i < size(); i++){
	  //string * l = matrix[i];
	  char c;
	  for(j=0; j< size(); j++){
	    c = getLabel(i,j);
	    cout << (c-EMPTY_EDGE) << " ";
	  }
	  cout << endl;
	}
	//cout << "occurrences size " << ocsize << endl;
	
	if( SAdjList){
		cout << "Adj list are " << endl;
		for(i=0; i< size(); i++){
			vector<GNSIZE> * t = SAdjList[i];
	
			cout << "at node " << i << "\t";
			for( j=0; j< t->size(); j++){
				cout << (*t)[j] << "  " ;
			}
			cout << "\n";
		}
		cout << endl;
	}

	
#if ( SAVERESULTS == 2)
	cout << "tid set " << endl;
	for( i=0; i<TIDs.size(); i++){
		cout << hex << TIDs[i] << " " << dec << endl;
	}
#endif
}

int AdjMatrix:: lastbit(){
  int s =size(), i;
  for( i=(size()-2); (i>=0 && getLabel(s-1, i) == EMPTY_EDGE); i--){
     ;
  }
  if( i< 0 ){
	  print();
	  error("i in last bit should not be negative\n");
  }
  return i;
}
//find the second to last bit's index
void AdjMatrix :: last2bitil(int & bi, int & bl)
{	
	int s = size();
	char c;
	int fi = 0, i;

	for( i= (s-2); i >= 0; i--){
		c = getLabel(s-1, i);
	    //c = lr[i];
		if( c!= EMPTY_EDGE ){
			if( ++fi == 2) { bi = i, bl = (s-1); return ;}
		}
	}

	//not in the last row
	if( fi == 1){
		for( i= (s-3); i >= 0; i--){
		    c = getLabel(s-2, i);
			if( c!= EMPTY_EDGE ){
				bi = i, bl = (s-2); 
				return;
			}
		}
	}
	
	cout << " i " << i ;
	error("could not find the second to last bit");
}

void AdjMatrix::addNode(char label)
{
	int s = size(), i;

	for( i= 0; i < s; i++){
		pushLabel(EMPTY_EDGE);
	}
	pushLabel(label);
	matrixSize++;
}


void AdjMatrix::addEdge(int i, int j, GLTYPE label)
{
#ifdef SAFEMODE
  if( i> size() || j > size() )
    cout << "i " << i << " j " << j << " size " << size() << endl;
  
	if( getLabel(i,j) != EMPTY_EDGE ) 
		error("add to a none empty edge\n");
#endif

	assignLabel(i, j, label);
	eSize++;
}

int  AdjMatrix :: esize(){

	//cout << "search esize" << endl;
	if( eSize ) return eSize;
        //cout << "search esize again" << endl;

	int k = 0; 
	for(int i=0; i<size(); i++){
		for(int j=0; j<i; j++){
			if( getLabel(i, j) != EMPTY_EDGE)
				k++;
		}
	}
	eSize = k;
	return k;
}
void AdjMatrix::removeEdge(int i, int j)
{
	assignLabel(i, j, EMPTY_EDGE);
}

string AdjMatrix:: canonicalForm()
{
	string s; 
	int s1 = size();

	for(int i=0; i< s1; i++){
		for(int j=0; j< i; j++){
			s += getLabel(i,j);
		}
	}
	return s;
}

/*
  return true if this is a subgraph of m
*/
bool AdjMatrix:: isSubgraphOf(AdjMatrix * m, int f){

#if (defined SAVERESULTS && SAVERESULTS == 3) 
	//cout << "why here\n" ; 
	if( m->esize() < esize() ){ /*cout << "===small size" << endl;*/  return false; }   //m's size must be large
	if( m->totalSup > totalSup ) { /*cout << "====larger suppport" << endl;*/ return false;  } //m's support value must be small
	if( m->TIDs.size() != TIDs.size() ) error("mismatching size");
	//check the embedding relation
	for(int i=0; i<TIDs.size(); i++){
		if( (m->TIDs[i] & TIDs[i]) != m->TIDs[i] ){
			/*cout << "====mis matching tids" << endl;*/
			return false;    //m's tid set is must be subsumed by subgraph's tid set
		}
	}

    UllmanIsomorphism UI = UllmanIsomorphism(this, m);
	//cout << "before returning" << endl;
    return UI.Ullman_isomorphism();
#else
        //cout << "should be here\n"; 
        if( m->esize() < esize() ){  return false; }
	UllmanIsomorphism UI = UllmanIsomorphism(this, m);
	if( f ){ //check automorphism 
	   //cout << "check automorphism" << endl;
	   if( m->esize() != esize() ){ return false; } 
           UllmanIsomorphism UI1 = UllmanIsomorphism(m, this);
	   bool flag =  UI.Ullman_isomorphism() && UI1.Ullman_isomorphism();
	   //cout << flag << endl;
	   return flag;
	}
	else return UI.Ullman_isomorphism();
#endif
}

//for validation checking, used in linux only
#ifdef VFLIB
void AdjMatrix:: setGraphForMatch(){

	if( myg == NULL){
		//cout << "buliding graph for iosmorphism testing" << endl;
		int l , s = size();
		int i, j;

		ARGEdit ed1;
		for(i=0; i<s; i++){
			l = (int) getLabel(i, i);
			ed1.InsertNode(new Point(l)); // The inserted node will have index i.                   
		}

		 // Insert the edges
		for(i=0; i<s; i++)
			for(j=0; j<s; j++){
				if (i!=j){
					l = (int) getLabel(i, j);
					//if( l != EMPTY_EDGE )
						ed1.InsertEdge(i, j, new Point(l)); // NULL stands for no sem. attribute.
			}
		}

		//set up now
		myg = new ARGraph<Point, void>(&ed1);
		myg->SetNodeDestroyer(new PointDestroyer());
		myg->SetNodeComparator(new PointComparator());
		myg->SetEdgeComparator(new PointComparator());
	}
}

/*
  return true if m is a subgraph of this
*/
bool AdjMatrix:: isSubgraph(AdjMatrix * m){

	setGraphForMatch();
	m->setGraphForMatch();
	
    // Now the Graph can be constructed...
	//UllSubState s0(&g2, myg);
	VF2SubState s0(m->getARGraph(), myg);

	int s = m->size(), i;
	int js = size() > s? size() : s;
	node_id ni1[js], ni2[js];
	int n;
	bool flag;

    if( match(&s0, &n, ni1, ni2) ){

#ifdef DEBUG5
		for( i=0; i<n; i++){
			cout << "graph small vs bigger mapping: " << ni1[i] << " -> " << ni2[i] <<endl;
		}
		cout << "end" << endl;
#endif
		return true;
	}
	else
		return false;

}

bool AdjMatrix:: isIsomorphic(AdjMatrix * m)
{
	setGraphForMatch();
	m->setGraphForMatch();

	VF2State s0(m->getARGraph(), myg);

	int s = m->size(), i;
	int js = size() > s? size() : s;
	node_id ni1[js], ni2[js];
	int n;
	bool flag;

    if( match(&s0, &n, ni1, ni2) ){

#ifdef DEBUG5
		for( i=0; i<n; i++){
			cout << "graph small vs bigger mapping: " << ni1[i] << " -> " << ni2[i] <<endl;
		}
		cout << "end" << endl;
#endif
		return true;
	}
	else
		return false;
}

bool AdjMatrix:: my_visitor(int n, node_id ni1[], node_id ni2[], void *usr_data)
{
	//FILE *f = (FILE *)usr_data;
     int * s = (int*)  usr_data;
     
	
    // Prints the matched pairs on the file
    int i;
	(*s)++;
    /*for(i=0; i<n; i++)
      fprintf(f, "(%hd, %hd) ", ni1[i], ni2[i]);
    fprintf(f, "\n");*/

    // Return false to search for the next matching
    return false;
}

int AdjMatrix:: findAllIsomophisms(AdjMatrix * m)
{

	int lsum =0;
	setGraphForMatch();
	m->setGraphForMatch();
	
	VF2SubState s0(m->getARGraph(), myg);
	match(&s0, &my_visitor, (void*) &lsum);
	
    //match(&s0, &my_visitor, NULL);
	//cout << "lsum " << lsum << endl;
	return lsum;
}

#endif

//setup the Adj list
void AdjMatrix :: buildAdjList()
{
	int s = size();
	SAdjList = new vector<GNSIZE> *[s];

	for(int i = 0; i< s; i++){
	
		vector<GNSIZE> * t = new vector<GNSIZE>;
		SAdjList[i] = t;

		for(int j =0; j< s; j++){
			if( getLabel(i,j) != EMPTY_EDGE && j!= i){			
				t->push_back(j);
			}
		}
	}

}

