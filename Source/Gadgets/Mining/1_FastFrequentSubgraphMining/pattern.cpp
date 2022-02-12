#pragma warning (disable:4786 4018 4267)
#include <vector>
#include <fstream>
#include <algorithm>
#include <set>
using namespace std;

#include "pattChild.h"

gBase    * pattern ::gb = NULL;
MyTimer  * pattern ::Timers = NULL;
Register * pattern ::Registers = NULL;
DLONG pattern ::totalPatterns = 0;

bool pattern :: denseEnough(int d, int sizeup, bool final){
    
    int s = M->size();
    if( s == 1 || s== 2) return true;
    return M->size() <= sizeup && M->denseEnough(d, final);
}

void pattern :: initPattern(gBase *gb1, MyTimer & t, Register & mr){
		gb = gb1; Timers = &t; Registers = &mr; 	
}

int pattern:: isCarnonicalForm(void) { 

	stage = 1;
	
#ifdef TIMERFLAG
	Timers->start_timer(TIMER4);
#endif

	fixM(); 

#ifdef TIMERFLAG
	Timers->stop_timer(TIMER4);
#endif
	

#ifdef TIMERFLAG
	Timers->start_timer(TIMER2);
#endif

	int flag = M->isCanonicalForm();  
	
#ifdef TIMERFLAG
	Timers->stop_timer(TIMER2);
#endif

	return flag;

}

void pattern:: sprint(int level)
{
	cout << "the adj matrix is" << endl;
	if( M )
		M->print();
	else
		cout <<"matrix is not initialized" <<endl;
	
	cout << "pattern optimal "  << optimal << " cost " << hex << mycost.intValue() << dec << endl;

	if( level & 2  ){
		cout << "the coccs structures are: size: " << coccs->size() << endl;
		
		COCCS :: iterator ip1;
		for(ip1 = coccs->begin(); ip1 != coccs->end(); ip1++){
			DLONG key = *ip1;
			cout << "index is " << (FIRSTKEY(key) ) << " node is " << 
				 (THIRDKEY(key)) << " graph " << (SECONDKEY(key)) << endl;
		}
	}

	if( level & 4  ){
		cout << "there are total: " << child.size() << " children in this pattern " << endl;
		
		vector<pattern *> :: iterator ip1 ;
		for(ip1 = child.begin(); ip1 != child.end(); ip1++){
			(*ip1)->sprint(5);
		}
	}
}

void pattern:: print()
{
	cout << "the adj matrix is" << endl;
	if( M )
		M->print();
	else
		cout <<"matrix is not initialized" <<endl;

	//int i;
	cout << "the index structure are" << endl;
	COCCS :: iterator ip1;

	for(ip1 = coccs->begin(); ip1 != coccs->end(); ip1++){
		DLONG key = *ip1;
		cout << "index is " <<  (int)(FIRSTKEY(key) ) << " node is " << (int) (THIRDKEY(key)) << endl;
	}

	cout << "cost is " << hex << mycost.intValue() << endl;
}

void pattern:: fsm1( vector<AdjMatrix * > & result, int freq, int level, int density, string nodef,
					string edgef, string header, int sizelimit, int sizeup, map<int, int> & gIndex) {

	CHLTYPE :: iterator ip;
	//cout << "start search total frequent node :" << child.size() << endl;

	vector<pattern*> wlist;
	int totalnum = 0;

	for( ip = child.begin(); ip!=child.end(); ip++){
		 //cout << "search with support: " << (*ip)->sup() << endl;

#ifdef DEBUG3
		cout << "start search" << endl;
		(*ip)->sprint(3);
#endif

		//(*ip)->print();
		(*ip)->fsm2(result, freq, level+1, totalnum, density, nodef, edgef, 
			header,  sizelimit, sizeup, gIndex);
		
		//delete the edge;
		delete (*ip);

	}

	//cout << "=====FFSM Results Are =====" << endl;
        cout << "Total frequent patterns: " << dec << totalnum << "\t";

#ifdef DEBUG8
	edbLocal->print();
#endif
	
}


void pattern:: fsm2( vector<AdjMatrix * > & result, int freq, int level, int& totalnum, 
					int density, string nodef, string edgef, 
					string header, int sizelimit, int sizeup, map<int, int> & gIndex){

	
#ifdef DEBUG2
	cout << "search at level " << level << endl;
	sprint(3);	
#endif

	stage = 4;
	if( level >= 2){
		totalnum++;
		if( totalnum % 1000 == 0 ) cout <<".";
		if( totalnum % 10000 == 0 ) cout << totalnum  << "\n";

#ifdef REGISTERFLAG
		if( M->isPath() ) Registers->peg(REGISTER1);
#endif

		M->setTotalSup(sup());
#ifdef SAVERESULTS
/*
#if      (SAVERESULTS == 3)
			//vector<occur*> loccs;
			//recordOccs(loccs);
			//M->setGraphs(loccs);
			*/
#if (SAVERESULTS == 2)
		M->setTIDs(gids);
#endif
		result.push_back(M);
#endif
	}

	//should output the results
	if( result.size() >= BUFFERSIZE){
	   cout << "output results " << totalnum << endl;
	           
	   for(int ri=0; ri< result.size(); ri++){
			bool dflag = true;

#ifdef DENSITY
			if(!result[ri]->denseEnough(density, true) && M->size() <= sizeup) dflag = false;
#endif
#ifdef SAVERESULTS
			if( dflag &&  result[ri]->size()>= sizelimit) 		 
	   			result[ri]->print(totalPatterns++, nodef, edgef, header, gIndex);
			delete result[ri];
#endif

	   }	
	   result.clear();
	}
       	
	pattern * pi, *pj;
    bool maximal;

	for( int i = 0; i< child.size(); i++){
		pi = child[i];
		vector<pattern*> outChild;

#ifdef DEBUG3
	    for( int lii = 0; lii< level; lii++)
			cout << "  ";
		cout << "cost and child size " << hex << pi->mycost.intValue() << dec << "  " <<  pi->child.size() << endl;
#endif

		if( (pi->optimal & 1) && pi->denseEnough(density, sizeup, false) ){  //the pattern must be in carnonical form for subsequent search

			maximal = true;

			if( size != pi->size ) {	//only one edge in the last row

#ifdef TIMERFLAG
	            Timers->start_timer(TIMER3);
#endif
#ifdef TIMERFLAG
		        Timers->start_timer(TIMER6);
#endif
				maximal = pi->proposeOuterChild(freq, outChild);
				
#ifdef TIMERFLAG
		       Timers->stop_timer(TIMER6);
#endif
#ifdef TIMERFLAG
	          Timers->stop_timer(TIMER3);
#endif
			}

#if (defined MAXIMALPATTERN && (MAXIMALPATTERN & 1) )
			if( maximal ){
#endif

#ifdef DEBUG2
				cout << "children are (proposed by outer): " << endl;
				int k = 0;
				for( ; k< pi->outChild.size(); k++)	(pi->child[k])->sprint();
#endif

				for(int j = i; j< child.size(); j++){
					pj = child[j];
					pi->join(pj, freq);
				}

#ifdef DEBUG2
				cout << "completed child proposing" << endl;
				cout << "children are: " << endl;
				for( k=0; k< pi->child.size(); k++)		(pi->child[k])->sprint();			
#endif
				//post processing
				(pi->child).insert((pi->child).end(), (pi->newchild).begin(), (pi->newchild).end());
				(pi->child).insert((pi->child).end(), outChild.begin(), outChild.end());
				pi->fsm2( result,freq, level+1 , totalnum, density, nodef, edgef,
					header, sizelimit, sizeup, gIndex);

#if (defined MAXIMALPATTERN && (MAXIMALPATTERN & 1) ) 
			}
#endif
			//if( pi->size == 2) gb->removeEdge(pi->mycost.intValue());
		}
		delete pi; /*definitely need to be taken care of*/
	}
}

	
void pattern:: constructM (void){
	char node1 = (char) mycost.index();
	char node2 = mycost.nlabel();
	char edge  = mycost.elabel();
	M = new AdjMatrix();
	M->addNode(node1);
	M->addNode(node2);
	M->addEdge(0, 1, edge);
	M->buildAdjList();
}

void pattern:: fixM(void){

	int index = mycost.index();
	GLTYPE nlabel;

	int s1 = par->getSize();
	if( s1 == size){
		index -= s1;
		nlabel = EMPTY_EDGE;
	}
	else
		nlabel = mycost.nlabel();

	index = size-2-index;
	//M = new AdjMatrix(M, index, mycost.elabel(), nlabel );
	M = new AdjMatrix(M, size-1, index, mycost.elabel(), nlabel );
	stage &= 1;
}	

void pattern::setChild(EDGES & cm, EDGES & nm)
{
	EDGES::iterator ip;
	
	vector<int> u;

	ip = cm.begin();
	int v = COSTINDEX(ip->first), v1, k;
	u.push_back(ip->first);
	ip++;

	//search for those edges sharing the start node
	for( ; ip!= cm.end() ; ip++){
		k = ip->first;
		v1 = COSTINDEX(k);	
		if( v1 == v)
			u.push_back(k);
		else{
			child.push_back(setUpChild(u, cm, nm, v));
			u.clear();
			u.push_back(k);
			v = v1;
		}
	}
	child.push_back(setUpChild(u, cm, nm, v));
}

pattern * pattern :: setUpChild(vector<int> u, EDGES & edges, EDGES & gnodes, int nlabel)
{
	pattern * p1 = new pattern(this, 0, 1), * p2;
	CHLTYPE & c1 = p1->child;

	TRA lmapper;
	//populate the lmapper
	OCCS * es = gnodes[nlabel];
	OCCS::iterator oip = es->begin();
	int index =0;
	for(; oip != es->end(); oip++){
		occur * oc = *oip;
		DLONG key = MAKEKEY(INULLINDEX, oc->getGraph(), oc->getNodes(0));
		if( lmapper.count(key) ) error("duplicated key");
		lmapper[key] = index++;
		p1->coccs->push_back(key);
	}

	int i;
	for(i=0; i< u.size(); i++){
		es = edges[u[i]];
		p2 = new pattern(p1, u[i], 2);
		p2->constructM();
		setUpEdges(p2, es, lmapper);

#ifdef DEBUG5
		cout << "p2 " << endl;
		p2->print();
#endif
		c1.push_back(p2);
	}


#ifdef DEBUG5
	cout << "p1 " << endl;
	p1->print();
#endif

	return p1;
}

void pattern :: setUpEdges(pattern * p2, OCCS * es, TRA & tra1)
{
	//set up for the p2(edge)
	COCCS * cos2 = p2->coccs;	
	GRAPHS1 & gid2 = p2->gids;

	int gi, d1, d2, t1, j;
	DLONG key;

	//every edge is frequent and in carnonical form
	p2->optimal = 1;

	for(j=0; j< es->size(); j++){
		occur * noc = (*es)[j];
		gi = noc->getGraph();
		d1 = noc->getNodes(0), d2 = noc->getNodes(1);

		key = MAKEKEY(INULLINDEX, gi, d1);
		
		if( !tra1.count(key) ) error("the hash at p1 is not setting up");
		t1 = tra1[key];
		
		cos2->push_back(MAKEKEY(t1, gi, d2));

		//for(int m=0; m< gid2.size(); m++){
		//	if( gid2[m] > gi) error("why redundant graphs");
		if( !gid2.size() || gi != gid2[gid2.size()-1] ) gid2.push_back(gi);
		//}
		//if( !gid2.size() ) gid2.push_back(gi);
	}
}

void pattern :: join(pattern * p2, int freq){
	cost1 c2 = p2->mycost;
	int index1 = mycost.index(),node1 = mycost.nlabel();
	int index2 = c2.index(),	node2 = c2.nlabel();
	
	//check the intersection
	int sum=0;
	GRAPHS1 ::iterator iip1 = gids.begin(), iip2 = p2->gids.begin();
    while( iip1 != gids.end() && iip2 != p2->gids.end() ){
		if( *iip1 == * iip2){ iip1++; iip2++; sum++;}
		else if( * iip1 > *iip2) *iip2++;
		else iip1++;
	}

	//no way to be frequent
	if( sum < freq ) return;

	//test carnonical form first or occurrence first
	bool mf = true;
	if( size > 5) mf = false;

	//get the type of the join
	pattern * p3 = NULL;;
	int t = par->getSize(), type;

	//special treatment only for the first edge
	if( size == 2) index1 = index2 = 0;

	if( index2 >= t) {  //inner + inner  //notice it is okay for look back
		if( index1 != index2  && node1  == node2 ){
			p3 = new pattern(this, c2.intValue(), size);
			type = 1;		
		}
	}
	else if( index1 >= t){ //inner + outer
		p3 = new pattern(this, c2.intValue(), size+1);
		type = 2;
	}
	else{
		if( index1 != index2 ){
#ifndef TREEBASEDSCAN
			if( node1 == node2 ){
				int c = ((size+index2) << 16 ) | (c2.elabel() << 8 ) | node2;
				p3 = new pattern(this, c, size);
				type = 3;
				if ( join1(p2, p3, type, freq, mf) ) child.push_back(p3);
				else {
					delete p3 ; /*definitely need to be taken care of*/
				}
			}
#endif
		}
		else{
			if( mycost.intValue()> c2.intValue() && p2->optimal){
				//the only case we might switch the order
				p2->join(this, freq);
			}
		}
		
		
		int c = ((1+index2) << 16 ) | (c2.elabel() << 8 ) | node2;
		p3 = new pattern(this, c, size+1);
		type = 4;
	}

	//test for inner join of occurrences
	if ( p3 )
		if(join1(p2, p3, type, freq, mf) ){
			if( type == 4)
				newchild.push_back(p3);
			else
				child.push_back(p3);
		}
		else{
			delete p3; /*definitely need to be taken care of*/
		}
	}

bool pattern :: join1(pattern * p2, pattern * p3, int type, int freq, bool mf)
{


	//test carnonical form is mf is marked
	int t;
	bool carno = false;
	bool support = false;

#ifdef REGISTERFLAG
	if( !support){
		Registers->peg(REGISTER1);
	}
#endif

	if( mf) {
		t = p3->isCarnonicalForm();
		if( t< 0) {
			carno = false;
		}
		else{
			if( t> 0 ) p3->optimal = 1;
			carno = true;   //for optimal and suboptimal
		}
	}

	//get the inner join of occurrences
#ifdef TIMERFLAG
	Timers->start_timer(TIMER3);
#endif

	if( !mf || carno ){
		if( (type & 1)  ){
			support = propose_inner(p2, p3, freq);
		}
		else if( type  == 2){
			support = propose_inner_outer(p2, p3, freq);
		}
		else{
			support = propose_outer2(p2, p3, freq);
		}
	}
	
#ifdef REGISTERFLAG
	if( !support){
		Registers->peg(REGISTER2);
	}
#endif

#ifdef TIMERFLAG
	Timers->stop_timer(TIMER3);
#endif

	//if frequent, to see whether in carnonical form or not
	if( !mf && support) {
		t = p3->isCarnonicalForm();

		if( t< 0){ carno = false;  optimal |= 2; } //added by Luke for Maximal
		else{
			if( t> 0 ) p3->optimal = 1;
			carno = true;
		}
	}

#ifdef DEBUG2
		cout << "type : carnonical form and support test: " << type << "  " << t << " " << support << endl;
#endif

	//passed everything
	if( !carno || !support ) return false;
	else{
		return true;
	}
}

//for inner + inner or outer + outer `case 1
bool pattern :: propose_inner(pattern * p2, pattern * p3, int freq)
{
	//set up
	COCCS :: iterator ip1, ip2;
	COCCS * cos2 = p2->coccs;
	ip1 = coccs->begin(); ip2 = cos2->begin();

	COCCS * cos3 = p3->coccs;		
	GRAPHS1 & gid3 = p3->gids;

	int index =0; 
	bool flag = true;

	while(  ip1 != coccs->end()  && ip2 != cos2->end() ){
	
		if( *ip1 == *ip2 ){
		
			int gi = SECONDKEY(*ip1);
			cos3->push_back(MAKEKEY(index, gi, GNULLINDEX));
			//graph3->insert(gi);
			
			//update
			if( !gid3.size()  || gi != gid3[gid3.size()-1]) gid3.push_back(gi);
			//end

			ip1++; index++;	ip2++;
		}
		else if( *ip1 > *ip2 ){
			ip2++; flag = false;
		}
		else{
			ip1++; index++;
		}
	
	}

#if (defined MAXIMALPATTERN && (MAXIMALPATTERN & 1) )
	if( flag ){
		p2->optimal=0; //added by Luke for Maximal
#ifdef DEBUG11
		cout << "this is one absorbing "<< endl;
#endif
	}
#endif

	if( p3->sup() >= freq) {
		p3->stage = p3->stage | 2;
		return true;
	}
	else return false;
}

//for inner + outer 
bool pattern :: propose_inner_outer(pattern * p2, pattern * p3, int freq)
{
	//set up
	COCCS :: iterator ip1, ip2, ips, ipe;
	COCCS * cos2 = p2->coccs;
	ip1 = coccs->begin(); ip2 = cos2->begin();

	COCCS * cos3 = p3->coccs;		
	GRAPHS1 & gid3 = p3->gids;

	int index = 0, t, k;
	bool flag = true;

	while( ip1 != coccs->end()  && ip2 != cos2->end() ){

		t = FIRSTKEY(*ip1);
		k = FIRSTKEY(*ip2);

		if( t == k ){

			//find the ipe
			ips = ipe = ip2;
			while( ipe != cos2->end() &&  ( FIRSTKEY(*ipe) == t ) ) ipe++;

			int gi = SECONDKEY(*ip1);
			//graph3->insert(gi);
			
			//update
			if( !gid3.size() || gi != gid3[gid3.size()-1]) gid3.push_back(gi);
			//end

			for( ; ip2 != ipe; ip2++){
				cos3->push_back(MAKEKEY(index, gi, THIRDKEY(*ip2)));
			}

			ip1++;	index++; 
			
		}
		else if( t > k ){
			ip2++;  flag = false; 
		}
		else{
			ip1++; index++;
		}
	
	}

#if ( defined MAXIMALPATTERN && (MAXIMALPATTERN  & 1) )
	if( flag ){
		p2->optimal=0; //added by Luke for Maximal
#ifdef DEBUG11
		cout << "this is one absorbing in inner_outer "<< endl;
#endif
	}
#endif

	if( p3->sup() >= freq) {
		p3->stage = p3->stage | 2;
		return true;
	}
	else return false;
}

//for outer + outer case 2
bool pattern :: propose_outer2(pattern * p2, pattern * p3, int freq)
{
	COCCS :: iterator ip1, ip2, ips, ipe;
	COCCS * cos2 = p2->coccs;
	ip1 = coccs->begin(); ip2 = cos2->begin();

	COCCS * cos3 = p3->coccs;
	GRAPHS1 & gid3 = p3->gids;

	int index = 0, t, k;

	while( ip1 != coccs->end()  && ip2 != cos2->end() ){

		t = FIRSTKEY(*ip1);
		k = FIRSTKEY(*ip2);

		if( t == k ){

			//find the ipe
			ips = ipe = ip2;
			while( ipe != cos2->end() &&  ( FIRSTKEY(*ipe) == t ) ) ipe++;

			int gi = SECONDKEY(*ip1);

			bool flag = false;
			for( ; ip2 != ipe; ip2++){
				if( *ip1 !=  *ip2){
					cos3->push_back(MAKEKEY(index, gi, THIRDKEY(*ip2)));
					flag = true;
				}
			}

			//if( flag ) graph3->insert(gi);
			if( flag ) 	
				//update
				if( !gid3.size() || gi != gid3[gid3.size()-1]) gid3.push_back(gi);
				//end
	
			ip1++;	index++; ip2 = ips;
			
		}
		else if( t > k ){
			ip2++; 
		}
		else{
			ip1++; index++;
		}
	
	}

	if( p3->sup() >= freq) {
		p3->stage = p3->stage | 2;
		return true;
	}
	else return false;
}

bool pattern:: proposeOuterChild(int freq, vector<pattern*> & outChild){
	 	 
	COCCS :: iterator ip = coccs->begin();
	int i, u, gi, x, key, level;
	map<int, pattern * , greater<int> > children;
	int elabel, nlabel;
	AdjMatrix * gbi;
	pattern * p3;

    //for update
	int edgeSize = M->esize(), j;
	pattern * currp = this;
	vector<COCCS*> coccslist;
	vector<int> gindices(edgeSize, INULLINDEX);
	vector<short> checkpoint(edgeSize+1, 0);

//	DLONG mykey;
	int fkey, skey, fkey1;

	//cout << "start search" << endl;
#if  (defined MAXIMALPATTERN && (MAXIMALPATTERN & 1) )
	initialEdgeCan();
#endif

	for( i=0; i<= edgeSize; i++){
		//update
		if( THIRDKEY((*(currp->coccs))[0]) != GNULLINDEX )
			checkpoint[i] = 1;
		//end

		coccslist.push_back(currp->coccs);
		currp = currp->par;
	}

	vector<GNSIZE> filter(size, 0);
	
	for( i=0; ip != coccs->end(); i++ , ip++){

		fkey = FIRSTKEY((*coccs)[i]), skey = THIRDKEY((*coccs)[i]);
		gi = SECONDKEY((*coccs)[i]);
		u = skey;
		level = 1;
		filter[0]= skey;

		for( j=0; j< edgeSize && fkey != gindices[j]; j++){
			gindices[j] = fkey;
			fkey1 = FIRSTKEY((*(coccslist[j+1]))[fkey]);
			if( checkpoint[j+1])  
				filter[level++]= THIRDKEY((*(coccslist[j+1]))[fkey]);
			fkey = fkey1;
		}


#ifdef DEBUG10
		vector<GNSIZE> :: iterator  fip = filter.begin();
		cout << "the nodes in the filter " << i << endl;
		for( fip = filter.begin(); fip != filter.end(); fip++){
			cout << (int) *fip << "\t" ;
		}
		cout << endl;
#endif

		//get the neighbors
		gbi = gb->graph(gi);
		
#if (defined MAXIMALPATTERN && (MAXIMALPATTERN & 1) )
		fillMoreEdge(filter, gbi);
#endif
		//for debugging only
		vector<GNSIZE> * nei = gbi->getNeighbors(u);
		
		//propose child
		for( j=0; j< nei->size(); j++){
			x = (*nei)[j];

			if( (find(filter.begin(), filter.end(), x) == filter.end() ) && 
					gbi->getLabel(u, x) != EMPTY_EDGE){

				elabel = gbi->getLabel(u, x);
				nlabel = gbi->getLabel(x, x);
				key = elabel << 8 | nlabel;

			
#ifdef REGISTERFLAG 
				//cout << "* " << endl;
				Registers->peg(REGISTER3);
#endif
				if( children.count(key) ){
					p3 = children[key];
				}
				else{
					p3 = new pattern(this, key, size+1);
					children[key] = p3;
				}
		
				if( !p3->gids.size() || gi != p3->gids[p3->gids.size()-1]) p3->gids.push_back(gi);
				p3->coccs->push_back(MAKEKEY(i, gi, x));
			}
		}
	 }


#if (defined MAXIMALPATTERN && (MAXIMALPATTERN & 1) )
	if( !isInnerMaximal() ){
		map<int, pattern * , greater<int> > ::iterator cit = children.begin();
		for( ; cit != children.end(); cit++)  delete cit->second;
//#ifdef DEBUG12
		//cout << "remove one "<< endl;
//#endif
		return false;
	}
#endif


	//scan children to see whether it is frequent and optimal 
	scan(children,freq, outChild);
	return true;
}

void pattern:: scan(map<int, pattern*, greater<int> > & cand, int f, vector<pattern*> & result)
{

	map<int, pattern*, greater<int> > :: iterator ip;
	
	for(ip = cand.begin(); ip!= cand.end(); ip++){
		pattern *  p = ip->second;

#ifdef REGISTERFLAG 
		Registers->peg(REGISTER8);
#endif
		if( p->sup() < f ) { 

#ifdef REGISTERFLAG 
			Registers->peg(REGISTER7);
#endif

#ifdef DEBUG2
			cout << "failed support test " << endl;
#endif
			delete p;  /*definitely need to be taken care of*/
		}
		else {
			int t = p->isCarnonicalForm();

#ifdef DEBUG2
			cout << "the scan test and t value: " << t << endl;
#endif
			if( t < 0 ) {

#ifdef REGISTERFLAG 
				Registers->peg(REGISTER6);
#endif
				delete p;	 /*definitely need to be taken care of*/
				optimal |= 2; /*added by Luke for maixmal*/
			}
			else {
				p->optimal = t;
				//if( t == 1)   p->fixOccs();
				result.push_back(p);

#ifdef REGISTERFLAG 
				if( !t )
					Registers->peg(REGISTER5);
				else			
					Registers->peg(REGISTER4);
#endif
			}
		}

	}
}

//save results
#ifdef SAVERESULTS
void pattern:: recordOccs(vector<occur*> & loccs)
{
	//for update
	//M->print();
	int edgeSize = M->esize(), j;
	pattern * currp = this;
	vector<int> gindices(edgeSize, INULLINDEX);
	vector<short> checkpoint(edgeSize+1, 0);
	vector<COCCS*> coccslist;

//	DLONG mykey;
	int fkey, skey, fkey1, gi, u, level, i;

	//cout <<"edge size " << edgeSize << endl;

	for( i=0; i<= edgeSize; i++){
		if( THIRDKEY((*(currp->coccs))[0]) != GNULLINDEX )
			checkpoint[i] = 1;
		coccslist.push_back(currp->coccs);
		currp = currp->par;
	}

	/*
	cout << "check points" << endl;
	for(i=0; i< checkpoint.size(); i++)
		cout << checkpoint[i] << " ";
	cout << endl;
	*/

	vector<GNSIZE> filter(size, 0);
	COCCS :: iterator ip = coccs->begin();

	for( i=0; ip != coccs->end(); i++ , ip++){

		fkey = FIRSTKEY((*coccs)[i]), skey = THIRDKEY((*coccs)[i]);
		gi = SECONDKEY((*coccs)[i]);
		u = skey;
		level = 0;

		if( checkpoint[0] ){
			filter[0]= skey;
			level++;
		}

		for( j=0; j< edgeSize && fkey != gindices[j]; j++){

			gindices[j] = fkey;
			fkey1 = FIRSTKEY((*(coccslist[j+1]))[fkey]);
			if( checkpoint[j+1]) 
				filter[level++]= THIRDKEY((*(coccslist[j+1]))[fkey]);
			fkey = fkey1;
		}

		loccs.push_back(new occur(filter, gi));
	}
}
#endif


#ifdef GSPANFLAG
void pattern :: findAllIsomorphisms(){
	GRAPHS ::iterator ip = graphs->begin();
	
	int sum=0;

	for(; ip!= graphs->end(); ip++){
		//cout << *ip << endl;
		AdjMatrix * tg = gb->graph(*ip);
		sum += tg->findAllIsomophisms(M);
	}

	if( sum != occs->size()){
		cout << "sum:" << sum << "\t" << occs->size();
		M->print();
		error("couting error");
	}
}
#endif

#ifdef MAXIMALPATTERN
void pattern :: fillMoreEdge(vector<GNSIZE> & nodes, AdjMatrix * gbi )
{
	int key;
	short n1, n2;
	GLTYPE c;
	list<int>::iterator sit = edgeCan.begin(), sit1;
	list<GLTYPE> ::iterator lit = edgeLabels.begin(), lit1;
	bool flag;
	int n = nodes.size();

	//cout << "new round of updating" << endl;
	for( sit; sit != edgeCan.end(); ){
		key = *sit;
		n1 = key >> 16;
		n2 = (short) (key & 0x0000ffff);
		c = gbi->getLabel(nodes[n-n1-1], nodes[n-n2-1]);
		
		//cout << "n1 n2 c " << n1 << " " << n2  << " " << c << " " << (*lit) << endl;

		flag = false;
		if( c != EMPTY_EDGE ){
			if( (*lit) == EMPTY_EDGE )    (*lit) = c;  //initialization
			else if( c != (*lit))	flag = true;
		}
		else flag = true;
		
		sit1 = sit++; lit1 = lit++; 
		if( flag ){
			edgeCan.erase(sit1);
			edgeLabels.erase(lit1);
		}
	}
}

bool pattern :: isInnerMaximal(){
	return ( !edgeLabels.size() );
}
void pattern :: initialEdgeCan(){
	int key;
	int index = mycost.index();
	int s1 = par->getSize();
	if( s1 == size) index -= s1;
	index = size-2-index;
	//M->print();

	//cout <<"initialization" << endl;
	for(short i=1; i< (short) size; i++){
		for( short j=0; j< i; j++){
			if( M->getLabel(i, j) == EMPTY_EDGE &&  ( i != (size-1) || j< index) ){	
				key = (((int) i) << 16 ) | ((int) j);
				edgeCan.push_back(key);
				edgeLabels.push_back(EMPTY_EDGE);
			}
		}
	}
}

bool pattern :: isMaximal(int threshold){

	CGOCC localCounter;

	/*
	//obtain each occurrence
	int edgeSize = M->esize(), j;
	pattern * currp = this;
	vector<int> gindices(edgeSize, INULLINDEX);
	vector<short> checkpoint(edgeSize+1, 0);
	vector<COCCS*> coccslist;

	int fkey, skey, fkey1, gi, level, i;
	for( i=0; i<= edgeSize; i++){
		if( THIRDKEY((*(currp->coccs))[0]) != GNULLINDEX )
			checkpoint[i] = 1;
		coccslist.push_back(currp->coccs);
		currp = currp->par;
	}
	vector<GNSIZE> filter(size, 0);
	COCCS :: iterator ip = coccs->begin();

	for( i=0; ip != coccs->end(); i++ , ip++){

		fkey = FIRSTKEY((*coccs)[i]), skey = THIRDKEY((*coccs)[i]);
		gi = SECONDKEY((*coccs)[i]);
		//u = skey;
		level = size-1;

		if( checkpoint[0] )	filter[level--]= skey;

		for( j=0; j< edgeSize && fkey != gindices[j]; j++){

			gindices[j] = fkey;
			fkey1 = FIRSTKEY((*(coccslist[j+1]))[fkey]);
			if( checkpoint[j+1]) 
				filter[level--]= THIRDKEY((*(coccslist[j+1]))[fkey]);
			fkey = fkey1;
		}
		*/

	SCANOCCURRENCES
		//update the counting procedure
		if( countSuperGraph(filter, gi, localCounter, threshold) ) return false;
		//loccs.push_back(new occur(filter, gi));
	}

	//is maximal 
	return true;
}

bool pattern :: countSuperGraph(vector<GNSIZE> & nodes, int gi, CGOCC & counter, int threshold){

	AdjMatrix * gbi = gb->graph(gi);

	//inner checking
	GLTYPE c;
	int n = nodes.size()-1;
#ifdef SAFEMODE
	if( n >= 127 ) error("pattern size is too large (> byte)");
#endif

	int key, value, i, j, k;
	vector<GNSIZE> * nei;
	vector<GNSIZE> :: iterator nip;

	for( i=1; i<= n; i++){
		nei = gbi->getNeighbors(nodes[i]);
		for( k=0; k< nei->size(); k++){
			if( (j = (int) (find(nodes.begin(), nodes.end(), (*nei)[k]) - nodes.begin() ) ) < i &&
				M->getLabel(i,j) == EMPTY_EDGE){
				c = gbi->getLabel(nodes[i], nodes[j]);
				key = (i << 24) | (j << 16) | 0 | c;
				if( !counter.count(key) ){
					counter[key] = ( (gi<<16) | 1);
					//if( 1 >= threshold) return true;
				}
				value = counter[key];
				if( (value >> 16) != gi) {
					counter[key] =  (gi<<16) | ( (++value) & 0x0000ffff ); 
					if( (value & 0x0000ffff)  >= threshold ) return true;
				}
			}
		}
	}

	//outer checking
	int x;
	for( i=0; i<= n; i++){
		 nei = gbi->getNeighbors(nodes[i]);
		
		//propose child
		for( j=0; j< nei->size(); j++){
			x = (*nei)[j];
			if( (find(nodes.begin(), nodes.end(), x) == nodes.end() ) && 
					(c=gbi->getLabel(nodes[i], x)) != EMPTY_EDGE  ){
				key = 	((n+1) << 24) | (i<<16) | (gbi->getLabel(x,x) << 8 ) | c;
				if( !counter.count(key) ){
					counter[key] = ( (gi<<16) | 1);
					//if( 1 >= threshold) return true;
				}
				value = counter[key];
				if( value >> 16 != gi){
					counter[key] =  (gi<<16) | (( ++value ) & 0x0000ffff );
					if( (value & 0x0000ffff) >= threshold ) return true;
				}
			}
		}
	}

	return false;
}
#endif

void pattern:: frequentElement(vector<pattern*> & results, int freq){
	
	//getElements
	FREELE candidates;  
	scanElements(candidates,  freq);

	//add to true patterns
	enuFreqEle(candidates, results,  freq);
}

void pattern:: enuFreqEle(FREELE & candidates, vector<pattern*> & results , int freq)
{

}
	
void pattern:: scanElements(FREELE & candidates , int freq)
{
	ELEFRENQC  counter;

	SCANOCCURRENCES
		addInstances(counter, filter, gi, i);
	}
	
    //post processing
    ELEFRENQC ::iterator eip = counter.begin();
	IIDTYPE * tids;
	int sup;

	for( ; eip != counter.end(); eip++){
		tids = eip->second;
		tids->push_back(eip->first);
		if( (sup = getSupport(tids) ) > freq )
		{
			tids->push_back(sup);
			candidates.push_back(tids);
		}
		else
			delete tids;
	}

	//potentially sorting 
}
	
void pattern :: addInstances(ELEFRENQC & counter, vector<GNSIZE> & nodes, int gi, int instanceID)
{
	AdjMatrix * gbi = gb->graph(gi);
	int i, k,j, key;
	GLTYPE c;
	vector<GNSIZE> * nei;
	int n = nodes.size()-1;
	IIDTYPE *  value;

	for( i=1; i<= n; i++){
		nei = gbi->getNeighbors(nodes[i]);
		for( k=0; k< nei->size(); k++){
			if( (j = (int) (find(nodes.begin(), nodes.end(), (*nei)[k]) - nodes.begin() ) ) < i &&
				M->getLabel(i,j) == EMPTY_EDGE){
				c = gbi->getLabel(nodes[i], nodes[j]);
				key = (i << 24) | (j << 16) | 0 | c;
				if( !counter.count(key) )  //initialization
					counter[key] = new IIDTYPE( (coccs->size()+31)/32 );  //initialize with zero
				value = counter[key];
				TURNON(*value, instanceID);
			}
		}
	}
}

/*
bool pattern :: intersectInstances(FREELE & ele1, FREELE & ele2){


	return true //is closed
}
*/

int  pattern:: getSupport( IIDTYPE * iid)
{
	int sup = 0;
	int gi=-1;
	for(int i=0; i< coccs->size(); i++){
		if( GETTID(*iid, i) ){
			if( gi != SECONDKEY((*coccs)[i]) ){
				gi = SECONDKEY((*coccs)[i]);
				sup++;
			}
		}
	}
	return sup;
}
