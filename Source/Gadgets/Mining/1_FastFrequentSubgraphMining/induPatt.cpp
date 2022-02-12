#pragma warning (disable:4786 4018 4267)
#include <vector>
#include <fstream>
#include <algorithm>
#include <set>
using namespace std;

#include "induPatt.h"

int iPattern:: isCarnonicalForm(void) { 

	stage = 1;
	
#ifdef TIMERFLAG
	Timers.start_timer(TIMER4);
#endif

	fixM(); 

#ifdef TIMERFLAG
	Timers.stop_timer(TIMER4);
#endif
	

#ifdef TIMERFLAG
	Timers.start_timer(TIMER2);
#endif

	int flag = M->isCanonicalForm();  
	
#ifdef TIMERFLAG
	Timers.stop_timer(TIMER2);
#endif

	return flag;

}

void iPattern:: sprint(int level)
{
	cout << "the adj matrix is" << endl;
	if( M )
		M->print();
	else
		cout <<"matrix is not initialized" <<endl;
	
	//cout << "pattern optimal "  << optimal << " cost " << hex << mycost.intValue() << dec << endl;
	if( level & 2  ){
		cout << "the coccs structures are: size: " << coccs->size() << endl;
		
		COCCS :: iterator ip1;
		for(ip1 = coccs->begin(); ip1 != coccs->end(); ip1++){
			DLONG key = *ip1;
			cout << "index is " << (int)(int)(key >> HALFDLONG ) << " node is " << (SECONDKEY(key)) << endl;
		}
	}

	if( level & 4  ){
		cout << "there are total: " << child.size() << " children in this pattern " << endl;
		
		//vector<pattern *> :: iterator ip1 ;
		//for(ip1 = child.begin(); ip1 != child.end(); ip1++){
		//	(*ip1)->sprint(5);
		//}
	}
}

void iPattern:: print()
{
	cout << "the adj matrix is" << endl;
	if( M )
		M->print();
	else
		cout <<"matrix is not initialized" <<endl;

	cout << "the index structure are" << endl;
	COCCS :: iterator ip1;

	for(ip1 = coccs->begin(); ip1 != coccs->end(); ip1++){
		DLONG key = *ip1;
		cout << "index is " << (int) (key >> HALFDLONG ) << " node is " << (int) (SECONDKEY(key)) << endl;
	}

	//cout << "cost is " << hex << mycost.intValue() << endl;
}

void iPattern:: fsm1( vector<AdjMatrix * > & result, int freq, int level){

	iCHLTYPE :: iterator ip;
	cout << "start search total frequent edge :" << child.size() << endl;

	vector<iPattern*> wlist;
	int totalnum = 0;

	for( ip = child.begin(); ip!=child.end(); ip++){
		 //cout << "search with support: " << (*ip)->sup() << endl;
#ifdef DEBUG3
		cout << "start search" << endl;
		(*ip)->sprint(3);
#endif

		(*ip)->fsm2(result, freq, level+1, totalnum);

	}

	cout << "Total frequent patterns: " << dec << totalnum << "\t";

#ifdef DEBUG8
	edbLocal->print();
#endif
	
}


void iPattern:: fsm2( vector<AdjMatrix * > & result, int freq, int level, int& totalnum){

	stage = 4;
	if( level >= 2){
	    totalnum++;
	    if( totalnum % 10000 == 0) cout << "." ;
	    if( totalnum % 100000 == 0) cout << endl;
#ifdef SAVERESULTS
		if( M )
			//M->setGraphs(occs);
			//M->print(totalnum-1, PFILE, GFILE, g_reindex);	 
			result.push_back(M);
#endif
	}

	iPattern * pi, *pj;

	for( int i = 0; i< child.size(); i++){
		pi = child[i];

#ifdef DEBUG3
	    for( int lii = 0; lii< level; lii++)
			cout << "  ";
		cout << "cost and child size " << hex << mycost.getCost() << "  size:" <<dec << child.size() << endl;
#endif

		if( pi->optimal ){  //the pattern must be in carnonical form for subsequent search
		
			for(int j = i; j< child.size(); j++){
				pj = child[j];
				pi->join(pj, freq);
			}

#ifdef DEBUG2
			int k;
			cout << "children are: " << endl;
			for( k=0; k< pi->child.size(); k++){
				(pi->child[k])->sprint();
			}			
#endif


#ifdef TIMERFLAG
			Timers.start_timer(TIMER3);
#endif

#ifdef TIMERFLAG
			Timers.start_timer(TIMER6);
#endif
			//propose outer child
			pi->proposeOuterChild(freq);
		
#ifdef TIMERFLAG
			Timers.stop_timer(TIMER6);
#endif

#ifdef TIMERFLAG
			Timers.stop_timer(TIMER3);
#endif

		
#ifdef DEBUG2
			cout << "completed child proposing" << endl;
			cout << "children are (proposed by outer): " << endl;
			for( ; k< pi->child.size(); k++){
				(pi->child[k])->sprint();
			}
#endif
			pi->fsm2( result,freq, level+1 , totalnum);

			//if( pi->size == 2) gb->removeEdge(pi->mycost.intValue());
		}
		delete pi; /*definitely need to be taken care of*/
	}
}

/*
void iPattern:: fixOccs(void){
	for(int i=0; i< occs->size(); i++){
		occur * occ = (*occs)[i];
		//occ->fixOccs();
	}
}
*/
	
void iPattern:: constructM (int mc){
	M = new AdjMatrix();

	if( size == 2){
		M->addNode(COSTINDEX(mc));
		M->addNode(COSTNL(mc));
		M->addEdge(0, 1, COSTEL(mc));
	}
	else
		M->addNode(COSTNL(mc));
}

void iPattern:: fixM(void){

	M = new AdjMatrix(M, mycost);
	stage &= 1;
}	

void iPattern::setChild(EDGES & cm, EDGES & nm)
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

iPattern * iPattern :: setUpChild(vector<int> u, EDGES & edges, EDGES & gnodes, int nlabel)
{
	iPattern * p1 = new iPattern(this, nlabel, 1), * p2;
	iCHLTYPE & c1 = p1->child;

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
		p2 = new iPattern(p1, u[i], 2);
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

void iPattern :: setUpEdges(iPattern * p2, OCCS * es, TRA & tra1)
{
	//set up for the p2(edge)
	COCCS * cos2 = p2->coccs;	GRAPHS1 * gs = p2->gids;

	int gi, d1, d2, t1,j;
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
		if( !gs->size() || (*gs)[gs->size()-1] != gi) gs->push_back(gi);
		
	}
}

void iPattern :: join(iPattern * p2, int freq){
	
	int sum=0;
	GRAPHS1 ::iterator iip1 = gids->begin(), iip2 = (p2->gids)->begin();
    while( iip1 != gids->end() && iip2 != p2->gids->end() ){
		if( *iip1 == * iip2){ iip1++; iip2++; sum++;}
		else if( * iip1 > *iip2) *iip2++;
		else iip1++;
	}

	//no way to be frequent
	if( sum < freq) return ;

	//get the type of the join
	propose_outer2(p2, freq);

}


//for outer + outer case 2
void iPattern :: propose_outer2(iPattern * p2, int freq)
{
	COCCS :: iterator ip1, ip2, ips, ipe;
	COCCS * cos2 = p2->coccs;
	ip1 = coccs->begin(); ip2 = cos2->begin();
	CANDIDATE  cand;

	GRAPHS1 * gid3 ;
	int index =0, t, k;
	
	while( ip1 != coccs->end()  && ip2 != cos2->end() ){

		t = FIRSTKEY(*ip1);
		k = FIRSTKEY(*ip2);

		if( t == k ){

			//find the ipe
			ips = ipe = ip2;
			while( (ipe != cos2->end()) &&  ( FIRSTKEY(*ipe) == t) ) ipe++;

			bool flag = false;
			for( ; ip2 != ipe; ip2++){
				if( *ip1 !=  *ip2){
					int k1 = THIRDKEY(*ip1), k2 = THIRDKEY(*ip2);
					int gi = SECONDKEY(*ip1);

					int kc = MAKECOST(0, gb->graph(gi)->getLabel(k1, k2), gb->graph(gi)->getLabel(k2, k2));
					
					//cost need to be updated
					iPattern * p3;  

					if( !cand.count(kc) ){
						Cost lc = p2->mycost;
						lc.insert(size-1, gb->graph(gi)->getLabel(k1, k2)); 
						cout << "lc is " << lc.getCost() << endl;
						p3 = new iPattern(this, lc, size+1);
						cand[kc] = p3;
					}
					else
						p3 = cand[kc];

					p3->coccs->push_back(MAKEKEY(index, gi, THIRDKEY(*ip2))); /*need work*/
					gid3 = p3->gids;
					if( !gid3->size() || gi != (*gid3)[gid3->size()-1]) gid3->push_back(gi);
				}
			}
			
			ip1++;	index++; ip2 = ips;
			
		}
		else if( t > k ){
			ip2++; 
		}
		else{
			ip1++; index++;
		}
	
	}

	scan(cand, freq, child);
}

void iPattern:: proposeOuterChild(int freq){

		 	 
	COCCS :: iterator ip = coccs->begin();
	int i, u, gi, x, key, level;
	CANDIDATE cand;
	int elabel, nlabel;
	AdjMatrix * gbi;
	iPattern * p3;

    //for update
	int j;
	iPattern * currp = this;
	vector<COCCS*> coccslist;
	vector<int> gindices(size, INULLINDEX);
	vector<short> checkpoint(size+1, 0);

//	DLONG mykey;
	int fkey, skey, fkey1;

	for( i=0; i<= size; i++){
		coccslist.push_back(currp->coccs);
		currp = currp->par;
	}

	vector<GNSIZE> filter(size, 0);
	
	for( i=0; ip != coccs->end(); i++ , ip++){

		//update
		 
		//vector<GNSIZE> filter;
		//fkey = i, skey = THIRDKEY((*coccs)[i]);
		fkey = FIRSTKEY((*coccs)[i]), skey = THIRDKEY((*coccs)[i]);
		gi = SECONDKEY((*coccs)[i]);
		u = skey;
		level = 1;

		//update
		filter[0]= skey;
		//end

		for( j=0; j< size && fkey != gindices[j]; j++){

			gindices[j] = fkey;
			fkey1 = FIRSTKEY((*(coccslist[j+1]))[fkey]);
			//if( skey != GNULLINDEX)
				//filter.push_back(skey);
			if( checkpoint[j+1]) 
				//filter[level++]= skey;
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
		
		//for debugging only
		vector<GNSIZE> * nei = gbi->getNeighbors(u);
		
		//propose child
		for( j=0; j< nei->size(); j++){
			x = (*nei)[j];

			bool rflag = true;

			for( int m =0; m< filter.size(); m++){
				if( filter[m] == x){
					rflag = false;
					break;
				}
			}

			if( rflag && gbi->getLabel(u, x) != EMPTY_EDGE){
			//if( !filter1.count(x) ) {

			elabel = gbi->getLabel(u, x);
			nlabel = gbi->getLabel(x, x);
			key = elabel << 8 | nlabel;

			
#ifdef REGISTERFLAG 
			Registers.peg(REGISTER3);
#endif
			if( cand.count(key) ){
				p3 = cand[key];
			}
			else{
				p3 = new iPattern(this, key, size+1);
				cand[key] = p3;
			}
		
			//update
			if( !p3->gids->size() || (gi != (*(p3->gids))[(p3->gids)->size()-1]) )
				p3->gids->push_back(gi);
			//end

			p3->coccs->push_back(MAKEKEY(i, gi, x));
			}
		}
	 }
	//scan children to see whether it is frequent and optimal 
	scan(cand,freq, child);
}

void iPattern:: scan(CANDIDATE & cand, int f, vector<iPattern*> & result)
{

	map<int, iPattern*, greater<int> > :: iterator ip;
	
	for(ip = cand.begin(); ip!= cand.end(); ip++){
		iPattern *  p = ip->second;

#ifdef REGISTERFLAG 
		Registers.peg(REGISTER8);
#endif
		if( p->sup() < f ) { 

#ifdef REGISTERFLAG 
			Registers.peg(REGISTER7);
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
				Registers.peg(REGISTER6);
#endif
				delete p;	 /*definitely need to be taken care of*/
			}
			else {
				p->optimal = t;
				/*if( t == 1)	p->fixOccs();*/
				result.push_back(p);

#ifdef REGISTERFLAG 
				if( !t )
					Registers.peg(REGISTER5);
				else			
					Registers.peg(REGISTER4);
#endif
			}
		}

	}
}

#ifdef GSPANFLAG
void iPattern :: findAllIsomorphisms(){
	GRAPHS1 ::iterator ip = gids->begin();
	
	int sum=0;
	//cout << "graph population size " << graphs->size() << endl;
	//M->print();

	for(; ip!= gids->end(); ip++){
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


