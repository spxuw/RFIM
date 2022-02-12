#pragma warning (disable:4786 4018 4267)
//#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <map>
#include <ctime>
#include <iomanip>
#include <errno.h>
#include <stdio.h>

using namespace std;

#include "common.h"

#include "myTimer.h"
#include "register.h"
#include "pattChild.h"
#include "induPatt.h"
#include "gBase.h"

void setTimer(MyTimer & t);
void setRegister(Register & myregister);
void statistics(MyTimer & Timers, Register & myregister);
//void populateChild(children * & c, map<int, vector<occur*> *> & cm);
#ifdef VALIDATION
typedef map<int, vector<AdjMatrix*> *> INNER;
typedef map<int, INNER * > OUTER;
void frequency_checking(vector<AdjMatrix*> & r, gBase * gb );
void count(AdjMatrix * m, gBase * gb, set<int> & oc);
void redudancy_checking(vector<AdjMatrix*> & r);
void verify(OUTER & out);
void verify_inner(vector<AdjMatrix *> * v, int &k );
#endif

void validation( vector<AdjMatrix*> & r, gBase * gb );
void result_output(vector<AdjMatrix*> & result, string f1, string f2, string f3,
				   map<int, int> & g_reindex, int d, int s, int sizeup, DLONG ind);
void obtainMaximalPattern(vector<AdjMatrix*> & result);
void obtainExeclusiveMaximalPattern(vector<AdjMatrix*> & result);

int main (int ARGC, char * ARGV[]){
	
	string f1, f2, f3, f4, f5;
	int freq, density, sizelimit, sizeup;

	density   = 0;
	sizelimit = PATTERN_SIZE_LOW;
	sizeup    = PATTERN_SIZE_UP;

	if (ARGC < 7) {
	   cout <<  "FFSM nodeFile edgeFile outNodeFile outEdgeFile" 
		   << " outFeatureFile support density sizelimit sizeuplimit\n";
           //cout << "FFSM nodeFile edgeFile support\n";
       return 1;
	}

#ifdef _WIN32

	f1 = ARGV[1];
	f2 = ARGV[2];
	f3 = ARGV[3];
	f4 = ARGV[4];
	f5 = ARGV[5];
	freq = atoi( (char*)ARGV[6]);
	//optional parameters
	if( ARGC >= 8  ) density   = atoi( ARGV[7]);
	if( ARGC >= 9  ) sizelimit = atoi( ARGV[8] );
	if( ARGC >= 10 ) sizeup    = atoi( ARGV[9] );

	//cout << f1 << " " << f2 << " " << f3 << " " << f4 << endl;
	/*
	f1 = "n1.txt";
    f2 = "e1.txt";
	f3 = f1 + "_out.txt";
	f4 = f2 + "_out.txt";
	f5 = f3 + "_f.txt";
    freq = 1;
	density = 1000;
	sizelimit = 4;
	sizeup    = 15;
	*/
#else
	f1 = ARGV[1];
	f2 = ARGV[2];
	f3 = ARGV[3];
	f4 = ARGV[4];
	f5 = ARGV[5];
	freq = atoi( (char*)ARGV[6]);
	//optional parameters
	if( ARGC >= 8  ) density   = atoi( ARGV[7]);
	if( ARGC >= 9  ) sizelimit = atoi( ARGV[8] );
	if( ARGC >= 10 ) sizeup    = atoi( ARGV[9] );	
#endif

        //clean f3 f4 f5
    ofstream of3(f3.c_str()); of3.close();
	ofstream of4(f4.c_str()); of4.close(); 
	ofstream of5(f5.c_str()); of5.close();  		
    cout << "=========Welcome using FFSM: fast frequent subgraph mining=========" << endl;
	cout << "Developed by University of North Carolina at Chaple Hill" << endl;
 

	//set up timer and registers
	MyTimer Timers;
	setTimer(Timers);
	Register myregister;
	setRegister(myregister);

	gBase * gb = new gBase(f1, f2, freq);
	cout << "graph database size " << gb->size();
	
	vector<AdjMatrix * > result;
	int ts = gb->getThreshold();
	if( ts < 1) ts = 1;
	cout << " threshold " << ts << endl;

	/*for( int i =0; i < gb->size(); i++){
		gb->graph(i)->print();
	}*/


#ifdef INDUCEDS
	iPattern * p = new iPattern(gb, Timers, myregister, cm, nm);
	//gb->deleteMap();
	//gb->deleteNodeMap();
	Timers.start_timer(TIMER1);
	p->fsm1(result, ts, 0 , density, f3, f4, f5, sizelimit, gb->getGIndex());
	
#else
	//pattern * p = new pattern(gb, Timers, myregister, cm, nm);
	//construct the initial pattern
	map<int, vector<occur*> *, greater<int> > cm = gb->getMap();
	map<int, vector<occur*> *, greater<int> > nm = gb->getNodeMap();
	pattern * p = new pattern();
	map<int, int> gMap;
	gb->getGIndex(gMap);

	p->initPattern(gb, Timers, myregister);
	p->setChild( cm,  nm);
	
	Timers.start_timer(TIMER1);
	p->fsm1(result, ts, 0, density /*,  gb->getGIndex()*/, f3, f4, f5, sizelimit, sizeup, gMap);
	
#endif

	Timers.stop_timer(TIMER1);
	cout << "Search Time "  << Timers.value_timer(TIMER1) << " seconds" << endl;

	 //get the maximal patterns
	 Timers.start_timer(TIMER1);

	 int totalTree = 0;
	 for(int i=0; i< result.size(); i++){
		totalTree += result[i]->isTree();
	 }
	 //cout << "total number of trees " << totalTree << endl;

	 //output statistics
	 statistics(Timers, myregister);

	 //obtainMaximalPattern(result);
	 Timers.stop_timer(TIMER1);
	 //cout << "Post processing time "  << Timers.value_timer(TIMER1) << " seconds" << endl;
	 //obtainExeclusiveMaximalPattern(result);

	//for result output
	DLONG ind = p->getTotalPatterns() ;
	result_output(result, f3, f4, f5, gMap, density, sizelimit, sizeup, ind);

	//validation
	//validation(result, gb);
	return 0;

}

void statistics(MyTimer & Timers, Register & myregister)
{
	//timer related
	cout << "statistics: (if TIMERFLAG or REGISTERFLAG is turned on)" << endl;
#ifdef TIMERFLAG
	Timers.print();
#endif

#ifdef REGISTERFLAG
	myregister.print();
#endif
}

void setRegister(Register & myregister)
{
	
	myregister.allocate_register(REGISTER1, "total_tried_children_proposing");
	myregister.allocate_register(REGISTER2, "total_children_instance_proposed");
	myregister.allocate_register(REGISTER3, "total_children");
	myregister.allocate_register(REGISTER4, "passed");
	myregister.allocate_register(REGISTER5, "okay");
	myregister.allocate_register(REGISTER6, "failed");
	myregister.allocate_register(REGISTER7, "no_frequent");
	myregister.allocate_register(REGISTER8, "extendedd");
	myregister.allocate_register(REGISTER9, "ext_total_children");
	myregister.allocate_register(REGISTER10, "ext_passed");
	myregister.allocate_register(REGISTER11, "ext_okay");
	myregister.allocate_register(REGISTER12, "extfailed");
	myregister.allocate_register(REGISTER13, "ext_no_frequent"); 
	
	//myregister.print();
}

void setTimer(MyTimer & Timers){

	//automatic initialize the timer after allocation
	Timers.allocate_timer(TIMER1, "overall time");
	Timers.allocate_timer(TIMER2, "subgraph_testing");
	Timers.allocate_timer(TIMER3, "children_proposing");
	Timers.allocate_timer(TIMER4, "handling_data_structure");
	Timers.allocate_timer(TIMER5, "inner_proposing");
	Timers.allocate_timer(TIMER6, "outer_proposing");
	Timers.allocate_timer(TIMER7, "isomorphism");
	 
}

void validation( vector<AdjMatrix*> & r, gBase * gb)
{

	
#ifdef VALIDATION
    cout << "start frequency check: " << r.size() << endl;
	frequency_checking(r, gb);
    cout << "complted"  << endl;

	return;

	cout << "start redundancy check" << endl;
	redudancy_checking(r);
	cout << "complted"  << endl;

#endif
  

}

#ifdef VALIDATION
void frequency_checking(vector<AdjMatrix*> & r, gBase * gb )
{
	for(int i =0; i< r.size(); i++){
		AdjMatrix * m = r[i];
		set<int> oc;
		count(m, gb, oc);
		set<int> & oc1 = m->getGraphs();

		if( oc.size() != oc1.size() ){
			cout << "matrix size" << m->size() << " size 1 " << oc.size() << " size 2 " << oc1.size() << endl;
                        m->print();			
			log("not match");
                }
		else{
			if( (i+1) % 10 == 0) {
				cout << ".";
				if ( (i+1) % 200 == 0 ) cout << endl;
				cout.flush();
			}
			set<int> :: iterator ip;
			for( ip = oc.begin(); ip != oc.end(); ip++){
				if( !oc1.count(*ip) )
					cout << "not in graph" << endl;
			}
		}

	}
}

void redudancy_checking(vector<AdjMatrix*> & r)
{

	//set up the data structure
	OUTER out;
	int i=0; 
	for( i=0; i< r.size(); i++){
		AdjMatrix * m = r[i];
		vector<AdjMatrix *> * v;
		INNER * inr; 
		int s = m->size(), es = m->esize();

		if( !out.count(s) ){
			v = new vector<AdjMatrix*>();
			v->push_back(m);
			inr = new INNER();
			(*inr)[es] = v;
			out[s] = inr;
		}
		else{
			inr = out[s];
			if( !inr->count(es) ){
				v = new vector<AdjMatrix*>();
				v->push_back(m);
				(*inr)[es] = v;
			}
			else{
				v = (*inr)[es];
				v->push_back(m);
			}
		}
	}

	//verify the structure
	verify(out);
}

void verify(OUTER & out)
{
	typedef map<int, vector<AdjMatrix*> *> INNER;
	typedef map<int, INNER * > OUTER;

	OUTER::iterator ip ;

	int k =0;
	for( ip = out.begin(); ip != out.end(); ip++){
		INNER * inr = ip->second;
		INNER::iterator ip1;

		for( ip1 = inr->begin(); ip1 != inr->end(); ip1++){
			vector<AdjMatrix *> * v = ip1->second;

		
			verify_inner(v, k);
		}
	}
}

void verify_inner(vector<AdjMatrix *> * v, int &k )
{

	for(int i=0; i<v->size(); i++){		
		
		if( ++k % 10 == 0){
			cout <<"." ;
			if( k % 200 == 0) cout << endl;
			cout.flush();
		}

		AdjMatrix * m1 = (*v)[i];
		for(int j=i+1; j<v->size(); j++){
			AdjMatrix * m2 = (*v)[j];
			if( m1->isIsomorphicOf(m2) ){
				m1->print();
				m2->print();
				error("redundancy");
			}
		}
	}
}

void count(AdjMatrix * m, gBase * gb, vector<int> & oc){
	
	for( int i = 0; i< (gb->size()+31)/32*32; i++){
        if( i < gb.size() & m->isSubgraphOf(gb->graph(i), 0) )
			oc.insert(i);
	}
}

#endif

void result_output(vector<AdjMatrix*> & result, string f1, string f2, string f3, 
				   map<int,int> & g_reindex, int d, int s, int sizeup, DLONG ind)
{

	int i;

	//cout << GFILE << PFILE << endl;
	cout <<"start output final result: " << f1 << endl;

	for( i=0; i< result.size(); i++){
		AdjMatrix * m = result[i];

#ifdef SAVERESULTS
		bool flag = true;
#ifdef DENSITY
		if(!m->denseEnough(d, true) ) flag = false;
#endif
		if(m->size() >=s && flag && m->size() <= sizeup)
	   	  m->print(ind++, f1, f2, f3, g_reindex); 

	/*	else if( !flag){
		  cout << " element is not dense enough\n";
		}*/
		
#else
		m->print();
#endif
	}

        //flush results
        ofstream of1(f1.c_str(), ios::app);
        of1.flush();
        ofstream of2(f2.c_str(), ios::app);
        of2.flush();
	cout << "end" << endl;

}

void obtainMaximalPattern(vector<AdjMatrix*> & result)
{
	vector<AdjMatrix*> lb;

	if( result.size() > 0)	lb.push_back(result[0]);

	cout << "total pattern " << result.size() << endl;
	//check the maximal pattern
	for(int i=1; i< result.size(); i++){
		AdjMatrix * t = result[i];

		if( !( (i+1) % 100 ) ) cout << ".";
		if( !( (i+1) % 1000 ) ) cout << i << endl;

		bool flag = true;
		for(int j=0; j< lb.size() && flag; j++){
			//lb[j]->print();
			flag = !(t->isSubgraphOf(lb[j], 0));
#if (defined SAFEMODE && defined VFLIB)
			//cout << "double check " << lf << endl;
			if( flag == lb[j]->isSubgraph(t) ){
			  //cout << "is subgraph " << lf << "\t0: no";
				error("subgraph mismatch");
			}
#endif
		}
		if( flag ){ /*cout << "t is maximal " << endl; */ lb.push_back(t); }
		else{ /* cout << "t is not maximal " << endl;*/ delete t; }
	}
	
    result.clear();
	swap(result, lb);
	cout << "total maximal pattern " << result.size() << endl;
}

void obtainExeclusiveMaximalPattern(vector<AdjMatrix*> & result)
{

	vector<AdjMatrix*> lb;
	cout << "total pattern " << result.size() << endl;
	for(int i=0; i< result.size(); i++){
		AdjMatrix * t = result[i];
		//t->print();
		bool flag = true;
		for(int j=0; j < result.size(); j++){
			if( i != j && t->isSubgraphOf(result[j],0) ){
				flag = false; break;
			}
		}
		if( flag) lb.push_back(t);
	}
	result.clear();
	swap(result, lb);
	cout << "total maximal pattern " << result.size() << endl;
}

