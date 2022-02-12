/********************************************************************
 * gBase.cpp
 * Implementation of the class gBase
 *
 ********************************************************************/
#pragma warning (disable:4786 4018 4267)
#include <iostream>
#include <fstream>
#include <assert.h>
#include <vector>
#include <set>
#include <map>
#include <algorithm>
#ifdef _WIN32
	#include <strstream>
#else
	#include <sstream>
#endif

using namespace std;   

#include "gBase.h"
#include "occur.h"
#include "common.h"

void gBase::setAdditionalInfo(string &f)
{
	ifstream fadd(f.c_str());
	char buff[100001];
    string t1, t2;
    int sup, bf;
	int size, density;
	if( !fadd ) error ("can not open additional information\n");

	//istrstream s(buff);
	int i=0;
	while ( fadd.getline(buff, 1000000) ){
#ifdef _WIN32
	    istrstream s(buff);
#else
		istringstream s(buff);
#endif
	    //cout << buff << endl;
        s >> t1  >> sup >> t1 >> bf >> t1 >> size >> t1 >> density >> t1;
	    //cout << "sup " << sup << endl;
        AdjMatrix * M = graph(i++);
        M->setTotalSup(sup);	    
	}
}

void gBase::scan() 
{
	ifstream fnode(nodefile.c_str());
	ifstream fedge(edgefile.c_str());

	if (!fnode)
		error("cannot open input file", nodefile.c_str());
	if (!fedge)
		error("cannot open input file", edgefile.c_str());

	vector< vector<char> *> nodes_graphs;
	vector< vector<Edge *>* > edges_graphs;
	map<char, set<int>*> nodes_freq;

	// scan I
        scanNode(fnode, nodes_graphs, nodes_freq);	
        //cout << "complete node\n";
	scanEdge(fedge, nodes_graphs, edges_graphs);
	//cout << "complete node/edge\n";

	//scanII
	scanII(nodes_graphs, edges_graphs, nodes_freq);
         //cout << "complete scan 2\n";

	scanIII();
	//cout << "complete scan 3\n";

#ifdef DEBUG
	cout << "\n****** Scan II - edge-occur map ******" <<endl;
	cout << "threshold: "<< threshold << endl;
	map<int, vector<occur *>* >::iterator it = edges_occur.begin();
	for (; it != edges_occur.end(); ++it) {
		cout << "key: ";
		printKey(it->first);
		cout << endl;
		for (int k = 0; k < it->second->size(); k++) {
			cout << "occur "<< k << ":  ";
			(*it->second)[k]->print();
		}
	}
	cout << "\n****** ScanII - gb ******" << endl;
	print();
#endif

	// clear dynamic allocated memory
	int i, j;
	for (i = 0; i < nodes_graphs.size(); i++) 
		delete nodes_graphs[i];
	for (i = 0; i < edges_graphs.size(); i++) {
		vector <Edge *> *e_vec = edges_graphs[i];
		for (j = 0; j < e_vec->size(); j++)
			delete (*e_vec)[j];
		delete edges_graphs[i];
	}
	map<char, set<int>* >::iterator ip = nodes_freq.begin();
	for (; ip != nodes_freq.end(); ++ip) 
		delete ip->second;
}

    
void gBase::scanNode(ifstream &fnode, vector< vector<char> *> &nodes_graphs, 
					 map<char, set<int>*> &nodes_freq)
{
	int current_graph = -1;
	while (true) {
		string snode;
		int graph_no, node_no, label_no;

		if (!(fnode >> snode)) // EOF
			break;
			
		if (snode.compare(NODE_TOKEN)){
	          cout << "snode " << snode << endl;
		  error("wrong format in node files");
		}
		fnode >> graph_no >> node_no >> label_no;
		char label = label_no + '0';

		// add node to vector nodes_graphs
		if (graph_no > current_graph) {
			vector<char> *g = new vector<char>;
			nodes_graphs.push_back(g);
			current_graph = graph_no;
		}
		nodes_graphs[graph_no]->push_back(label);

		// add node and graph_no to map nodes_freq
		if (!nodes_freq.count(label))
			nodes_freq[label] = new set<int>;
		nodes_freq[label]->insert(graph_no);
	}

	//threshold = (int)((current_graph + 0.1) * freq)
		;
#ifdef DEBUG
	cout << "\n****** ScanI - Node ******" << endl; 
	int i, j;
	for(i = 0; i < nodes_graphs.size(); i++) {
		vector<char> * v = nodes_graphs[i];
		cout << "graph " << i << endl;
		for(j = 0; j < v->size(); j++) 
			cout << (*v)[j] << " ";
		cout << endl;
	}

	map<char, set<int>* >::iterator iter = nodes_freq.begin();
	for (; iter != nodes_freq.end(); ++iter) {
		cout << iter->first <<": ";
		set<int>::iterator p = (iter->second)->begin();
		for (; p != (iter->second)->end(); ++p)
			cout << *p << " ";
		cout << endl;
	}
#endif
}
	

void gBase::scanEdge(ifstream &fedge, vector< vector<char> *> &nodes_graphs, 
					 vector< vector<Edge *>* > &edges_graphs)
{
	int current_graph = -1;
	while (true) {
		string sedge;
		int graph_no, node1, node2, elabel_no;

		if (!(fedge >> sedge)) // EOF
			break;
			
		if (sedge.compare(EDGE_TOKEN))
			error("wrong format in edge file");

		fedge >> graph_no >> node1 >> node2 >> elabel_no;
		char elabel = elabel_no + '0';

		// check if encounter a new graph
		if (graph_no > current_graph) {
			vector<Edge *> *g = new vector<Edge *>;
			edges_graphs.push_back(g);
			current_graph = graph_no;
		}

		// node1's label always <= node2's label;
		int label1 = (*nodes_graphs[graph_no])[node1];
		int label2 = (*nodes_graphs[graph_no])[node2];
		if (label1 > label2 ) {
			swap(&node1, &node2);
			swap(&label1, &label2);
		}
		
		// construct a new edge and add it to vector edges_graphs
		Edge * e = new Edge(node1, node2, elabel, label1, label2);
		edges_graphs[graph_no]->push_back(e);

		// add edge and graph_no to map edges_freq
		if (!edges_freq.count(e->key)) 
			edges_freq[e->key] = new set<int>;
		edges_freq[e->key]->insert(graph_no);
	}

#ifdef DEBUG
	cout << "\n****** ScanI - Edge ******" << endl; 
	int i, j;
	for(i = 0; i < edges_graphs.size(); i++) {
		vector<Edge *> * v = edges_graphs[i];
		cout << "graph " << i << endl;
		for(j = 0; j < v->size(); j++) 
			printEdge((*v)[j]);
	}

	map<int, set<int>* >::iterator iter = edges_freq.begin();
	for (; iter != edges_freq.end(); ++iter) {
		printKey(iter->first);
		cout <<": ";
		set<int>::iterator p = (iter->second)->begin();
		for (; p != (iter->second)->end(); ++p)
			cout << *p << " ";
		cout << endl;
	}
#endif
}	


// remove infrequent nodes and edges from the input
// construct edges_occur and gb.

void gBase::scanII(vector< vector<char> *> &nodes_graphs, vector< vector<Edge *>* > &edges_graphs,
				map<char, set<int>*> &nodes_freq)
{
	int i, j;

    // scan nodes
	int g_current = 0;
	for (i = 0; i < nodes_graphs.size(); i++) {
		// for each graph, set infreq. nodes to -1 and reindex freq. nodes
		vector<int> *index_vec = new vector<int>;
		v_reindex.push_back(index_vec);
		int n_index = 0;
		bool graph_added = false;

		for (j = 0; j < nodes_graphs[i]->size(); j++) {
			char n_label = (*nodes_graphs[i])[j];
			if (nodes_freq[n_label]->size() >= threshold) {
				index_vec->push_back(n_index++);
	
				// add node to gb
				if (!graph_added) {
					gb.push_back(new AdjMatrix);
					graph_added = true;
				}
				gb[g_current]->addNode(n_label);

				//build node_occur
				if( !node_occur.count(n_label) )		
					node_occur[n_label] = new vector<occur*>();
				   //(node_occur[n_label])->push_back(new occur(new vector<GNSIZE>(1, (n_index-1)), g_current) );
				vector<GNSIZE> v;
				v.push_back(n_index-1);
				(node_occur[n_label])->push_back(new occur(v, g_current) );

			}
			else
				index_vec->push_back(-1);
		}

		if (n_index != 0) 
			g_reindex.push_back(g_current++);
		else
			g_reindex.push_back(-1);
	}

	//scan edges
	for (i = 0; i < edges_graphs.size(); i++) {
		vector<Edge *> *eg = edges_graphs[i];
		int gi = g_reindex[i];
		for (j = 0; j < eg->size(); j++) {
			Edge *e = (*eg)[j];
			if (edges_freq[e->key]->size() >= threshold) {

				// update edge's nodes to new index
				int ne1 = (*v_reindex[i])[e->node1_no];
				int ne2 = (*v_reindex[i])[e->node2_no];

				if( ne1 < 0 || ne2 < 0)	error("node reindex error");

				// add to gb
				if( ne1 == 1033 && ne2 == 41)
				  cout << "gi " << gi << endl;
				gb[gi]->addEdge(ne1, ne2, e->label);

				// add to edges_occur
				if (!edges_occur.count(e->key))
					edges_occur[e->key] = new vector<occur *>;

				vector<GNSIZE> v;
				v.push_back(ne2);
				v.push_back(ne1);
				edges_occur[e->key]->push_back(new occur(v, gi));

				int label1 = (*nodes_graphs[i])[e->node1_no];
				int label2 = (*nodes_graphs[i])[e->node2_no];

				if( label1 == label2 ){
					v.clear();
					v.push_back(ne1);
					v.push_back(ne2);
					edges_occur[e->key]->push_back(new occur(v, gi));
				}

				//update node index in edge
				e->node1_no = ne1;
				e->node2_no = ne2;
				
			}
		}
	}
}


void gBase::scanIII()
{

	//set up the adjacency list representation
	int i;
	for ( i = 0; i < gb.size(); i++) {
		gb[i]->buildAdjList();
	}


	//sort the edge occur and node occur vector
	map<int, vector<occur*>*, greater<int> > :: iterator ip;	
	for(ip = node_occur.begin(); ip!= node_occur.end(); ip++){
		//vector<occur*> * t = ip->second;
		//cout << "sort node " << ip->first << endl;
		sort(ip->second->begin(), ip->second->end(), Comp());
	}
	for(ip= edges_occur.begin(); ip!= edges_occur.end(); ip++){
		//vector<occur*> * t = ip->second;
		sort(ip->second->begin(), ip->second->end(), Comp());
	}

	//scan to make sure no redundancy
	for(ip = node_occur.begin(); ip!= node_occur.end(); ip++){
		//vector<occur*> * t = ip->second;
		vector<occur*> * t= ip->second;
		//cout << "node " << hex << ip->first << dec << endl;

		for(i=0; i<t->size()-1; i++){
			occur * oc1 = (*t)[i], * oc2 = (*t)[i+1];
			if( oc1->getGraph() == oc2->getGraph() && oc1->getNodes(0) == oc2->getNodes(0) )
				error("redudant node ocurrence");
		}
	}

	for(ip = edges_occur.begin(); ip!= edges_occur.end(); ip++){
		//vector<occur*> * t = ip->second;
		vector<occur*> * t= ip->second;
		//cout << "node " << hex << ip->first << dec << endl;

		for(i=0; i<t->size()-1; i++){
			occur * oc1 = (*t)[i], * oc2 = (*t)[i+1];
			if( oc1->getGraph() == oc2->getGraph() && oc1->getNodes(0) == oc2->getNodes(0)
				&& oc1->getNodes(1) == oc2->getNodes(1) ) {
                                cout << "graph " << oc1->getGraph() << " node 1: " << oc1->getNodes(0) <<
                                       " node 2: " << oc1->getNodes(1) << endl;
				error("redudant edge ocurrence");
                        }
		}
	}

}


void gBase::swap(int *x, int *y) {
	int tmp = *x;
	*x = *y;
	*y = tmp;
}


void gBase::print()
{
	for (int i = 0; i < gb.size(); i++) {
		cout << "Graph " << i << endl;
		gb[i]->print();
	}
}

void gBase::printEdge(Edge *e)
{
	cout << "node1: " << e->node1_no <<" node2: " << e->node2_no
		<< " edge label: " << e->label << " key: ";
	printKey(e->key);
	cout << endl;
}

void gBase::printKey(int k)
{
	cout << (char)(k>>16) << "-" << (char) ((k>>8) & 0xff0000ff) 
		<< "-" << (char)(k & 0xff0000ff);
}

void gBase::removeEdge(int key)
{
	if( edges_occur.count(key) ){
		vector<occur*> * v = edges_occur[key]; 

		for( int i =0; i<v->size(); i++){
			occur * oc = (*v)[i];
			int gi = oc->getGraph();
			int n1 = oc->getNodes(0);
			int n2 = oc->getNodes(1);
#ifndef VALIDATION
			gb[gi]->removeEdge(n1, n2);
#endif
		}
	}
	else{
		error("no such key");
	}
	//

}

gBase::~gBase() 
{
	// clear dynamic allocated memory
	int i;
	for (i = 0; i < v_reindex.size(); i++) 
		delete v_reindex[i];
	
	map<int, set<int>* >::iterator ip = edges_freq.begin();
	for (; ip != edges_freq.end(); ++ip) 
		delete ip->second;

	map<int, vector<occur *>*, greater<int> >::iterator it = edges_occur.begin();
	for (; it != edges_occur.end(); ++it) 
		delete it->second;  // occur* deallocated by lh

	for(i=0; i< gb.size(); i++)	
		delete gb[i];
}

bool Comp :: operator()( occur * const  & oc1,  occur * const & oc2)
{
	
//      pattern ** p1 = (pattern**) ip1;^M
//      pattern ** p2 = (pattern**) ip2;^M

        int g1 = oc1->getGraph(), g2 = oc2->getGraph();
		int n11 = oc1->getNodes(0), n21 = oc2->getNodes(0);

		if( g1 != g2 ){
			if( g1 < g2 ) return true;
			else return false;
		}
		else if( n11 != n21 ){
			if( n11 < n21 ) return true;
			else return false;
		}
		else{
			if( oc1->size() < 2 ){
				return true;
			}
			else{
				int n12 = oc1->getNodes(1), n22 = oc2->getNodes(1);
				return (n12 < n22);
			}
		}
     
}

