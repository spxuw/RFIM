#ifndef _PARSER_H_
#define _PARSER_H_

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <ctime>
#include <limits.h>
#include <iterator>
#include <sstream>

#include <list>
#include <vector> 
#include <queue> 
#include <stack> 
#include <map> 
#include <algorithm>

using namespace std;

typedef list<int> Nbl;


void Parser_s2n(char* fname1);
void Parser_s2n(char* fname1, map<string,int>& MAP);
void Parser_s2n_TF(char* fname1);
void Parser_s2n_Weighted(char* fname1);
void Parser_n2n_Undirected_Weighted(char* fname1);

void Parser_gml(char* fname1);
void Parser_paj(char* fname1);
void Parser_t2b(char* fname1);
void Parser_S2A(char* fname1);


void Parser_BiGG(char* fname, bool option_compartment);
void Parser_BiGG_includeBiomass(char* fname, bool option_compartment);
void Parser_BiGG_wo_cofactors(char* fname, bool option_compartment);

void Parser_BiGG_perturbation(string path, string file, bool option_compartment, double p, int seed);


void Parser_BiGG_to_StandardReactionFile(char* fname);
void Parser_BiGG_to_StandardReactionFile_includeBiomass(char* fname);
void Parser_SRF_to_SERF(char* fname);


void Parser_BiGG_to_JMassBalance(char* fname);
void Parser_JMassBalance(char* fname);

void Parser_R2E(char* fname1);


void readsubstrate(string& String, list<string>& substratelist, bool option_compartment, string compartment);

bool GetStoicheometricCoefficient(string reaction, int M, map<string,int>& MAP, vector<vector<double> >& S);



void Parser_R2ODE(char* fname1, map<string,int>& MAP, vector<int>& Receiver);
void Parser_R2ODE(char* fname1, map<string,int>& MAP, vector<string>& Xdot, int ParameterizeType, vector<string>& Parameters);

void Parser_R2ODE_ExtendedSystem(char* fname, map<string,int>& MAP, vector<string>& Xdot, vector<Nbl>& Aoutlist, vector<vector<double> >& S,   vector<string>& NodeName);

void Parser_s2n_FC(char* fname1);

#endif
