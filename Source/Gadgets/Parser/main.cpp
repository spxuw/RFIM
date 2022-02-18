#include <iostream>
#include <fstream>
#include "Parser.h" 
#include "timer.h"

using namespace std;



int main(int argc, char** argv)
{
  if(argc!=3) {
    cout << "\nNumber of arguments is incorrect.\n";
    cout << "\nCommand format: ./parser_s2n path file \n";
    exit(0);
  }

  string  path           = argv[1];
  string  file           = argv[2];

  timer time;
  time.start();

  char elistfile[256];
  sprintf(elistfile, "%s/%s", path.c_str(), file.c_str());
  
  Parser_s2n(elistfile); 
  
  cout << "In total it takes " << time << " s.\n\n";

  exit(1);
}


