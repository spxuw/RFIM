#ifndef _COLOR_H_
#define _COLOR_H_

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <string>
#include <list>
#include <set>
#include <stack>
#include <deque>
#include <limits>
#include <vector>



using namespace std;

/*
const char WHITE    = 0;
const char GRAY     = 1;
const char BLACK    = 2;
const char GREEN    = 3;
const char BLUE     = 4;
const char RED      = 5;
const char PURPLE   = 6;
const char GREEN1   = 31;
const char GREEN2   = 32;
const char GREEN3   = 33;
*/

const char Normal[] = { 0x1b, '[', '0', ';', '3', '9', 'm', 0 };
const char Black[] = { 0x1b, '[', '0', ';', '3', '0', 'm', 0 };
const char Red[] = { 0x1b, '[', '0', ';', '3', '1', 'm', 0 };
const char Green[] = { 0x1b, '[', '0', ';', '3', '2', 'm', 0 };
const char Yellow[] = { 0x1b, '[', '0', ';', '3', '3', 'm', 0 };
const char Blue[] = { 0x1b, '[', '0', ';', '3', '4', 'm', 0 };
const char Purple[] = { 0x1b, '[', '0', ';', '3', '5', 'm', 0 };
const char Cyan[] = { 0x1b, '[', '0', ';', '3', '6', 'm', 0 };
const char Lgray[] = { 0x1b, '[', '0', ';', '3', '7', 'm', 0 };
const char Dgray[] = { 0x1b, '[', '0', ';', '3', '8', 'm', 0 };
const char Bred[] = { 0x1b, '[', '1', ';', '3', '1', 'm', 0 };
//for bold colors, just change the 0 after the [ to a 1

#endif /* _COLOR_H_ */
