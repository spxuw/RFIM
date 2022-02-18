#if !defined (LINK_H)
#define LINK_H

#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <vector>

using namespace std;

class Link
{
 public:
  Link() { head = tail = color = 0;}
  Link(int a, int b) {head = a; tail = b; color = 1;}
  Link(int a, int b, bool c) {head = a; tail = b; color = c;}
  Link(int a, int b, int i, double w) {head = a; tail = b; index = i; weight = w;}

  Link(int a, int b, int i) {head = a; tail = b; index = i;}


  int  GetHead() const {return head;}
  int  GetTail() const {return tail;}
  int  GetIndex() const {return index;}
  double GetWeight() const {return weight;}

  bool GetColor() const {return color;}

  void ChangeHead(int h) {head = h;}
  void ChangeTail(int t) {tail = t;}
  void SwitchColor() {color = !color;}

  void Flip() {swap(head, tail);}

  void AddWeight(double w0) {weight += w0;}


  bool operator==(const Link& other) const {
    return (head==other.GetHead() && tail==other.GetTail());
  }


  friend ostream& operator<<(ostream& output, const Link& e) {
    output << "(" <<  e.GetHead() << "," << e.GetTail() <<")";
    return output;  
  }

 private:
  int head; 
  int tail;   
  int index;  
  
  double weight; 

  bool color;  
  
};


typedef vector<Link>::iterator vLinkitr;


#endif
