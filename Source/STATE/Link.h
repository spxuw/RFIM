#ifndef LINK_H_
#define LINK_H_
#include <iostream>
using namespace std;


//////////////////////////////////////////////////////////////////////////////////////////////////////////
typedef list<int> Nbl; // neighbor list
typedef Nbl::iterator Nbl_itr;

const string numbers="0123456789";

class Link
{
 public:
  Link() { head = tail = color = 0;}
  Link(int a, int b) {head = a; tail = b; color = 1;}
  Link(int a, int b, bool c) {head = a; tail = b; color = c;}
  Link(int a, int b, int i, double w) {head = a; tail = b; index = i; weight = w;}
  Link(int a, int b, int i) {head = a; tail = b; index = i;}

  //Link(const Link&);                                         // copy constructor  
  //Link& operator=(const Link&);                              //copy assignment    

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
  void SetWeight(double w0) {weight = w0;}

  //bool operator<(const Link& e2) const {
  //return (head < e2.GetHead() && tail < e2.GetTail());
  //}

  // for directed edge
  bool operator==(const Link& other) const {
    return (head==other.GetHead() && tail==other.GetTail());
  }

  //bool operator() (const Link& e1, const Link& e2) const {
  //return (e1.GetHead() < e2.GetHead() && e1.GetTail() < e2.GetTail()) ;
  //}

  friend ostream& operator<<(ostream& output, const Link& e) {
    output << "(" <<  e.GetHead() << "," << e.GetTail() <<")";
    return output;  
  }

 private:
  int head; // the starting node of the link
  int tail;   // the ending node of the link
  int index;  // index of this edge
  double weight; // weight of this edge

  bool color;  // 1 : head activates tail
  // 0 : head represses tail
};
//////////////////////////////////////////////////////////////////////////////////////////////////////////


#endif /* ! LINK_H_ */
