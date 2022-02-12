#ifndef _SPLITTER_H_
#define _SPLITTER_H_
#include <string>
#include <vector>
using namespace std;

// Maintains a collection of substrings that are delimited by a string of one or more characters
class Splitter {
  // Contains the split tokens
  vector<string> _tokens;
public:
  // Subscript type for use with operator[]
  typedef  vector<string>::size_type size_type;
public:
  // Create and initialize a new Splitter
  //
  // param[in] src The string to split
  // param[in] delim The delimiter to split the string around
  Splitter(const  string& src, const  string& delim ) {
    reset(src, delim );
  }
  // Retrieve a split token at the specified index
  // param[in] i The index to search for a token
  // return The token at the specified index
  // throw  out_of_range If the index is invalid
  string& operator[] (size_type i ) {
    return _tokens.at(i);
  }
  // Retrieve the number of split tokens
  // return The number of split tokesn
  size_type size() const {
    return _tokens.size();
  }
  // Re-initialize with a new soruce and delimiter
  //
  // param[in] src The string to split
  // param[in] delim The delimiter to split the string around
  void reset (const string& src, const string& delim) {
    vector< string> tokens;
    string::size_type start = 0;
    string::size_type end;
    for ( ; ; ) {
      end = src.find ( delim, start );
      tokens.push_back ( src.substr ( start, end - start ) );
      // We just copied the last token
      if ( end ==  string::npos )
	break;
      // Exclude the delimiter in the next search
      start = end + delim.size();
    }
    _tokens.swap ( tokens );
  }
};


#endif
