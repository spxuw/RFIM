#ifndef _POINT_H_
#define _POINT_H_

#include <cmath>
#include <iostream>
#include <vector>
using namespace std;

///////////////////////////////////////////////////////
// Note that this point class is only valid for D<=3. 
// and I limit it useness in integers.
class Point 
{
public:
    Point(int x=0, int y=0, int z=0)                    // constructor
	{ 
	    m_x = x; 
	    m_y = y; 
	    m_z = z; 
	} 
    
   // Pay attention to that Coords[D-1] is the x component!
    Point(int dimension, int* Coords)                    // constructor
	{ 
	  
	    switch(dimension)
	    {
		case 1: 
		    m_x = Coords[0];////size[0]   = W;    stride[0] = 1;  //x
		    m_y = 0;
		    m_z = 0;
		    break;

		case 2:
		    m_x = Coords[1];// size[1]   = W;    stride[1] = 1;  //x
		    m_y = Coords[0];//size[0]   = L;    stride[0] = W;  //y
		    m_z = 0; 
		    break;
		    
		case 3: 
		    m_x = Coords[2];//size[2]   = W;    stride[2] = 1;  //x
		    m_y = Coords[1];//size[1]   = L;    stride[1] = W;  //y
		    m_z = Coords[0];//size[0]   = H;    stride[0] = W*L;//z
		    break;
	    }
  
	} 

    void Move(int, int, int);                          // Move to position (r,s)
    void Move(Point);                                            // Move to position (r,s)

    void Shift(int, int, int);                          // Move to position (r,s)
    void Shift(int, vector<int>);  

    int x() const;
    int y() const;
    int z() const;

    int coord(int, int) const;

    int GetIxx() const;
    int GetIxy() const;
    int GetIxz() const;
    int GetIyy() const;
    int GetIyz() const;
    int GetIzz() const;

    void Draw() const;                                          // a const member 
    friend double Norm(Point);                                   // a friend
	 
private:
    int m_x;
    int m_y;
    int m_z;
};


inline void Point::Shift(int x=0, int y=0, int z=0)  
{
    m_x -= x; 
    m_y -= y; 
    m_z -= z; 
}

inline void Point::Move(int x=0, int y=0, int z=0)  
{
    m_x = x; 
    m_y = y; 
    m_z = z; 
}

inline void Point::Move(Point P)
{
    m_x = P.m_x; 
    m_y = P.m_y; 
    m_z = P.m_z; 
}


inline int Point::x() const
{
    return m_x; 
}

inline int Point::y() const
{
    return m_y; 
}

inline int Point::z() const
{
    return m_z; 
}

inline int Point::coord(int D, int i) const
{
    int temp=0;

    if(D==1)
    {
	temp = m_x;
    }
    else if(D==2)
    {
	switch(i) {
	    case 1: temp = m_x; break;
	    case 0: temp = m_y; break; 
	}
    }
    else if(D==3)
    {
	switch(i) {
	    case 2: temp = m_x; break;
	    case 1: temp = m_y; break;
	    case 0: temp = m_z; break; 
	}
    }

    return temp; 
}



inline void Point::Shift(int D, vector<int> deltCoords)  
{
    switch(D)
    {
	case 1:
	    m_x -= deltCoords[0];
	    break;
	case 2:
	    m_x -= deltCoords[1];// size[1]=W;  stride[1]=1; 
	    m_y -= deltCoords[0];// size[0]=L;  stride[0]=W;
	    break;

	case 3:
	    m_x -= deltCoords[2];// size[2]=W; stride[2]=1;
	    m_y -= deltCoords[1];// size[1]=L; stride[1]=W;
	    m_z -= deltCoords[0];// size[0]=H  stride[0]=W*L;
	    break;
    }
}







inline int Point::GetIxx() const
{
    return (m_y*m_y + m_z*m_z); 
}


inline int Point::GetIxy() const
{
    return (-m_x*m_y); 
}


inline int Point::GetIxz() const
{
    return (-m_x*m_z); 
}


inline int Point::GetIyy() const
{
    return (m_x*m_x + m_z*m_z); 
}


inline int Point::GetIyz() const
{
    return (-m_y*m_z); 
}


inline int Point::GetIzz() const
{
    return (m_x*m_x + m_y*m_y); 
}


inline void Point::Draw() const 
{
    cout << '(' << m_x << ',' << m_y << ',' << m_z << "),";
}

inline double Norm(Point P)
{
    return sqrt(double(P.m_x*P.m_x + P.m_y*P.m_y + P.m_z*P.m_z ));
}


#endif // _POINT_H_
