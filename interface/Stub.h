#ifndef _STUB_H_
#define _STUB_H_

#include <vector>

using namespace std;

class Stub
{
public:
	Stub();

    Stub(
  		unsigned int layerId,
  		float phi,
  		float r,
  		float z,
      int tp
	);

    unsigned int layerId() 		{ return layerId_; }
    float z()			 		{ return z_;       }
    float r()					{ return r_; 	     }
    float phi()			  { return phi_;	   }
    int tp()          { return tp_;      }
    float zErr()      { return (layerId_-5 < 3) ? 0 : 2.5;}

	~Stub();
	
private: 
	unsigned int layerId_;
  float phi_;
  float r_;
	float z_;
  int   tp_;
};

#endif
