#ifndef _CELL_H_
#define _CELL_H_

#include <vector>
#include "../interface/Stub.h"

using namespace std;

const unsigned int minStubLayers_ = 5;


class Cell
{
public:
	Cell() {};
	~Cell() {};

	// Add stub to this cell in HT array.
	void store (Stub* stub) { vStubs_.push_back(stub); }
	void end();
	unsigned int numStubs() const { return vStubs_.size(); }
	unsigned int numLayers() const { return numLayers_; }

	bool trackCandFound() const { return (numLayers_ >= minStubLayers_); }
	std::vector<Stub*> stubsInCell() const { return vStubs_;}
	
private: 
	unsigned int numLayers_;
	std::vector<Stub*> vStubs_;
};

#endif
