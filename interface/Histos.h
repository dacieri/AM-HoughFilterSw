#ifndef _HISTOS_H_
#define _HISTOS_H_

#include <vector>
#include "../interface/Stub.h"

using namespace std;

class Histos
{
public:
	Histos() {};
	~Histos() {};

	// Add stub to this cell in HT array.
	void store (Stub* stub) { vStubs_.push_back(stub); }
	void end();
	unsigned int numStubs() const { return vStubs_.size() }
	unsigned int numLayers() const { return numLayers_; }

	bool trackCandFound() const { return (numLayers_ >= minStubLayers_) }
	
private: 
	unsigned int numLayers_;
	std::vector<Stub*> vStubs_;

};

#endif
