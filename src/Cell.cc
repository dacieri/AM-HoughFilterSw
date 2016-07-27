#include "../interface/Cell.h"

void Cell::end(){
	const int maxLayerID(30);
	vector<bool> foundLayers(maxLayerID, false);

	for (Stub* stub: vStubs_) {
		foundLayers[stub->layerId()] = true;
	}

	unsigned int ncount = 0;
  	for (const bool& found: foundLayers) {
    	if (found) ncount++;
  	}

  	numLayers_ = ncount;
}