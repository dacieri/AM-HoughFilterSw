#ifndef _SETTINGS_H_
#define _SETTINGS_H_

#include <vector>
#include <fstream>

using namespace std;

class Settings
{
public:
	Settings();
	Settings(unsigned htType, unsigned xbins, unsigned ybins, unsigned averageRes, unsigned debug);
	~Settings() {};

	unsigned HoughType() { return htType_; } // Return the type of Hough Transform to be performed (0: rphi, 1: rz)
	unsigned numXbins() { return numXbins_; } // Number of x bins of HT array
	unsigned numYbins() { return numYbins_; } // Number of x bins of HT array
	unsigned averageRes() { return averageRes_; } // If 1 use average resolution 
	double   phiResolution() { return resPhi_;} // Average resolution in phi
	double   z0Resolution()  { return resZ0_; } // Average resolution in z0
	double   thetaResolution() { return resTheta_; } // Average resolution in theta
	unsigned debug() { return debug_; } // Print some debug information
	double chosenRofPhi() {return chosenRofPhi_; }

private: 

	// definition of cfg parameters
	unsigned htType_;
	unsigned numXbins_;
	unsigned numYbins_;
	unsigned averageRes_;
	double   resPhi_;
	double   resZ0_;
	double   resTheta_;
	unsigned debug_;
	double chosenRofPhi_;	


};

#endif
