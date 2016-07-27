#include "../interface/Settings.h"

Settings::Settings(){
	htType_ = 0;
	numXbins_ = 0;
	numYbins_ = 0;
	averageRes_ = 0;
	resPhi_ = 0;
	resZ0_ = 0;
	resTheta_ = 0;
	debug_ = 0;
	chosenRofPhi_ = 0.;
}


Settings::Settings(unsigned htType, unsigned xbins, unsigned ybins, unsigned averageRes, unsigned debug){
	htType_ = htType;
	numXbins_ = xbins;
	numYbins_ = ybins;
	averageRes_ = averageRes;
	// resPhi_ = 0.035;
	resPhi_ = 0.035;
	resZ0_ = 0;
	resTheta_ = 0;
	debug_ = debug;
	chosenRofPhi_ = 65.;
}