#ifndef _HOUGH_H_
#define _HOUGH_H_

#include <vector>
#include "boost/numeric/ublas/matrix.hpp"
#include <utility>
#include "../interface/Cell.h"
#include "../interface/Stub.h"
#include "../interface/Settings.h"
#include <math.h> 
#include <algorithm>
#include <iostream>
#include <string>


using  boost::numeric::ublas::matrix;
using namespace std;

class Hough
{
public:
	Hough() {};
	~Hough() {};

	void init(Settings* settings, double xmin, double xmax, double ymin, double ymax, double charge);
	void store(Stub* stub);
	std::pair<unsigned int, unsigned int> iPhiRange(Stub* stub, unsigned int xbin);
	int PhiBin(Stub* stub, unsigned int xbin);
	std::pair<unsigned int, unsigned int> iZ0Range(Stub* stub, unsigned int xbin);

	std::vector<Stub*> InputStubs()		{return InputStubs_;}
	std::vector<Stub*> FilteredStubs() { return FilteredStubs_;}
	unsigned int maxNumStubs() 			{ return FilteredStubs_.size();}
	bool trackCandFound();

	double xMin() { return xmin_;}
	double xMax() { return xmax_;}
	double yMin() { return ymin_;}
	double yMax() { return ymax_;}
private: 
	matrix<Cell> htArray_; 
	Settings* settings_;
	
	double ymin_;
	double ymax_;
	double ymean_;
	double xmin_;
	double xmax_;
	
	unsigned int xbins_;
	unsigned int ybins_;
	double xbinsize_;
	double ybinsize_;
	double bField_;
	double invPtToDphi_;
	double chosenRofPhi_;
	
	std::vector<Stub*> InputStubs_;
	std::vector<Stub*> FilteredStubs_;

};

#endif
