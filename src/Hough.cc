#include "../interface/Hough.h"

using namespace std;
void Hough::init(Settings* settings, double xmin, double xmax, double ymin, double ymax, double charge){
	

	settings_ = settings;
	bField_ = 3.8112;
	double minTrackPt = 3.;

	// No Average Resolution
	if(settings->averageRes()==0){
		// Rphi Hough Transform
		if(settings->HoughType()==0){
			if(charge == -1){
				if(xmin > minTrackPt)
					xmin_ = charge/xmin ;
				else 
					xmin_ = charge/minTrackPt;
				xmax_ = charge/xmax ;
			} else if(charge == 1){
				xmin_ = charge/xmax ;
				if(xmin > minTrackPt)
					xmax_ = charge/xmin;
				else
					xmax_ = charge/minTrackPt;
			} else{
				if(xmin > minTrackPt){
					xmin_ = -1/xmin ;
					xmax_ = 1/xmin ;
				} else{
					xmin_ = -1/minTrackPt;
					xmax_ = 1/minTrackPt;
				}
			}
		} // Rz Hough Transform
		else {
			xmin_ = 1/tan(xmax);
			xmax_ = 1/tan(xmin);
		}
		if(ymin > -15.)
			ymin_ = ymin ;
		else
			ymin_ = -15.;

		if(ymax < 15.)
			ymax_ = ymax ;
		else
			ymax_ = 15.;
	} // Average Resolution 
	else{
		double xmean = (xmin+xmax)/2.;
		ymean_ = (ymin+ymax)/2.;
		// Rphi Hough Transform 
		if(settings->HoughType()==0){
			if(charge == -1){
				xmin_ = charge/(xmin);
				xmax_ = 0.;
			} else if(charge == 1){
				xmin_ = 0.;
				xmax_ = charge/xmin;
			} else{
				xmin_ = -1/xmin;
				xmax_ = 1/xmin;
			}
			ymin_ = ymean_ - settings->phiResolution();
			ymax_ = ymean_ + settings->phiResolution();
		} 
		else{
			xmin_ = 1./tan(xmax);
			xmax_ = 1./tan(xmin);
			ymin_ = ymean_ - settings->z0Resolution();
			ymax_ = ymean_ + settings->z0Resolution();
			// cout << "xmin_ "<< xmin_ << " xmax_ "<< xmax_ << endl;
			// cout << "ymin_ "<< ymin_ << " ymax_ "<< ymax_ << endl;
		}
	}

	xbins_ = settings->numXbins();
	ybins_ = settings->numYbins();
	ybinsize_ = (ymax_ - ymin_)/ybins_;
	xbinsize_ = (xmax_ - xmin_)/xbins_;

	if(settings_->HoughType()==0 && settings_->averageRes()==1)
		xbinsize_ = (2./3.)/xbins_;
	else if(settings_->HoughType()==1 && settings_->averageRes()==1)
		xbinsize_ = ((1./tan(0.8))-(1./tan(2.)))/xbins_;

	invPtToDphi_ = bField_*(3.0E8/2.0E11);
	chosenRofPhi_ = settings_->chosenRofPhi();
	htArray_.resize(xbins_, ybins_, false);
}

void Hough::store(Stub* stub){

	InputStubs_.push_back(stub);

	for(unsigned int xbin = 0; xbin<xbins_ ; ++xbin){
		// In this q/Pt bin, find the range of phi bins that this stub is consistent with.
		bool goodXbin = true;
		// rphi Hough Transform with average Resolution
		if(settings_->HoughType()==0 && settings_->averageRes()==1){
			float qOverPtBin = -1./3 + (xbin + 0.5) * xbinsize_;
			if(qOverPtBin < xmin_ || qOverPtBin > xmax_)
				goodXbin = false;
		} // rz Hough Transform with average Resolution
		else if(settings_->HoughType()==1 && settings_->averageRes()==1){
			float TanThetaBin = 1/tan(2.) + (xbin + 0.5) * xbinsize_;
			if(TanThetaBin < xmin_ || TanThetaBin > xmax_){
				goodXbin = false;
			}
		}

		pair<unsigned int, unsigned int> iRange;
		// int ybin = this->PhiBin(stub, xbin);
		if(settings_->HoughType() == 0)
			iRange = this->iPhiRange( stub, xbin );
		else
			iRange = this->iZ0Range( stub, xbin );

		unsigned int yBinMin = iRange.first;
		unsigned int yBinMax = iRange.second;

		// cout << " yBinMin "<< yBinMin << " yBinMax "<< yBinMax << endl;

		if(goodXbin){
			for(unsigned int ybin = yBinMin; ybin <= yBinMax; ++ybin){
				htArray_(xbin, ybin).store(stub);
				// cout << "stub stored" << endl;
			}
		}
	}



}

bool Hough::trackCandFound(){
	unsigned int numStubs = 0;
	bool candFound = false;
	for(unsigned int xbin = 0; xbin<xbins_ ; ++xbin){
		for(unsigned int ybin = 0; ybin<ybins_ ; ++ybin){
			htArray_(xbin,ybin).end();
			if(htArray_(xbin,ybin).numStubs() > numStubs){
				FilteredStubs_ = htArray_(xbin,ybin).stubsInCell();
				numStubs = htArray_(xbin,ybin).numStubs();
				if(htArray_(xbin,ybin).trackCandFound()){
					candFound = true;
					// cout << "candidate found "<< endl;
				}
			}
			if(htArray_(xbin,ybin).trackCandFound() && settings_->debug()==1 && InputStubs_.size() > 15){
				cout << "xbin "<<xbin<< " ybin "<< int(ybin) - int(ybins_/2)<< endl;
				for(Stub* stub : htArray_(xbin,ybin).stubsInCell()){
					cout << "stub phi "<< stub->phi() <<" phiR "<< stub->phi()-ymean_ << " stub rT "<< stub->r()-chosenRofPhi_ << endl;
				}
			}
			

		}
	}
	// cout << settings_->debug() << endl;
	if(candFound && settings_->debug()==1 && InputStubs_.size() > 15){
		cout << InputStubs_.size() << " input stubs, phi Road " << ymean_ << " mMin " << floor((xmin_+1/3.)/xbinsize_) - xbins_/2 << " mMax "<< floor((xmax_+1/3.)/xbinsize_) - xbins_/2 << endl;
		for(Stub* stub : InputStubs_){
			cout << "stub phi "<< stub->phi() <<" phiR "<< stub->phi()-ymean_ << " stub rT "<< stub->r()-chosenRofPhi_;
			cout << " integers phi "<< floor((stub->phi()-ymean_)*14628.57143) << " r: "<<floor((stub->r()-chosenRofPhi_)*6.969051429)<< endl;
		}

		for(Stub* stub : InputStubs_){
		cout << "inputStub.r <= std_logic_vector(to_signed("<< floor((stub->r()-chosenRofPhi_)*6.969051429)<<", rWidth));"<< endl << "inputStub.phi <= std_logic_vector(to_signed("<<floor((stub->phi()-ymean_)*14628.57143)<<",phiWidth));"<< endl << "inputStub.layerId <= std_logic_vector(to_unsigned("<<stub->layerId()-5<<", layerIdWidth));"<< endl <<"wait for clk_period;"<< endl;
		}
	}

	return candFound;
}

std::pair<unsigned int, unsigned int> Hough::iPhiRange(Stub* stub, unsigned int xbin){
  // Note q/Pt value corresponding to centre of this bin.
  float qOverPtBin;
  if(settings_->averageRes() == 0)
  	qOverPtBin = xmin_ + (xbin + 0.5) * xbinsize_;
  else
  	qOverPtBin = -1./3 + (xbin + 0.5) * xbinsize_;
  // Note change in this q/Pt value needed to reach either edge of the bin. 
  float qOverPtBinVar = 0.5*xbinsize_;

  // Reducing effective bin width can reduce fake rate.
  //qOverPtVar = 0.4*binSizeQoverPtAxis_;

  // Calculate range of track-phi that would allow a track in this q/Pt range to pass through the stub.
  float phiTrk = stub->phi() + invPtToDphi_ * qOverPtBin * (stub->r() - chosenRofPhi_);
  // The next line does the phiTrk calculation without the usual approximation, but it doesn't 
  // improve performance.
  // float phiTrk    = stub->phi() + asin(invPtToDphi_ * qOverPtBin * stub->r()) - asin(invPtToDphi_ * qOverPtBin * chosenRofPhi_);
  float phiTrkVar =  invPtToDphi_ * qOverPtBinVar * fabs(stub->r() - chosenRofPhi_);
  float phiTrkMin = phiTrk - phiTrkVar;
  float phiTrkMax = phiTrk + phiTrkVar;
  

  int phiBinMin = floor( ( phiTrkMin - ymin_ ) / ybinsize_ );
  int phiBinMax = floor( ( phiTrkMax - ymin_ ) / ybinsize_ );	
  

  // Limit range to dimensions of HT array.
  phiBinMin = max(phiBinMin, 0);
  phiBinMax = min(phiBinMax, int(ybins_) - 1);

  // If whole range is outside HT array, flag this by setting range to specific values with min > max.
  if (phiBinMin > int(ybins_) - 1 || phiBinMax < 0) {
    phiBinMin = int(ybins_) - 1;
    phiBinMax = 0;
  }

  std::pair<unsigned int, unsigned int> phiRange(phiBinMin, phiBinMax);

  return phiRange;	
}

int Hough::PhiBin(Stub* stub, unsigned int xbin){
	float qOverPtBin;
	if(settings_->averageRes() == 0)
    	qOverPtBin = xmin_ + (xbin + 0.5) * xbinsize_;
	else
	  	qOverPtBin = -1./3 + (xbin + 0.5) * xbinsize_;

    float phiTrk = stub->phi() + invPtToDphi_ * qOverPtBin * (stub->r() - chosenRofPhi_);
    int phiBin = floor( (phiTrk - ymin_)/ybinsize_ );
    if(phiBin > int(ybins_) -1 || phiBin < 0){
    	phiBin = -1; 
    }

    return phiBin;
}

std::pair<unsigned int, unsigned int> Hough::iZ0Range(Stub* stub, unsigned int xbin){
	// Note tan(theta) value corresponding to centre of this bin.
	float cotanThetaBin;
	if(settings_->averageRes() == 0)
  		cotanThetaBin = xmin_ + (xbin + 0.5) * xbinsize_;
  	else 
  		cotanThetaBin = 1/tan(2.) + (xbin + 0.5) * xbinsize_;
  // Note change in this q/Pt value needed to reach either edge of the bin. 
  float cotanThetaBinVar = 0.5*xbinsize_;

  // Reducing effective bin width can reduce fake rate.
  //qOverPtVar = 0.4*binSizeQoverPtAxis_;

  // Calculate range of track-phi that would allow a track in this q/Pt range to pass through the stub.
  float z0 = stub->z() - cotanThetaBin * stub->r();
  // The next line does the phiTrk calculation without the usual approximation, but it doesn't 
  // improve performance.
  // float phiTrk    = stub->phi() + asin(invPtToDphi_ * qOverPtBin * stub->r()) - asin(invPtToDphi_ * qOverPtBin * chosenRofPhi_);
  float z0Var = cotanThetaBinVar * fabs(stub->r()) + stub->zErr();
  float z0Min = z0 - z0Var;
  float z0Max = z0 + z0Var;
  
  int z0BinMin = floor( ( z0Min - ymin_ ) / ybinsize_ );
  int z0BinMax = floor( ( z0Max - ymin_ ) / ybinsize_ );	

  // Limit range to dimensions of HT array.
  z0BinMin = std::max(z0BinMin, 0);
  z0BinMax = std::min(z0BinMax, int(ybins_) - 1);

  // If whole range is outside HT array, flag this by setting range to specific values with min > max.
  if (z0BinMin > int(ybins_) - 1 || z0BinMax < 0) {
    z0BinMin = int(ybins_) - 1;
    z0BinMax = 0;
  }

  // cout << "z0BinMin "<< z0BinMin << " z0BinMax "<< z0BinMax << endl;


  std::pair<unsigned int, unsigned int> z0Range(z0BinMin, z0BinMax);

  return z0Range;	
}
