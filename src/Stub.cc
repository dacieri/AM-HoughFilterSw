#include "../interface/Stub.h"

Stub::Stub(){
	layerId_ = 0;
	phi_ = 0;
	r_ = 0;
	z_ = 0;
	tp_ = 0;
}

Stub::Stub(unsigned int layerId, float phi, float r, float z, int tp) : 
	layerId_(layerId),
	phi_(phi),
	r_(r),
	z_(z),
	tp_(tp)
{

}

Stub::~Stub() { }