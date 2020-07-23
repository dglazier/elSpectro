//////////////////////////////////////////////////////////////
///
///Class:		PhaseSpaceDecay
///Description:
///             Calculate intensity for phase space decay (=1 )!

#pragma once

#include "DecayModel.h"

namespace elSpectro{

 
  class PhaseSpaceDecay : public DecayModel {

  public:
    
    PhaseSpaceDecay()=default;
    //only declaring default constructor
    //so other 5 constructors also defaulted(rule of 5)
    //constructor to decay into particles
    PhaseSpaceDecay( particle_ptrs , const std::vector<int> pdgs );

    // Each model must define its intensity
    double Intensity() const final{
      if(CheckThreshold())
	return 1.;
      else
	return 0.;
    }
    bool RegenerateOnFail() const  noexcept final {return false;}

    
  private:
 
    mutable CurrentEventInfo _myInfo;//!

    ClassDef(elSpectro::PhaseSpaceDecay,1); //class PhaseSpaceDecay
    
  };//class PhaseSpaceDecay

}//namespace elSpectro
