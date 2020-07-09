//////////////////////////////////////////////////////////////
///
///Class:		DecayModelQ2W
///Description:
///             Control behaviour of Particle decay to Particle products
///             Defined by
///             1) list of Particle products
///             2) Intensity function given the produt LorentzVectors
///
///            Note drevied classes should include a constructor to initialise
///            DecayModelQ2W( particle_ptrs , const std::vector<int> pdgs );
#pragma once

#include "DecayModel.h"

namespace elSpectro{

 
  class DecayModelQ2W : public DecayModel {

  public:
    
    DecayModelQ2W()=delete;
    //delete default constructor so have to use threshold
    //so other 5 constructors also defaulted(rule of 5)

    //constructor giving W theshold, just produces scatted electron kinematics
    DecayModelQ2W(  double thresh  );
    //constructor giving W theshold and subsequent primary decay of Nucl+gamma* system
    DecayModelQ2W(  double thresh, DecayModel* primary  );
    
    // Each model must define its intensity
    const CurrentEventInfo* Intensity(const CurrentEventInfo* info=nullptr) const override;

    
  private:
    void Init();
    
    mutable PhotoProdInfo _myInfo;//!
    double _threshold = {0};
    
    Particle* _electron={nullptr};
    Particle* _gstarNuc={nullptr};
      
    ClassDefOverride(elSpectro::DecayModelQ2W,1); //class DecayModelQ2W
    
  };//class DecayModelQ2W

}//namespace elSpectro
