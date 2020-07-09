//////////////////////////////////////////////////////////////
///
///Class:		DecayManager
///Description:
///            Class to manage decays
///            Should only be accessed via Manager::Decays();
///            Decays can be of many types (derived classes) so store unique_ptr
///             1) Keep ownership of all decays
#pragma once

#include "DecayModel.h"

namespace elSpectro{

  class Manager;
  using decay_model_uptr = std::unique_ptr<DecayModel>;

  
  class DecayManager{

  public:
    DecayModel* Take(DecayModel* p){
      _decays.push_back(std::move(decay_model_uptr{p}));
      return _decays.back().get();
    }
    
  private:
    
    friend Manager; //only Manager can construct a DecayManager
    DecayManager()=default;
    
    std::vector<decay_model_uptr> _decays;
    
    ClassDef(elSpectro::DecayManager,1); //class DecayManager
  };

}//namespace elSpectro
