//////////////////////////////////////////////////////////////
///
///Class:		MassPhaseSpace
///Description:
///            Class to manage full reaction phase space
///            Should only be accessed via Manager::PhaseSpace;
///            Needs to be given the primary decay to calculate
///            the Mass phase space element for all children
///            and allocate masses for all particles in the chain
#pragma once

#include "DecayModel.h"

namespace elSpectro{

  class Manager;
  
  class MassPhaseSpace{

  public:

  private:
    
    friend Manager; //only Manager can construct a MassPhaseSpace
    MassPhaseSpace()=default;
    
    
    ClassDef(elSpectro::MassPhaseSpace,1); //class MassPhaseSpace
  };

}//namespace elSpectro
