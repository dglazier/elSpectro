//////////////////////////////////////////////////////////////
///
///Class:		NBodyPhaseSpace
///Description:
///            Class to manage full reaction phase space
///            Should only be accessed via Manager::PhaseSpace;
///            Needs to be given the primary decay to calculate
///            the Mass phase space element for all children
///            and allocate masses for all particles in the chain
#pragma once

#include "DecayModel.h"
#include <TRandom.h>
#include <TBenchmark.h>

namespace elSpectro{


  void SampleNBodyPhaseSpace(double parentM, DecayModel* model);
  
}//namespace elSpectro
