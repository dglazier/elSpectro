//////////////////////////////////////////////////////////////
///
///Class:		ProductionProcess
///Description:
///            Interface to different production processes
///            e.g Electroproduction, Photoproduction
///            1) model required to generate
///               reaction CoM state(from DecayingParticle)
#pragma once

#include "CurrentEventInfo.h"
#include "DecayModel.h"
#include "DecayingParticle.h"

namespace elSpectro{
  
 
  class ProductionProcess : public DecayingParticle {

  public:
    ProductionProcess()=delete;
    
    //only construct via and take ownership of model 
    ProductionProcess(DecayModel* model);
    ProductionProcess(int pdg,DecayVectors* decayer, DecayModel* model);

    
  

    
  private:

 
    
  };



}

