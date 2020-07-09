#include "ProductionProcess.h"

namespace elSpectro{

   //////////////////////////////////////////////////////////////////
  ProductionProcess::ProductionProcess(DecayModel* model):
    DecayingParticle{model}
  {


  }
  ProductionProcess::ProductionProcess(int pdg,DecayVectors* decayer, DecayModel* model):
    DecayingParticle{pdg,decayer,model}{

  }
}
