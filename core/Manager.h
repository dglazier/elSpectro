//////////////////////////////////////////////////////////////
///
///Class:		Manager
///Description:
///            Class to manage manangers!
///           1) Access ParticleManager via Manager::Instance()->Particles()
///           2) Access DecayManager via Manager::Instance()->Decays()
///           3) Access ProductionProcess via Manager::Instance()->Process()
#pragma once

#include "ParticleManager.h"
#include "DecayManager.h"
#include "ProductionProcess.h"
#include "MassPhaseSpace.h"
#include <TRandom3.h>

namespace elSpectro{

  class Manager{

  public:
    
     static Manager& Instance() { static Manager instance; return instance; }

     ParticleManager& Particles() noexcept{return _particles;}
     DecayManager& Decays() noexcept{return _decays;}
     
     void Reaction(ProductionProcess* prod){
       _process.reset(prod);
     }
     
     ProductionProcess* const Reaction(){return _process.get();}

     void SetSeed(ULong_t seed = 0){gRandom->SetSeed(seed);}


     void SetModelForMassPhaseSpace(DecayModel* amodel){_massPhaseSpace.SetModel(amodel);}
     void  FindMassPhaseSpace(double parentM,const  DecayModel* amodel) {
       _massPhaseSpace.Find(parentM,amodel);
     }
  private:

    ParticleManager _particles;
    DecayManager _decays;

    std::unique_ptr<ProductionProcess> _process;

    MassPhaseSpace _massPhaseSpace;
      
    ClassDef(elSpectro::Manager,1); //class Manager
  };

}
