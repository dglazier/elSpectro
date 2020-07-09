///////////////////////////////////////////////////
///
///    Functions to simplify interfacing with managers
#pragma once

#include "Manager.h"
//#include "DecayModelWQ2.h"
//#include "ParticleManager.h"
//#include "ParticleManager.h"

namespace elSpectro{
  
   //////////////////////////////////////////////////////////////
  Manager& generator(){return Manager::Instance();}
  ParticleManager& particles(){return generator().Particles();}

  //////////////////////////////////////////////////////////////
  Particle* particle(int pdg){
    return particles().Take(new Particle{pdg});
  }
  Particle* particle(int pdg,DecayModel* const model){
    return particles().Take(new DecayingParticle{pdg,model});
  }
  void mass_distribution(int pdg,Distribution *dist){
    particles().RegisterMassDistribution(pdg,dist);
  }
   //////////////////////////////////////////////////////////////
  //note there is a std::decay function, so call this model
  //const =>can change the model but not th pointer
  DecayModel* const model(DecayModel* model){
    return generator().Decays().Take(model);
  }
  
  //////////////////////////////////////////////////////////////
  ProductionProcess* const reaction(ProductionProcess* prod){
    generator().Reaction(prod);
    return  generator().Reaction();
  }
  ProductionProcess* const eic(double ep,double ionp,DecayModelQ2W *totalXsec=nullptr,int ionpdg=2212){
    //computational limits seem to be 1E-6
    auto dist =  new DistVirtPhotFlux_xy(ep,1E-6,1,1E-6,1);
    //ElectronScattering(electronP,ionP,electronAngle,ionAngle)
    if(totalXsec!=nullptr)
      generator().Reaction(new ElectronScattering(ep,ionp,TMath::Pi(),0,new ScatteredElectron_xy(dist),totalXsec) );
    else
      generator().Reaction(new ElectronScattering(ep,ionp,TMath::Pi(),0,new ScatteredElectron_xy(dist)) );

     
    return  generator().Reaction();
  }
  /*  ProductionProcess* const eic(double ep,double ionp,DecayModelWQ2 *totalXsec,int ionpdg=2212){
    auto dist =  new DistVirtPhotFlux_xy(ep,1E-10,1,1E-4,1);   
    generator().Reaction(new ElectronScattering(ep,0,0,TMath::Pi(),new ScatteredElectron_xy(dist)),totalXsec);
    return  generator().Reaction();
    }*/
 
  ProductionProcess* const mesonex(double ep,DecayModelQ2W *totalXsec=nullptr,int ionpdg=2212){
    auto dist =  new DistVirtPhotFlux_xy(ep,1E-10,1,1E-4,1);  
    if(totalXsec!=nullptr)
      generator().Reaction(new ElectronScattering(ep,0,0,TMath::Pi(),new ScatteredElectron_xy(dist),totalXsec ));
    else	
      generator().Reaction(new ElectronScattering(ep,0,0,TMath::Pi(),new ScatteredElectron_xy(dist)));
	
    return  generator().Reaction();
  }
  
}//namespace elSpectro
