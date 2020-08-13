///////////////////////////////////////////////////
///
///    Functions to simplify interfacing with managers
#pragma once

#include "Manager.h"
#include "DecayModelQ2W.h"
#include "Distribution.h"
#include "DistVirtPhotFlux_xy.h"
#include "ElectronScattering.h"
#include "ScatteredElectron_xy.h"
#include <TDatabasePDG.h>

//#include "ParticleManager.h"
//#include "ParticleManager.h"

namespace elSpectro{
  
   //////////////////////////////////////////////////////////////
  inline Manager& generator(){return Manager::Instance();}
  //////////////////////////////////////////////////////////////

  inline void nextEvent(){
    generator().Clear();
    generator().Reaction()->GenerateProducts();
    generator().Write();
  }
  //////////////////////////////////////////////////////////////
  inline ParticleManager& particles(){return generator().Particles();}

  //////////////////////////////////////////////////////////////
  inline Particle* particle(int pdg){
    return particles().Take(new Particle{pdg});
  }
  //////////////////////////////////////////////////////////////
  inline Particle* particle(int pdg,DecayModel* const model){
    return particles().Take(new DecayingParticle{pdg,model});
  }
  //////////////////////////////////////////////////////////////
  inline void mass_distribution(int pdg,Distribution *dist){
    particles().RegisterMassDistribution(pdg,dist);
  }
   //////////////////////////////////////////////////////////////
  //note there is a std::decay function, so call this model
  //const =>can change the model but not th pointer
  inline DecayModel* model(DecayModel* model){
    return generator().Decays().Take(model);
  }
  //////////////////////////////////////////////////////////////
  inline ProductionProcess* reaction(ProductionProcess* prod){
    generator().Reaction(prod);
    return  generator().Reaction();
  }
  //////////////////////////////////////////////////////////////
  inline void writer(Writer* wr){
    generator().SetWriter(wr);
  }
  /////////////////////////////////////////////////////////////
  inline void initGenerator(){
    generator().InitGeneration();
  }
  //////////////////////////////////////////////////////////////
  inline ProductionProcess* eic(double ep,double ionp,DecayModelQ2W *totalXsec=nullptr,int ionpdg=2212){
    //computational limits seem to be 1E-6
    //auto dist =  new DistVirtPhotFlux_xy(ep,1E-6,1,1E-6,1);
    // auto dist =  new DistVirtPhotFlux_xy(ep,1E-16,1,1E-8,1);
  
    //ElectronScattering(electronP,ionP,electronAngle,ionAngle)
    if(totalXsec!=nullptr)
      generator().Reaction(new ElectronScattering(ep,ionp,TMath::Pi(),0,totalXsec) );
    else
      generator().Reaction(new ElectronScattering(ep,ionp,TMath::Pi(),0) );

     
    return  generator().Reaction();
  }
  /*  ProductionProcess* const eic(double ep,double ionp,DecayModelWQ2 *totalXsec,int ionpdg=2212){
    auto dist =  new DistVirtPhotFlux_xy(ep,1E-10,1,1E-4,1);   
    generator().Reaction(new ElectronScattering(ep,0,0,TMath::Pi(),new ScatteredElectron_xy(dist)),totalXsec);
    return  generator().Reaction();
    }*/
 
  inline ProductionProcess*  mesonex(double ep,DecayModelQ2W *totalXsec=nullptr,int ionpdg=2212){
  
    // auto dist =  new DistVirtPhotFlux_xy(ep,1E-10,1,1E-4,1);  
    if(totalXsec!=nullptr)
      generator().Reaction(new ElectronScattering(ep,0,0,TMath::Pi(),totalXsec ));
    else	
      generator().Reaction(new ElectronScattering(ep,0,0,TMath::Pi()));
	
    return  generator().Reaction();
  }
  
}//namespace elSpectro
