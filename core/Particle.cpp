#include "Particle.h"
#include <TDatabasePDG.h>
#include <iostream>

namespace elSpectro{

  Particle::Particle(int pdg):
    _pdg{pdg}{
    TDatabasePDG *pdgDB = TDatabasePDG::Instance();
    auto particle=pdgDB->GetParticle(pdg);

    if(particle==nullptr) {
      std::cerr<<"Particle::Particle pdg "<<pdg<<" does ntoe exist in table"<<std::endl;
      exit(0);
      
    }

    _pdgMass=particle->Mass();
    _dynamicMass=_pdgMass; //for stable particles
    
    SetXYZT(0,0,0,_pdgMass); //initialise at rest
    
  }
  ///////////////////////////////////////////
  void Particle::Print() const{

    std::cout<<"  Particle::Print() "<<std::endl;
 
    std::cout<<Pdg()<<" minimum mass = "<<MinimumMassPossible()<<" mass = "<<_dynamicMass<<"   PDG = "<<_pdgMass
	     <<"\n P4 ="<<P4().X()<<" "<<			\
      P4().Y()<<" "<<P4().Z()<<" "<<P4().T()<<" \n";
    
  }

}
