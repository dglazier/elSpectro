#include "CollidingParticle.h"

namespace elSpectro{

  ///////////////////////////////////////////////////////////
  //cannot default construct at least need a 4-momentum
  CollidingParticle::CollidingParticle(int pdg,LorentzVector lv):
    Particle(pdg){
    SetP4(lv);
    _interactingParticle=P4ptr();
  }
  /////////////////////////////////////////////////////////
  //or a model to generate a 4-momentum
  CollidingParticle::CollidingParticle(int pdg,int parentpdg,DecayModel* model,DecayVectors* decayer):
    Particle(parentpdg),_model{model},_decayer{decayer}
  {
    //Find the relevent particle pointer in the model
    for(auto& p:_model->Products()){
      if(p->Pdg()==pdg){
	if(_interactingParticle!=nullptr){
	  std::cerr<<"CollidingParticle::CollidingParticle, multiple particles with pdg =  "<<pdg<<std::endl; exit(0);
	}
	else	
	  _interactingParticle=p->P4ptr();

      }
    }
  }
  
}
