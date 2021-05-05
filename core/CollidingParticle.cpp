#include "CollidingParticle.h"
#include "FunctionsForGenvector.h"
#include "Manager.h"

namespace elSpectro{

  ///////////////////////////////////////////////////////////
  //cannot default construct at least need a 4-momentum
  CollidingParticle::CollidingParticle(int pdg,Double_t momentum):
    Particle(pdg){
    //Set interacting particle pdg
    _interactingPdg=pdg;
    auto mass= PdgMass();
    LorentzVector lv(0,0,momentum,TMath::Sqrt(momentum*momentum+mass*mass));
    SetP4(lv);
    _nominal=P4();
    _interactingParticle=P4ptr();
  }
  /////////////////////////////////////////////////////////
  //or a model to generate a 4-momentum
  CollidingParticle::CollidingParticle(int pdg,Double_t momentum,int parentpdg,DecayModel* model,DecayVectors* decayer):
    Particle(parentpdg),_model{model},_decayer{decayer}
  {
    //Set interacting particle pdg
    _interactingPdg=pdg;
    //set the parent LorentzVector
    auto mass= PdgMass();
    LorentzVector lv(0,0,momentum,TMath::Sqrt(momentum*momentum+mass*mass));
    SetP4(lv);
   _nominal=P4();
  
   //Find the relevent particle pointer in the model
   //this is the particle which is used in the production process
   UInt_t position=0;
   for(auto& p:_model->Products()){
       
     std::cout<<"CollidingParticle prod "<<p->Pdg()<<std::endl;
      if(p->Pdg()==pdg){
	
	if(_interactingParticle!=nullptr){
	  std::cerr<<"CollidingParticle::CollidingParticle, multiple particles with pdg =  "<<pdg<<std::endl; exit(0);
	}
	else{	
	  _interactingParticle=p->P4ptr();
	}
	//Not a stable final state particle!
	Manager::Instance().Particles().RemoveStable(p);
	_model->SwapProducts(0,position); //make sure interacting particle is first in products vector, must be done after remove
     }
      else{ //not a beam particle but "spectator"
	//remove from stable particle list
	//as these will be boosted to lab
	//whereas colliding particles and their
	//products are already in lab frame
	
	Manager::Instance().Particles().MoveStableToLab(p);

      }
      position++;
    }
    //move interacting particle to start of products vector
    
    //We need a "nominal" 4-momentum for our interacting particle
    //to do this we boost it from rest into lab frame of parent
    _decayer->BoostToParent(P4(),(*_interactingParticle));
  }
  /////////////////////////////////////////////////////////
  /// rotate beam angles
  void CollidingParticle::SetAngleThetaPhi(Double_t th,Double_t phi){
      _dirTheta=th;_dirPhi=phi;
      auto p4=_nominal; //copy 4-vector
      //and rotate it
      genvector::LorentzRotateY(p4,th);
      genvector::LorentzRotateZ(p4,phi);
      //Set the rotated vector
      SetP4(p4);
      _nominal=p4;
    }
  /////////////////////////////////////////////////////////
  void CollidingParticle::PostInit(ReactionInfo* info){
        //decay vertex position
    
    if(_model){
      auto& products=_model->Products();
      //same vertex as parent
      for(auto* prod: products){
	prod->SetVertex(VertexID(),VertexPosition());
      }
    
    _model->PostInit(info);
    }
    if(_decayer)_decayer->PostInit(info);
    
  }
  
}
