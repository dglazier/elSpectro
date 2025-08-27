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
    //_interactingParticle=P4ptr();
    _idxInteract = -1;
  }
  /////////////////////////////////////////////////////////
  //or a model to generate a 4-momentum
  CollidingParticle::CollidingParticle(int pdg,Double_t momentum,int parentpdg,decaymodel_ptr  model,decayer_ptr  decayer):
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
       
      if(p->Pdg()==pdg){
	
	if(_idxInteract!=-1){
	  std::cerr<<"CollidingParticle::CollidingParticle, multiple particles with pdg =  "<<pdg<<std::endl; exit(0);
	}
	else{	
	  //  _interactingParticle=p->P4ptr();
	  _idxInteract = position;
	}
      }
      
      position++;
    }
    
   //We need a "nominal" 4-momentum for our interacting particle
   //this is for integrations etc.
    //to do this we boost it from rest into lab frame of parent
   //Note theta= 0 from LorentzVector  lv(0,0,momentum,TMath::Sqrt(momentum*momentum+mass*mass));, so RandPhi does not matter
   if(_idxInteract!=-1){
     _nominal = GetInteracting4Vector();
     _decayer->BoostToParentWithRandPhi(P4(),_nominal);
     
   }
   
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
	prod->SetVertexID(VertexID());
      }
    
    _model->PostInit(info);
    }
    if(_decayer)_decayer->PostInit(info);
    
  }
  
}
