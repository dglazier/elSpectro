//////////////////////////////////////////////////////////////
///
///Class:		CollidingParticle
///Description:
///             1) Class for generating initial state LorentzVector
/// 


#pragma once

#include "Particle.h"
#include "DecayChannel.h"
//#include "DecayVectors.h"

namespace elSpectro{
  
 
  class CollidingParticle : public Particle {

  public:
    
    CollidingParticle()=delete;
    //cannot default construct at least need a 4-momentum
    CollidingParticle(int pdg,Double_t momentum);
    //or a model to generate a 4-momentum
    CollidingParticle(int pdg,Double_t momentum,int parentpdg,decaymodel_ptr  model,decayer_ptr  decayer);
 
    DecayModel*  Model()const {return _model.get();}

    double Generate(){
       if(_decayer.get()==nullptr) return 1.0;
      return _decayer->Generate(P4(),_model->Products());
    }

    const DecayVectors* Decayer() const {return _decayer.get();}
    void SetDecayer(DecayVectors* decayer){_decayer.reset(decayer);}
    
    Double_t  GenerateComponents(){
      if(_model==nullptr )return 1.0;
      
      Double_t weight=0;
      while((weight)==0){ //sample until physical
	weight=Generate();
	weight*=_model->Intensity();
      }
      // can add beam divergence etc. here
      // ...
      return weight;
    }
    DecayType IsDecay() const noexcept final {return DecayType::Production;}

    void PostInit(ReactionInfo* info) ;
    
    // const LorentzVector* GetInteracting4Vector() const {return _interactingParticle;}
    // const LorentzVector* GetNominal4Vector() const {return &_nominal;}

    //if decays return current decay P4, if not return my P4
    const LorentzVector& GetInteracting4Vector() const {
      return _idxInteract==-1 ? GetNominal4Vector() : _model->Product(_idxInteract)->P4();
    }
    //if decays return nominal decay p4, if not return my P4
    const LorentzVector& GetNominal4Vector() const {return _nominal;}
    
    Int_t GetInteractingPdg()const {return _interactingPdg;}
    
    /*void SetVertexXYZT(double x,double y,double z,double t){
      _productionVertex.SetXYZT(x,y,z,t);
      }*/

 
    void SetAngleThetaPhi(Double_t th,Double_t phi);
    void SetHorSize(Double_t s){_sizeHor=s;}
    void SetVerSize(Double_t s){_sizeVer=s;}
    void SetHorDivergence(Double_t s){_divHor=s;}
    void SetVerDivergence(Double_t s){_divVer=s;}

  
  private:

    decaymodel_ptr _model={nullptr}; //not owner
    
    decayer_ptr   _decayer={nullptr}; //owner

    //  LorentzVector *_interactingParticle={nullptr}; //not owner;
    
    LorentzVector _nominal; //nominal beam 4-momentum
    // LorentzVector _productionVertex;
    //int _prodVertexID={0};
    int _interactingPdg={0};
    Double_t _dirTheta={0};
    Double_t _dirPhi={0};
    Double_t _sizeHor={0};
    Double_t _sizeVer={0};
    Double_t _divHor={0};
    Double_t _divVer={0};

    Int_t _idxInteract={-1}; //index of interacting particle in model, -1 if no model
    
    ClassDefOverride(elSpectro::CollidingParticle,1); //class CollidingParticle
    
  };//class CollidingParticle

}//namespace elSpectro
