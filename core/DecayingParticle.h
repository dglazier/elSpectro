//////////////////////////////////////////////////////////////
///
///Class:		DecayingParticle
///Description:
///             1) take Initial e- scatter
///             2) Initiate decay chain


#pragma once

#include "CurrentEventInfo.h"
#include "Particle.h"
#include "DecayModel.h"
#include "DecayVectors.h"
#include "TwoBodyFlat.h"
#include "ReactionInfo.h"

namespace elSpectro{
  
  enum class DecayStatus{ Decayed, TryAnother, ReGenerate };

  
  class DecayingParticle : public Particle {

  public:
    
    DecayingParticle()=delete;
    //cannot default construct need
    DecayingParticle(DecayModel* const model);
    //or
    DecayingParticle(int pdg,DecayModel* model,DecayVectors* decayer = new TwoBodyFlat());
    //or if decayer already give required distribution
    DecayingParticle(int pdg,DecayVectors* decayer,DecayModel* model);

    DecayModel*  Model()const {return _decay;}

    double Decay(){return _decayer->Generate(P4(),_decay->Products());}

    const DecayVectors* Decayer() const {return _decayer.get();}
    void SetDecayer(DecayVectors* decayer){_decayer.reset(decayer);}
    //  protected:
    
    virtual DecayStatus GenerateProducts();
    
    // virtual const CurrentEventInfo* EventInfo() const {return nullptr;}

    double MinimumMassPossible() const  noexcept override {
      
      double minmass= _decay->MinimumMassPossible();
      if(MassDistribution()!=nullptr){
	if(MassDistribution()->GetMinX() > minmass)
	  minmass=MassDistribution()->GetMinX();
      }
      else if(Pdg()!=-2211)
	minmass = PdgMass();
      
      //std::cout<<"min masss "<<Pdg()<<" "<<minmass<<std::endl;
      return minmass;
    };
    void TakeMinimumMass(){
      SetP4M( MinimumMassPossible() );
    }
    void TakePdgMass(){
      SetP4M( PdgMass() );
    }
    void Print() const override;


    double  PhaseSpaceWeightSq(){
      return _decay->PhaseSpaceWeightSq(Mass());
    }
    virtual void PostInit(ReactionInfo* info);

    //temporary until deal with vertices properly i.e. non zero
    void GenerateVertexPosition()  noexcept{
      //auto old=VertexPosition();
      if( IsDecay()==DecayType::Detached){

      }
      // _decayVertex.SetXYZT(old.X(),old.Y(),old.Z(),old.T());
    }
    const LorentzVector& DecayVertexPosition()const noexcept{return _decayVertex;}
    int DecayVertexID()const noexcept{return _decayVertexID;}

    void DetermineProductMasses(){ //only want to call intially in Process
      _decay->DetermineProductMasses();
    }
  
    DecayType IsDecay() const noexcept override {return _decayType;}
    
  protected:
    
    DecayVectors* mutableDecayer() const {return _decayer.get();}


  private:
     
    DecayModel* _decay={nullptr}; //not owner
    
    std::unique_ptr<DecayVectors> _decayer={nullptr}; //owner
    
    LorentzVector _decayVertex;
    int _decayVertexID={0};
    DecayType _decayType;
    
    ClassDefOverride(elSpectro::DecayingParticle,1); //class DecayingParticle
    
  };//class DecayingParticle

}//namespace elSpectro
