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

namespace elSpectro{
    
 
  class DecayingParticle : public Particle {

  public:
    
    DecayingParticle()=delete;
    //cannot default construct need
    DecayingParticle(DecayModel* const model);
    //or
    DecayingParticle(int pdg,DecayModel* model,DecayVectors* decayer = new TwoBodyFlat());
    //or if decayer already give required distribution
    DecayingParticle(int pdg,DecayVectors* decayer,DecayModel* model);

    DecayModel* const Model()const {return _decay;}

    double Decay(){return _decayer->Generate(P4(),_decay->Products());}

    const DecayVectors* Decayer() const {return _decayer;} 
    //  protected:
    
    virtual double GenerateProducts(const CurrentEventInfo* parentInfo=nullptr);
    
    virtual const CurrentEventInfo* EventInfo() const {return nullptr;}

    double MinimumMassPossible() const override {
      double minmass= _decay->MinimumMassPossible();
      if(MassDistribution()!=nullptr)
	if(MassDistribution()->GetMinX() > minmass)
	  minmass=MassDistribution()->GetMinX();
      return minmass;
    };
    
    void Print() const override;


    double  PhaseSpaceWeightSq(){
      return _decay->PhaseSpaceWeightSq(Mass());
    }

   
  private:
     
    DecayModel* _decay={nullptr}; //not owner
    DecayVectors* _decayer={nullptr}; //not owner
    
 
    ClassDefOverride(elSpectro::DecayingParticle,1); //class DecayingParticle
    
  };//class DecayingParticle

}//namespace elSpectro
