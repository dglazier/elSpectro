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
#include "DistTF1.h"

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

    double MaximumMassPossible() const  noexcept override {

      Double_t maxMass=0;
      if(MassDistribution()!=nullptr){
	maxMass=MassDistribution()->GetMaxX();
      }
      else if(Pdg()!=-2211)
	maxMass = Particle::MaximumMassPossible();

      return maxMass;
    }
    
    double MinimumMassPossible() const  noexcept override {
      if(_minMass) return _minMass;
      
      auto minMass= _decay->MinimumMassPossible();
      if(MassDistribution()!=nullptr){
	if(MassDistribution()->GetMinX() > minMass)
	  minMass=MassDistribution()->GetMinX();
      }
      else if(Pdg()!=-2211)
	minMass = PdgMass();
      
      //std::cout<<"min masss "<<Pdg()<<" "<<minmass<<std::endl;
      return minMass;
    }
    
    void SetMinMass(double mass){_minMass=mass;}
    
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
    virtual void GenerateVertexPosition()  noexcept;
    
    const LorentzVector& DecayVertexPosition()const noexcept{return _decayVertex;}
    int DecayVertexID()const noexcept{return _decayVertexID;}

    void DetermineProductMasses(){ //only want to call intially in Process
      _decay->DetermineProductMasses();
    }
  
    DecayType IsDecay() const noexcept override {return _decayType;}

    void SetVertexXYZT(double x,double y,double z,double t){
      _decayVertex.SetXYZT(x,y,z,t);
    }
    
  protected:
    
    DecayVectors* mutableDecayer() const {return _decayer.get();}


  private:
     
    DecayModel* _decay={nullptr}; //not owner
    
    std::unique_ptr<DecayVectors> _decayer={nullptr}; //owner
    DistTF1* _decVertexDist=nullptr;//! needed if detached vertex
    
    double _minMass={0};
    LorentzVector _decayVertex;
    int _decayVertexID={0};
    DecayType _decayType;

    long _generateCalls={0};
    
    ClassDefOverride(elSpectro::DecayingParticle,1); //class DecayingParticle
    
  };//class DecayingParticle

}//namespace elSpectro
