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
#include "DecayChannel.h"
#include "TwoBodyFlat.h"
#include "ReactionInfo.h"
#include "DistTF1.h"

namespace elSpectro{
  
  enum class DecayStatus{ Decayed, TryAnother, ReGenerate };

  class DecayChannel;
  class ProductionProcess;
  
  class DecayingParticle : public Particle {

  public:
    
    DecayingParticle()=default;
    //cannot default construct need
    DecayingParticle(decaymodel_ptr model);
    //or
    DecayingParticle(int pdg,decaymodel_ptr model,decayer_ptr decayer);
    //or if decayer already give required distribution
    DecayingParticle(int pdg,decayer_ptr decayer,decaymodel_ptr model);
   //or if wish to add many decays
    DecayingParticle(int pdg):Particle{pdg}{};

    virtual ~DecayingParticle()=default;
    DecayingParticle(const DecayingParticle& other);
    DecayingParticle(DecayingParticle&& other);
    DecayingParticle& operator=(const DecayingParticle& other);
    DecayingParticle& operator=(DecayingParticle&& other);
    void CopyOther(const DecayingParticle& other);

    
    DecayModel*  Model() const {return _channels.CurrModel();}

    double Decay(){return Decayer()->Generate(P4(),Model()->Products());}

    DecayVectors* Decayer() const {return _channels.CurrDecayer();}
    
    void SetDecayer0(decayer_ptr decayer){
      _channels.SetDecayer(0,std::move(decayer));
      //_decayer = _channels.CurrDecayer();
    }
    void SetDecayer(decayer_ptr decayer){
      _channels.SetDecayer(_channels.CurrChannel(),std::move(decayer));
      //_decayer = _channels.CurrDecayer();
    }
    
    virtual DecayStatus GenerateProducts(const ProductionProcess* production);
    
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
      std::cout<<"DecayingParticle min masss "<<Pdg()<<" "<<Model()<<" "<<MassDistribution()<<std::endl;
      auto minMass= Model()->MinimumMassPossible();
     std::cout<<"DecayingParticle min masss "<<minMass<<std::endl;
      if(MassDistribution()!=nullptr){
	if(MassDistribution()->GetMinX() > minMass)
	  minMass=MassDistribution()->GetMinX();
      }
      else if(Pdg()!=-2211){
	minMass = PdgMass();
      }
      //std::cout<<"min masss "<<Pdg()<<" "<<minmass<<std::endl;
      return _minMass=minMass;
    }

    double IntegratedMass(double W, double M_other) const{
      //need to call this in TwoBodyProduction when full mass range
      //is not accessible due to low W high mass.
      //should weight the cross section by the ratio of
      //full integrated mass to this integrated mass
      double maxMass = W-M_other; //total invariant mass - mass of other
      //If no distribution return 1 if above threshold, 0 if below
      if(MassDistribution()==nullptr){
	return maxMass > MinimumMassPossible() ? 1 : 0 ;
      }
      
      //now need to integrate mass distribution from minimum to max
     return  MassDistribution()->Integrate1DX(_minMass,maxMass);
      
    }

    double MeanMass(double W, double M_other) const{
      //need to call this in TwoBodyProduction when full mass range
      //is not accessible due to low W high mass.
      double maxMass = W-M_other; //total invariant mass - mass of other
      //If no distribution return pdg 
      if(MassDistribution()==nullptr){
	return  PdgMass();
      }
      
      //now need to integrate mass distribution from minimum to max
     return  MassDistribution()->Mean1DX(_minMass,maxMass);
      
    }

    
    void SetMinMass(double mass) const {_minMass=mass;}
    
    void TakeMaximumMass(){
      SetP4M( MaximumMassPossible() );
    }
   void TakeMinimumMass(){
      SetP4M( MinimumMassPossible() );
    }
    void TakePdgMass(){
      SetP4M( PdgMass() );
    }
    void Print() const override;


    double  PhaseSpaceWeightSq(){
      return Model()->PhaseSpaceWeightSq(Mass());
    }
    virtual void PostInit(ReactionInfo* info);

    //temporary until deal with vertices properly i.e. non zero
    virtual void GenerateVertexPosition(const ProductionProcess* production)  noexcept;
    
    const LorentzVector& DecayVertexPosition()const noexcept{return _decayVertex;}
    int DecayVertexID()const noexcept{return _decayVertexID;}

    void DetermineProductMasses(){ //only want to call intially in Process
      Model()->DetermineProductMasses();
    }
  
    DecayType IsDecay() const noexcept override {return _decayType;}

    void SetVertexXYZT(double x,double y,double z,double t){
      _decayVertex.SetXYZT(x,y,z,t);
    }

    void AddDecay(double bratio,decaymodel_ptr  mod,decayer_ptr  dec){
      std::cout<<"DecayingParticle  AddDecay : "<<bratio<< std::endl;
      _channels.AddDecay(this,bratio,std::move(mod),std::move(dec));
      std::cout<<"DecayingParticle  AddDecay : "<<_channels.CurrModel()<< std::endl;
      // _channels.CurrModel()->SetParent(this);

    }
    void ChooseDecay(){
      //Randomly select a decay channel based on branching ratio
      _channels.ChooseDecay();
       //assign model and decayer for this event
      // _decay = _channels.CurrModel();
      // _decayer = _channels.CurrDecayer();
      //recurse daughter particles
      Model()->ChooseDecay();
    }

    void EventParticles(particle_ptrs& parts){
      Model()->EventParticles(parts);
    }

    const DecayChannel& Channels() const{return _channels;}
    
  protected:
    
    DecayVectors* mutableDecayer() const {return _channels.CurrDecayer();}


  private:

    DecayChannel _channels;
    // mutable DecayModel* _decay={nullptr}; //not owner
    //mutable DecayVectors* _decayer={nullptr}; //owner
    //ProductionProcess* _process={nullptr};//not owner
    
    DistTF1* _decVertexDist=nullptr;//! needed if detached vertex
    
    mutable double _minMass={0};
    LorentzVector _decayVertex;
    int _decayVertexID={0};
    DecayType _decayType;

    long _generateCalls={0};
    size_t _gtSample={0};
    
    ClassDefOverride(elSpectro::DecayingParticle,1); //class DecayingParticle
    
  };//class DecayingParticle


}//namespace elSpectro
