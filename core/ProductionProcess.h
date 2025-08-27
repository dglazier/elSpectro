//////////////////////////////////////////////////////////////
///
///Class:		ProductionProcess
///Description:
///            Interface to different production processes
///            e.g Electroproduction, Photoproduction
///            1) model required to generate
///               reaction CoM state(from DecayingParticle)
#pragma once

#include "CurrentEventInfo.h"
#include "DecayModel.h"
//#include "Manager.h"
#include "DecayingParticle.h"
#include "TwoBodyProduction.h"
#include "CollidingParticle.h"
#include "Distribution.h"
#include "DistConst.h"

namespace elSpectro{

  class Manager;
  
  class ProductionProcess : public DecayingParticle {

     
  public:
 
    ProductionProcess(const CollidingParticle& p1,const CollidingParticle& p2,decaymodel_ptr  model,decayer_ptr  decayer=nullptr);
  
    virtual ~ProductionProcess()=default;
    ProductionProcess(const ProductionProcess& other); //need the virtual destructor...so rule of 5
    ProductionProcess(ProductionProcess&&)=default;
    ProductionProcess& operator=(const ProductionProcess& other);
    ProductionProcess& operator=(ProductionProcess&& other) = default;

    virtual const ReactionInfo* GetReactionInfo() const=0;
    
    virtual void InitGen() =0;
    //the intermediate state produced by this production process
    //e.g. g* + N for electron scattering
    virtual const DecayingParticle& Product() =0;

    
    virtual double IntegrateCrossSection(TwoBodyProduction* model2body) = 0;
    virtual double IntegrateCrossSectionFast(TwoBodyProduction* model2body) = 0;
    double IntegrateCrossSections();

    void PostInit(ReactionInfo* info) override;
   
    void SetCombinedBranchingFraction(double branch){_branchFrac=branch;}
    double BranchingFraction()const noexcept {return _branchFrac;}
    
    const particle_ptrs InitialParticles()const {return _initialParticles;}
    const particle_ptrs FinalParticles()const {return _finalParticles;}
    
    void AddInitialParticlePtr(Particle* p){
      _initialParticles.push_back(p);
    }

    void SetFinalParticles(const particle_ptrs& parts){
      _finalParticles.clear();
      _finalParticles = parts;
    }
     
    DecayType IsDecay() const noexcept override {return DecayType::Production;}

    virtual double dsigma() const{return 1;}

    virtual void GenerateVertexPosition(const ProductionProcess* production)  noexcept override{
      SetDecayVertexXYZT(_xvertexDist->SampleSingle(),
		    _yvertexDist->SampleSingle(),
		    _zvertexDist->SampleSingle(),
		    _tvertexDist->SampleSingle());
    }

    void GiveXVertexDist(Distribution* dist){_xvertexDist.reset(dist);}
    void GiveYVertexDist(Distribution* dist){_yvertexDist.reset(dist);}
    void GiveZVertexDist(Distribution* dist){_zvertexDist.reset(dist);}
    void GiveTVertexDist(Distribution* dist){_tvertexDist.reset(dist);}

 
    void BoostToLab(LorentzVector& boostme) const noexcept{
      boostme=ROOT::Math::VectorUtil::boost(boostme,_boostToLab);
    }
    void SetBoostToLab(const elSpectro::BetaVector& boostv){
      _boostToLab=boostv;
    }
    const elSpectro::BetaVector& GetBoostToLab() const noexcept{
      return _boostToLab;
    }
  protected:

    CollidingParticle* Incident1() {return &_in1;}
    CollidingParticle* Incident2() {return &_in2;}

     void InitVertex();
 
  private:
    ProductionProcess()=delete;
    particle_ptrs _initialParticles;
    particle_ptrs _finalParticles;
    
    dist_uptr _xvertexDist=dist_uptr{new DistConst{0}};
    dist_uptr _yvertexDist=dist_uptr{new DistConst{0}};
    dist_uptr _zvertexDist=dist_uptr{new DistConst{0}};
    dist_uptr _tvertexDist=dist_uptr{new DistConst{0}};
    
    double _branchFrac={1};

    CollidingParticle _in1;
    CollidingParticle _in2;
    
    elSpectro::BetaVector _boostToLab;

  };



}

