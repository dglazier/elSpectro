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
#include "ParticleManager.h"
#include "DecayingParticle.h"
#include "CollidingParticle.h"
#include "Distribution.h"
#include "DistConst.h"

namespace elSpectro{

  
  class ProductionProcess : public DecayingParticle {

     
  public:
    //virtual ~ProductionProcess()=default;

    //only construct via and take ownership of model 
    ProductionProcess(DecayModel* model);
    ProductionProcess(CollidingParticle* p1,CollidingParticle* p2,DecayModel* model);
    ProductionProcess(int pdg, DecayVectors* decayer, DecayModel* model);

    virtual ~ProductionProcess()=default;
    ProductionProcess(const ProductionProcess& other); //need the virtual destructor...so rule of 5
    ProductionProcess(ProductionProcess&&)=default;
    ProductionProcess& operator=(const ProductionProcess& other);
    ProductionProcess& operator=(ProductionProcess&& other) = default;

 
    
    virtual void InitGen() =0;
    
    void PostInit(ReactionInfo* info) override;
    void Init();
   
    virtual double IntegrateCrossSection() = 0;
    virtual double IntegrateCrossSectionFast() = 0;
    
    void SetCombinedBranchingFraction(double branch){_branchFrac=branch;}
    double BranchingFraction()const noexcept {return _branchFrac;}
    
    const particle_constptrs InitialParticles()const {return _initialParticles;}
    void AddInitialParticlePtr(const Particle* p){
      _initialParticles.push_back(p);
    }
    
    DecayType IsDecay() const noexcept override {return DecayType::Production;}

    virtual double dsigma() const{return 1;}

    virtual void GenerateVertexPosition()  noexcept override{
      SetVertexXYZT(_xvertexDist->SampleSingle(),
		    _yvertexDist->SampleSingle(),
		    _zvertexDist->SampleSingle(),
		    _tvertexDist->SampleSingle());
    }

    void GiveXVertexDist(Distribution* dist){_xvertexDist.reset(dist);}
    void GiveYVertexDist(Distribution* dist){_yvertexDist.reset(dist);}
    void GiveZVertexDist(Distribution* dist){_zvertexDist.reset(dist);}
    void GiveTVertexDist(Distribution* dist){_tvertexDist.reset(dist);}

    CollidingParticle* Incident1() const {return _in1;}
    CollidingParticle* Incident2() const {return _in2;}

    void BoostToLab(LorentzVector& boostme) const{
      boostme=ROOT::Math::VectorUtil::boost(boostme,_boostToLab);
    }
    void SetBoostToLab(const elSpectro::BetaVector& boostv){
      _boostToLab=boostv;
    }
  protected:

   
    
  private:
    ProductionProcess()=default;
    particle_constptrs _initialParticles;
    
    dist_uptr _xvertexDist=dist_uptr{new DistConst{0}};
    dist_uptr _yvertexDist=dist_uptr{new DistConst{0}};
    dist_uptr _zvertexDist=dist_uptr{new DistConst{0}};
    dist_uptr _tvertexDist=dist_uptr{new DistConst{0}};
    
    double _branchFrac={1};

    CollidingParticle* _in1={nullptr};
    CollidingParticle* _in2={nullptr};
    
    elSpectro::BetaVector _boostToLab;

  };



}

