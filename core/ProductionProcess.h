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

namespace elSpectro{
  
  
  class ProductionProcess : public DecayingParticle {

     
  public:
    //virtual ~ProductionProcess()=default;

    //only construct via and take ownership of model 
    ProductionProcess(DecayModel* model);
    ProductionProcess(int pdg, DecayVectors* decayer, DecayModel* model);

    virtual ~ProductionProcess()=default;
    ProductionProcess(const ProductionProcess& other); //need the virtual destructor...so rule of 5
    ProductionProcess(ProductionProcess&&)=default;
    ProductionProcess& operator=(const ProductionProcess& other);
    ProductionProcess& operator=(ProductionProcess&& other) = default;

 
    
    virtual void InitGen() =0;

   
    
    const particle_constptrs InitialParticles()const {return _initialParticles;}
    void AddInitialParticlePtr(const Particle* p){
      _initialParticles.push_back(p);
    }
    
    DecayType IsDecay() const noexcept override {return DecayType::Production;}

  protected:

    void DetermineAllMasses(){
      DetermineProductMasses();
    }
    
  private:
    ProductionProcess()=default;
    
    particle_constptrs _initialParticles;
    
  };



}

