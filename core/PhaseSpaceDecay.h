//////////////////////////////////////////////////////////////
///
///Class:		PhaseSpaceDecay
///Description:
///             Calculate intensity for phase space decay (=1 )!

#pragma once

#include "DecayModel.h"

namespace elSpectro{

 
  class PhaseSpaceDecay : public DecayModel {

  public:
    
    PhaseSpaceDecay()=default;
    //only declaring default constructor
    //so other 5 constructors also defaulted(rule of 5)
    //constructor to decay into particles
    //    PhaseSpaceDecay( particle_ptrs , const std::vector<int> pdgs );
    PhaseSpaceDecay(const decaying_objs& decs, const particle_objs& stables);
    // Each model must define its intensity
    // Phase space intensity is handled by MassPhaseSpace
    double Intensity() const final{
      if(CheckThreshold())
	return 1.;
      else
	return 0.;
    }
    bool RegenerateOnFail() const  noexcept final {return false;}
    void SetParent(DecayingParticle* pa) override;
    void PostInit(ReactionInfo* info)  override;
    void SetMassMaster(DistFlatMassMaster* master);
    
   void Print() const override;
 private:
     
    void nBodyDecayer(DecayingParticle* parent, particle_objs& stable, decaying_objs& unstable );
    void SetParentAndProducts(DecayingParticle* pa, particle_objs& stable,  decaying_objs& unstable);
 
 private:

    DistFlatMassMaster* _massMaster = nullptr;
    
    ClassDefOverride(elSpectro::PhaseSpaceDecay,1); //class PhaseSpaceDecay
    
  };//class PhaseSpaceDecay



}//namespace elSpectro
