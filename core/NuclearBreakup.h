//////////////////////////////////////////////////////////////
///
///Class:		NuclearBreakup
///Description:
///          Break nucleus up into 2 pieces one of which used in collision

#pragma once

#include "DecayModel.h"

namespace elSpectro{

 
  class NuclearBreakup : public DecayModel {

  public:
    
    NuclearBreakup()=default;
    //only declaring default constructor
    //so other 5 constructors also defaulted(rule of 5)
    //constructor to decay into particles
    //NuclearBreakup( particle_ptrs , const std::vector<int> pdgs );
    NuclearBreakup( int pdg1, int pdg2 );

    // Each model must define its intensity
    // Phase space intensity is handled by MassPhaseSpace
    double Intensity() const final{
      return 1.;
    }
    bool RegenerateOnFail() const  noexcept final {return false;}
    ///void SetParent(DecayingParticle* pa);
    //void PostInit(ReactionInfo* info) override;
    
  private:
     ClassDefOverride(elSpectro::NuclearBreakup,1); //class NuclearBreakup
    
  };//class NuclearBreakup



}//namespace elSpectro
