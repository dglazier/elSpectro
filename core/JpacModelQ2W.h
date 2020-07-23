//////////////////////////////////////////////////////////////
///
///Class:		JpacModelQ2W
///Description:
///             Control behaviour of Particle decay to Particle products
///             Defined by
///             1) preconfigured jpacPhoto amplitude
///             2) it decay as a function of s and t JpacDecayst
///
///            Note derived classes should include a constructor to initialise
///            JpacModelQ2W( particle_ptrs , const std::vector<int> pdgs );
#pragma once

#include "DecayModelQ2W.h"
#include "JpacModelst.h"
#include "DistTH1.h"
#include "amplitudes/amplitude.hpp"

namespace elSpectro{
  
  using jpacAmp_ptr = jpacPhoto::amplitude*;
  
  class JpacModelQ2W : public DecayModelQ2W {
    
  public:
    
    JpacModelQ2W()=delete;
    //constructor giving jpac amplitude pointer (which we will now own)
    //and decay particles 
    JpacModelQ2W( jpacPhoto::amplitude* amp, particle_ptrs parts,
		  const std::vector<int> pdgs  );
    //constructor giving strengths of s and t distribution envelope
    JpacModelQ2W( jpacPhoto::amplitude* amp, particle_ptrs parts,
		  const std::vector<int> pdgs,
		  double s_strength,
		  double t_strength,double t_slope);
    
  
    double Intensity() const override;
    
    void PostInit(ReactionInfo* info) override;

    //regenerate => draw another photon candidate
    bool RegenerateOnFail() const noexcept final {return true;}
    
    const Particle* GetDecayBaryon() const noexcept{
      return dynamic_cast<const JpacModelst*>(GetGammaN()->Model())->GetBaryon();
    }
    const Particle* GetDecayMeson() const noexcept{
      return dynamic_cast<const JpacModelst*>(GetGammaN()->Model())->GetMeson();
    }
 
  private:
    
    jpacAmp_ptr _amp={nullptr}; //I am not the owner

    std::unique_ptr<DistTH1> _W_Dist;
    
    ClassDefOverride(elSpectro::JpacModelQ2W,1); //class JpacModelQ2W
    
  };//class JpacModelQ2W

}//namespace elSpectro
