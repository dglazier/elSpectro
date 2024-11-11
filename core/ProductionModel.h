//////////////////////////////////////////////////////////////
///
///Class:	ProductionModel
///Description:
///            Responsible for returning production cross sesction
///            as function of W. Derived classes need to integrate
///            over other variables
#pragma once

#include "DecayModel.h"
#include <TH1D.h>


namespace elSpectro{

 
  class ProductionModel : public DecayModel {

  public :
    
    ProductionModel()=default;
    ProductionModel( const decaying_objs& decs, const particle_objs& stables) :
      DecayModel(decs,stables) {}

    virtual  TH1D   CrossSectionW(TH1D hist) = 0;

    // const Particle* GetDecayBaryon()  = 0;
    // const Particle* GetDecayMeson() = 0;
    const Particle* GetMeson() const noexcept{return Products()[_idxMeson] ; }
    const Particle* GetBaryon() const noexcept{return Products()[_idxBaryon]; }

    void SetMesonIdx(uint idx){_idxMeson=idx;}
    void SetBaryonIdx(uint idx){_idxBaryon=idx;}

  protected:
    Particle* GetMutableMeson() noexcept{return Products()[_idxMeson] ; }
    Particle* GetMutableBaryon() noexcept{return Products()[_idxBaryon]; }

  private:
    uint _idxMeson=0;
    uint _idxBaryon=0;
  };

}
