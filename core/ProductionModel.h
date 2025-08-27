//////////////////////////////////////////////////////////////
///
///Class:	ProductionModel
///Description:
///            Responsible for returning production cross sesction
///            as function of W. Derived classes need to integrate
///            over other variables
#pragma once

#include "DecayModel.h"
#include "DecayingParticle.h"
#include "NBodyPhaseSpace.h"
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
    
    double get_W_FromParent()const {return Parent()->P4().M();}

    //Any preliminaries required
    bool ReadyForDecay() override{
      auto currW = get_W_FromParent();
      if(_cacheW==currW) return true;
      ChooseDecay();//Need full decay chain
      // auto idec = _channels.ChooseDecay(P4().M());
      //std::cout<<" ProductionModel::ReadyForDecay ChooseDecay()  "<<Pdg()<<" "<<idec<<std::endl;
 
      _cacheW = currW;
      SampleNBodyPhaseSpace(currW,this);
      //Need to check whether we should check threshold here
      //The W has been established at this point
      //so we should be above threshold anyway.
      //And not need to regenerate

      
      
      return true;
      //return CheckThreshold();
    }
 
  protected:
    Particle* GetMutableMeson() noexcept{return Products()[_idxMeson] ; }
    Particle* GetMutableBaryon() noexcept{return Products()[_idxBaryon]; }

  private:
    double _cacheW=-1.;
    uint _idxMeson=0;
    uint _idxBaryon=0;
  };

}
