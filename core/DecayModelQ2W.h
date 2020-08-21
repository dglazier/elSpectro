//////////////////////////////////////////////////////////////
///
///Class:		DecayModelQ2W
///Description:
///             Control behaviour of Particle decay to Particle products
///             Defined by
///             1) list of Particle products
///             2) Intensity function dependent on Q2 and W
///
///            Note derived classes should include a constructor to initialise
///            DecayModelQ2W( particle_ptrs , const std::vector<int> pdgs );
#pragma once

#include "DecayModel.h"
#include "DecayVectors.h"
#include "DecayingParticle.h"
#include "ReactionInfo.h"
#include "PhotonPolarisationVector.h"

namespace elSpectro{

 
  class DecayModelQ2W : public DecayModel {

  public:
    
    DecayModelQ2W()=delete;
    //delete default constructor so have to use threshold
    //so other 5 constructors also defaulted(rule of 5)

    //constructor giving W theshold, just produces scatted electron kinematics
    DecayModelQ2W(  double thresh  );
    //constructor giving W theshold and subsequent primary decay of Nucl+gamma* system
    //DecayModelQ2W(  double thresh, DecayModel* gNmodel,DecayVectors* gNdecayer=nullptr);
    DecayModelQ2W(  double thresh, DecayModel* gNmodel,DecayVectors* gNdecayer=new TwoBodyFlat());
    
    // Each model must define its intensity
    double Intensity() const override;

    void PostInit(ReactionInfo* info) override;

    bool RegenerateOnFail() const  noexcept override {return false;}


    const Particle*  GetScatteredElectron() const noexcept{return _electron; }
    const DecayingParticle* GetGammaN() const noexcept{return dynamic_cast<const DecayingParticle*>(_gstarNuc); }
    
  protected:
    
    //mutable PhotoProdInfo _myInfo;//!
    const ReactionElectroProd* ProdInfo() const noexcept {return _prodInfo;}
    mutable LorentzVector _gamma;
    mutable PhotonPolarisationVector _photonPol;


 
  private:
    void Init();
    
    double _threshold = {0};
    
    Particle* _electron={nullptr};
    Particle* _gstarNuc={nullptr};

    ReactionElectroProd* _prodInfo={nullptr};
    
    ClassDefOverride(elSpectro::DecayModelQ2W,1); //class DecayModelQ2W
    
  };//class DecayModelQ2W

}//namespace elSpectro
