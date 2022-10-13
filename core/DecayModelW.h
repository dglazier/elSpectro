//////////////////////////////////////////////////////////////
///
///Class:		DecayModelW
///Description:
///             Control behaviour of Particle decay to Particle products
///             Defined by
///             1) list of Particle products
///             2) Intensity function dependent W
///
///            Note derived classes should include a constructor to initialise
///            DecayModelW( particle_ptrs , const std::vector<int> pdgs );
#pragma once

#include "DecayModel.h"
#include "DecayVectors.h"
#include "DecayingParticle.h"
#include "ReactionInfo.h"
#include "PhotonPolarisationVector.h"
#include "DistTH1.h"
#include "DecayModelst.h"
#include <TH1D.h>

namespace elSpectro{

  static TH1D HistFromLargestBinContents(const TH1D& h1,const TH1D& h2);
  
  class DecayModelW : public DecayModel {

  public:
    
    DecayModelW()=delete;
    //delete default constructor so have to use threshold
    //so other 5 constructors also defaulted(rule of 5)

    //constructor giving W theshold, just produces scatted electron kinematics
    DecayModelW(  double thresh  );
    //constructor giving W theshold and subsequent primary decay of Nucl+gamma* system
    //DecayModelW(  double thresh, DecayModel* gNmodel,DecayVectors* gNdecayer=nullptr);
    DecayModelW(  double thresh, DecayModel* gNmodel,DecayVectors* gNdecayer=new TwoBodyFlat());
    
    // Each model must define its intensity
    double Intensity() const override;

    void PostInit(ReactionInfo* info) override;

    bool RegenerateOnFail() const  noexcept override {return false;}


    // const Particle*  GetPhoton() const noexcept{return _photon; }
    const DecayingParticle* GetGammaN() const noexcept{return dynamic_cast<const DecayingParticle*>(_gstarNuc); }

    const Particle* GetDecayBaryon()  noexcept{
      return dynamic_cast<const DecayModelst*>(GetGammaN()->Model())->GetBaryon();
    }
    const Particle* GetDecayMeson()  noexcept{
      return dynamic_cast<const DecayModelst*>(GetGammaN()->Model())->GetMeson();
    }

    double getW() const noexcept{ return _gstarNuc->Mass();}
    double getThreshold() const noexcept{return _threshold;}
    void setThreshold(double val) noexcept{
      if(val<_threshold) return;
      _threshold=val;
    }
    
    void FindExcitationSpectra();
    //DistTH1* GetApproxWDist() const {return _Wrealphoto_Dist.get();}
    Distribution* GetApproxWDist() const {return _Wrealphoto_Dist.get();}
 
    double dsigma() const override { return  dynamic_cast<DecayingParticle*>(_gstarNuc)->Model()->dsigma();}

  protected:
    
    //mutable PhotoProdInfo _myInfo;//!
    const ReactionPhotoProd* ProdInfo() const noexcept {return _prodInfo;}
    mutable LorentzVector _gamma;
    mutable PhotonPolarisationVector _photonPol;

  public:
    
     void ZeroPhoton(){
      _gamma.SetXYZT(0,0,0,0);
    }
    
  private:
    
    void Init();
    
    double _threshold = {0};
    
    LorentzVector* _photon={nullptr};
    Particle* _gstarNuc={nullptr};

    ReactionPhotoProd* _prodInfo={nullptr};

    TH1D _hWPhaseSpace;
    //  std::unique_ptr<DistTH1> _Wrealphoto_Dist;
    std::unique_ptr<Distribution> _Wrealphoto_Dist;

    ClassDefOverride(elSpectro::DecayModelW,1); //class DecayModelW
    
  };//class DecayModelW

}//namespace elSpectro
