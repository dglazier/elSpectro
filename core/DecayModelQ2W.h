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
#include "DistTH1.h"
#include "DecayModelst.h"
#include <TH1D.h>

namespace elSpectro{

  static TH1D HistFromLargestBinContents(const TH1D& h1,const TH1D& h2);
  
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

    const Particle* GetDecayBaryon()  noexcept{
      return dynamic_cast<const DecayModelst*>(GetGammaN()->Model())->GetBaryon();
    }
    const Particle* GetDecayMeson()  noexcept{
      return dynamic_cast<const DecayModelst*>(GetGammaN()->Model())->GetMeson();
    }

    double getQ2() const noexcept{return -_gamma.M2();}
    double getW() const noexcept{ return GetGammaN()->Mass();}
    double getThreshold() const noexcept{return _threshold;}
    void setThreshold(double val) noexcept{
      if(val<_threshold) return;
      _threshold=val;
    }
    
    void FindExcitationSpectra();
    
    //Q2 dependence from The H1 Collaboration: Elastic electroproduction of ρ mesons at HERA eqn 49 https://link.springer.com/content/pdf/10.1007/s100520000150.pdf
    constexpr double Q2H1RhoAt0() const  noexcept {return 3.0610097;}//1./TMath::Power((0.77549000*0.77549000),2.2); 
    double Q2H1Rho() const noexcept {return 1./TMath::Power((getQ2() + 0.77549000*0.77549000),2.2)/Q2H1RhoAt0(); }


    double dsigma() const override { return  dynamic_cast<DecayingParticle*>(_gstarNuc)->Model()->dsigma();}// * Q2 factor }

  protected:
    
    //mutable PhotoProdInfo _myInfo;//!
    const ReactionElectroProd* ProdInfo() const noexcept {return _prodInfo;}
    mutable LorentzVector _gamma;
    mutable PhotonPolarisationVector _photonPol;

  public:
    
    double PhaseSpaceFactorToQ2eq0(double W, double targetM) const noexcept {
      auto cmBoost=_gstarNuc->P4().BoostToCM();
      auto p1cm=boost(_gamma,cmBoost);
      return p1cm.P()/kine::PDK(W,0,targetM);
    }
 
  private:
    
    void Init();
    
    double _threshold = {0};
    
    Particle* _electron={nullptr};
    Particle* _gstarNuc={nullptr};

    ReactionElectroProd* _prodInfo={nullptr};

    TH1D _hWPhaseSpace;
    std::unique_ptr<DistTH1> _Wrealphoto_Dist;

    ClassDefOverride(elSpectro::DecayModelQ2W,1); //class DecayModelQ2W
    
  };//class DecayModelQ2W

}//namespace elSpectro
