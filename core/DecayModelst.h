//////////////////////////////////////////////////////////////
///
///Class:		DecayModelst
///Description:
///             Control behaviour of Particle decay to Particle products
///             Defined by
///             1) Two body phase space as a fucntion of s and t
///             
///
///            Note derived classes should include a constructor to initialise
///            DecayModelst( particle_ptrs , const std::vector<int> pdgs );
#pragma once

#include "DecayModel.h"
#include "SDME.h"
#include "FunctionsForElectronScattering.h"
#include "DecayingParticle.h"
#include <TH1D.h>

namespace elSpectro{

 
  class DecayModelst : public DecayModel {

  public:
    
    DecayModelst()=delete;
    //constructor giving jpac amplitude pointer (which we will now own)
    //and decay particles 
    DecayModelst(  particle_ptrs parts,const std::vector<int> pdgs  );
    
    //
    double Intensity() const override; //this should perhaps be final so derived classes cannot overwrite...
    
    void PostInit(ReactionInfo* info) override;

    bool RegenerateOnFail() const noexcept override {return true;};
    bool HasAngularDistribution() override{return false; }

    const Particle* GetMeson() const noexcept{return _meson; }
    const Particle* GetBaryon() const noexcept{return _baryon; }

    void SetUseSDME(bool use=true){_useSDME=use;}
    

    /* double PgammaCMsq()const noexcept{*/
    double PgammaCMsq() const noexcept{
      if(_photon->M()==0) return kine::PDK2(_W,0,_target->M());
      auto  pgammaCM= PgammaCM();
      return  pgammaCM* pgammaCM; //for dt phase space factor
    }
    
    double PgammaCM()const noexcept{
      //in case no photon 4-vector yet
      if(_photon->M()==0) return kine::PDK(_W,0,_target->M());
      //else PDK does not qork for virtual photons
      auto cmBoost=Parent()->P4().BoostToCM();
      auto p1cm=boost(*_photon,cmBoost);
      
      //     std::cout<<"PgammaCMsq M"<<_photon->M()<<" PLAB "<<_photon->P()<<" PCM "<<p1cm.P()<<" or "<<kine::PDK(_W,_photon->M(),_target->M())<<" or "<<kine::PDK(_W,_photon->M2(),_target->M())<<" or "<<kine::PDK(_W,0,_target->M())<<std::endl;
      return p1cm.P();
    }
    
    /* double PgammaCMsq() const noexcept{
      //in case no photon 4-vector yet
      if(_photon->M()==0) return kine::PDK2(_W,0,_target->M());
      //else PDK does not qork for virtual photons
      auto cmBoost=Parent()->P4().BoostToCM();
      auto p1cm=boost(*_photon,cmBoost);
      
      //     std::cout<<"PgammaCMsq M"<<_photon->M()<<" PLAB "<<_photon->P()<<" PCM "<<p1cm.P()<<" or "<<kine::PDK(_W,_photon->M(),_target->M())<<" or "<<kine::PDK(_W,_photon->M2(),_target->M())<<" or "<<kine::PDK(_W,0,_target->M())<<std::endl;
      return p1cm.P()*p1cm.P();
      }*/
    
    const ReactionElectroProd* ProductionInfo() const { return _prodInfo; }
    
    void HistIntegratedXSection_ds(TH1D& hist);
    void HistIntegratedXSection(TH1D& hist);
    void HistMaxXSection(TH1D& hist);
    
    double PhaseSpaceFactor() const noexcept {
      /* auto fluxPhaseSpace = p1*_W;//eqn 47.28b https://pdg.lbl.gov/2019/reviews/rpp2019-rev-kinematics.pdf
	 auto ans =  1./fluxPhaseSpace
	 * kine::PhaseSpaceFactorDt(_W,p1,_meson->Mass(),_baryon->Mass())
	 * kine::PDK(_W,_meson->Mass(),_baryon->Mass())/_W
	 * PhaseSpaceNorm();//nbarn it
	 */ //Note above full calculation simplifies to
      //return PhaseSpaceNorm()/_s/kine::PDK2(_W,_photon->M(),_target->M());
      //Please note kine::PDK2(_W,_photon->M(),_target->M()) does not give
      //correct momentum
      return PhaseSpaceNorm()/_s/PgammaCMsq();
      //this would not be the case if the differential was dcosth rather than t
    }
    
    double PhaseSpaceFactorCosTh() const noexcept {
      return PhaseSpaceNormCosTh()* kine::PDK(_W,_meson->Mass(),_baryon->Mass())/_s/PgammaCM();
    }
    
  protected:
    
    virtual double MatrixElementsSquared_L() const {return 0; }
    virtual double MatrixElementsSquared_T() const {return 1; } //just real photon by default
    
    constexpr double  PhaseSpaceNorm() const {return 1./(2.56819E-6)/64/TMath::Pi();}// Convert from GeV^-2 -> nb
 
    constexpr double  PhaseSpaceNormCosTh() const {return 1./(2.56819E-6)/32/TMath::Pi();}// Convert from GeV^-2 -> nb

    
    SDME* const  GetMesonSDMEs() const {return _sdmeMeson;}
    SDME* const  GetBaryonSDMEs() const {return _sdmeBaryon;}
    
    virtual void CalcMesonSDMEs() const {};
    virtual void CalcBaryonSDMEs() const {};

    double FindMaxOfIntensity();

  public:
    double get_s() const noexcept{ return _s; }
    double get_t() const noexcept { return _t; }
    double get_W() const noexcept { return _W; }
    double get_Q2() const noexcept { return -_photon->M2(); }

  
    double kinCM_MesonP(double W) const {
      //  std::cout<<"kinCM_MesonP "<< kine::PDK(W,_meson->P4().M(),_baryon->P4().M()) <<W<<" "<<_meson->P4().M()<<" "<<_baryon->P4().M()<<std::endl;
       return kine::PDK(W,_meson->P4().M(),_baryon->P4().M());

    }
    double kinCM_MesonE(double W) const {
       auto m2_a =_meson->P4().M2();
       auto m2_b =_baryon->P4().M2();
       // std::cout<<"kinCM_MesonE "<< (W*W + m2_a - m2_b)/(2.0*W)<<std::endl;
       return (W*W + m2_a - m2_b)/(2.0*W);
    }
    double kin_tFromWCosTh(double cosTh) const{
      double W = Parent()->P4().M();
      auto cmBoost=Parent()->P4().BoostToCM();
      auto p1cm=boost(*_photon,cmBoost);
      // std::cout<<"kin_tFromWCosTh "<<p1cm.M2() + _meson->M2() - 2 * (p1cm.E()* kinCM_MesonE(W)-p1cm.P()* kinCM_MesonP(W)*cosTh)<<std::endl;

      return p1cm.M2() + _meson->M2() - 2 * (p1cm.E()* kinCM_MesonE(W)-p1cm.P()* kinCM_MesonP(W)*cosTh);
    }

    double dsigma_costh(double cosTh){
      //For integrating cross section
      _W = Parent()->P4().M();
      _s=_W*_W;
      _t = kin_tFromWCosTh(cosTh);
      // std::cout<<"dsigma_costh t "<<_t<<" "<<PhaseSpaceFactorCosTh()<<" "<<MatrixElementsSquared_T()<<std::endl;
  
      return PhaseSpaceFactorCosTh()*(MatrixElementsSquared_T()+(_photonPol->Epsilon()+_photonPol->Delta())*MatrixElementsSquared_L()) ;//2pi=>integrated over phi
    }
    double dsigma_costhW(double cosTh,double W){
     //For integrating cross section
      _W = W;
      _s=_W*_W;
      _t = kin_tFromWCosTh(cosTh);
      return PhaseSpaceFactorCosTh()*(MatrixElementsSquared_T()+(_photonPol->Epsilon()+_photonPol->Delta())*MatrixElementsSquared_L()) ;//2pi=>integrated over phi
    }

  private:

    double DifferentialXSect() const{//dont let others call this as need _s, _W and _t set
      //Note if your derived model already gives differential cross section
      //you will need to divide by PhaseSpaceFactor to get MatrixElementSquared from it
      return _dsigma=PhaseSpaceFactor() * ( MatrixElementsSquared_T() +
          				    (_photonPol->Epsilon()+_photonPol->Delta())*MatrixElementsSquared_L()); //eqn from Seyboth and Wolf
    }
       
    SDME* _sdmeMeson={nullptr};
    SDME* _sdmeBaryon={nullptr};
    PhotonPolarisationVector* _photonPol={nullptr};
 

    ReactionElectroProd* _prodInfo={nullptr};
 
    Particle* _baryon={nullptr};
    Particle* _meson={nullptr};
    LorentzVector* _photon={nullptr};
    LorentzVector* _target={nullptr};//{0,0,0,escat::M_pr()};
    LorentzVector* _ebeam={nullptr};//{0,0,0,escat::M_pr()};
 
    mutable double _max={0};
    mutable double _s={0};
    mutable double _t={0};
    mutable double _W={0};
    mutable double _dt={0};
    mutable double _dsigma={0};
    double _Wmax={0};
 
    bool _useSDME={false};

    ClassDefOverride(elSpectro::DecayModelst,1); //class DecayModelst
    
  };//class DecayModelst

}//namespace elSpectro
