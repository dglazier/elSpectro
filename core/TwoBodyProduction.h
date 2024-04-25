//////////////////////////////////////////////////////////////
///
///Class:		TwoBodyProduction
///Description:
///             Control behaviour of Particle decay to
///             Quasi-TwoBody meson + baryon final states
///             dependent on s and t
///             
#pragma once

#include "DecayModel.h"
#include "DecayingParticle.h"
#include "FunctionsForElectronScattering.h"
#include "TwoBodyEnvelope.h"
#include "DistTH1.h"
#include <TH1D.h>
#include <TH2D.h>

namespace elSpectro{

 
  class TwoBodyProduction : public DecayModel {

  public:
    
    TwoBodyProduction()=delete;
    //constructor giving jpac amplitude pointer (which we will now own)
    //and decay particles 
    TwoBodyProduction(  particle_ptrs parts,const std::vector<int> pdgs  );
    
 
    double Intensity() const override;
    
    void PostInit(ReactionInfo* info) override;

    bool RegenerateOnFail() const noexcept override {return true;};
    bool HasAngularDistribution() override{return false; }

    const Particle* GetMeson() const noexcept{return _meson; }
    const Particle* GetBaryon() const noexcept{return _baryon; }
    const LorentzVector* GetTarget() const noexcept{return _target; }
    const LorentzVector* GetPhoton() const noexcept{return _photon; }

    virtual double FindMaxOfIntensity();
    virtual void HistIntegratedXSection(TH1D& hist);
    virtual void HistMaxXSection(TH1D& hist);
    //std::unique_ptr<TwoBodyEnvelope> Build_WCosTh_Envelope(const LorentzVector *target, const LorentzVector *ebeam) const;
    void Use_WCosTh_Envelope();
    void Build_WCosTh_Envelope();
    
    virtual double dsigma_dcosth(double W,double cth) const;
    virtual double sigma(double W) const{ return 1.;}

    double GetUsedMassWeight() const {return _massWeight;}
    void SetMassWeight(double val) const{_massWeight=val;}
    
    double get_s() const noexcept{ return _s; }
    double get_t() const noexcept { return _t; }
    double get_W() const noexcept { return _W; }
    double get_cosThCM() const noexcept { return _cosThCM; }

    void set_W(double val) const{ _W=val;_s=_W*_W; }
    void set_s(double val) const{ _s=val;_W=TMath::Sqrt(_s);}
    void set_t(double val) const{ _t=val; }
    void set_cosThCM(double val) const{ _cosThCM=val; }
    
    double get_max() const noexcept{ return _max; }
    void set_max(double val) const { _max=val; }

    double get_W_FromParent()const {return Parent()->P4().M();}
  protected:
    
    virtual double PhaseSpaceFactor() const;
    virtual double MatrixElementsSquared_L() const {return 0; }
    virtual double MatrixElementsSquared_T() const {return 1; } //just real photon by default

    double Q2PhaseSpaceCorrect() const;
    constexpr double  PhaseSpaceNorm() const;
    double PgammaCMsq() const noexcept;
    double PgammaCM()const noexcept;
    double PhaseSpaceFactor_dCosTh() const noexcept ;
    double kinCM_MesonP(double W) const;
    double kinCM_MesonE(double W) const ;
    double kin_tFromWCosTh(double W, double cosTh) const;
    double CalcCosThCM() const;
    double kin_tFromWCosThatQ20(double W, double cosTh) const;
    
    double DiffXS() const{//dont let others call this as need _s, _W and _t set
      //dsigma/dcostheta

      // std::cout<<" DiffXS() "<<std::endl;
      //PhaseSpaceFactor();
      //std::cout<<" DiffXS() "<<std::endl;
      //MatrixElementsSquared_T();
      //std::cout<<" DiffXS() "<<std::endl;
      //Note if your derived model already gives differential cross section
      //you will need to divide by PhaseSpaceFactor to get MatrixElementSquared from it
      // std::cout<<"ps = "<<PhaseSpaceFactor() <<" mat = "<< MatrixElementsSquared_T()<<" mass "<<_meson->P4().M()<<std::endl;
      return _photonPol==nullptr ?
	PhaseSpaceFactor() * MatrixElementsSquared_T() :
	PhaseSpaceFactor() * ( MatrixElementsSquared_T()
	  + (_photonPol->Epsilon()+_photonPol->Delta())*MatrixElementsSquared_L());	
 

	// auto result = PhaseSpaceFactor() *
	// ( MatrixElementsSquared_T()
	//   + (_photonPol->Epsilon()+_photonPol->Delta())*MatrixElementsSquared_L());
      //eqn from Seyboth and Wolf
      
      //return result;
      
    }
    virtual void CalcKine() const{
      _W = Parent()->P4().M();
      _s=_W*_W;
      //calculate thetaCM for this event
      CalcCosThCM();
    }
    const ReactionPhotoProd* ProductionInfo() const {return _prodInfo;}

    bool IsSampling() const {return _isSampling;}
    
  protected:
    Particle* _baryon={nullptr};
    Particle* _meson={nullptr};
    LorentzVector* _photon={nullptr};
    LorentzVector* _target={nullptr};
    const LorentzVector* _ebeam={nullptr};
    PhotonPolarisationVector* _photonPol={nullptr};
    ReactionPhotoProd* _prodInfo={nullptr};
 
  private:
    //keep private so only called when s,W,cth,...are set
    //either through CalcKine or given values, e.g. dsigma_dcosth(cth)
    virtual double dsigma_dcosth() const{
     return  DiffXS() ;
    }

    TH1D _histHighXS;
 
    mutable double _max={0};
    mutable double _s={0};
    mutable double _t={0};
    mutable double _W={0};
    mutable double _cosThCM={0};

    mutable double _Ntries=0;
    mutable double _NhighWeight=0;
    mutable double _massWeight=1.;
    
    bool _isElProd={true};
    mutable bool _isSampling={false};
    bool _createMyEnvelope={false};

    std::unique_ptr<DistTH1> _distHighXS;
    std::unique_ptr<DistTH1>  _distMassFrac;
    ClassDefOverride(elSpectro::TwoBodyProduction,1); //class TwoBodyProduction
 
  };


  //Define inline functions
  inline double TwoBodyProduction::dsigma_dcosth(double W,double cth) const{
      //Parent()->P4().M();
    set_cosThCM(cth);
    set_W(W);
    
    // in case we need t make sure up-to-date
    set_t(kin_tFromWCosThatQ20(_W,_cosThCM));

    return dsigma_dcosth();
  }
  //Kinematic factors
  inline double TwoBodyProduction::Q2PhaseSpaceCorrect() const {
    return TMath::Sqrt(kine::PDK2(_W,0,_target->M())/PgammaCMsq());
  }

  inline constexpr double  TwoBodyProduction::PhaseSpaceNorm() const {return 1./(2.56819E-6)/32/TMath::Pi();}// Convert from GeV^-2 -> nb
     
     
  inline double TwoBodyProduction::PgammaCMsq() const noexcept{
    if(_photon->M()==0) return kine::PDK2(_W,0,_target->M());
    auto  pgammaCM= PgammaCM();
    return  pgammaCM* pgammaCM; //for dt phase space factor
  }
    
  inline double TwoBodyProduction::PgammaCM()const noexcept{
    //in case no photon 4-vector yet
    if(_photon==nullptr) return kine::PDK(_W,0,_target->M());
    if(_photon->M()==0) return kine::PDK(_W,0,_target->M());
    //else PDK does not qork for virtual photons
    auto cmBoost=Parent()->P4().BoostToCM();
    auto p1cm=boost(*_photon,cmBoost);
    //std::cout<<" TwoBodyProduction::PgammaCM() "<<1./p1cm.P()<<" "<< 1./kine::PDK(_W,0,_target->M())<<" Q2 "<<_photon->M2()<<std::endl;
    return p1cm.P();
  }
    
  inline double TwoBodyProduction::PhaseSpaceFactor_dCosTh() const noexcept {
    // p3/(p1*s)
    return PhaseSpaceNorm()* kinCM_MesonP(_W)/_s/PgammaCM();
  }
   
  inline double TwoBodyProduction::kinCM_MesonP(double W) const {
    //std::cout<<"kinCM_MesonP "<< kine::PDK(W,_meson->P4().M(),_baryon->P4().M()) <<" pdg mass "<<kine::PDK(W,_meson->PdgMass(),_baryon->P4().M())<<" ratio = "<< kine::PDK(W,_meson->P4().M(),_baryon->P4().M())/kine::PDK(W,_meson->PdgMass(),_baryon->P4().M())<<" meson mass diff "<<_meson->P4().M()-_meson->PdgMass()<<" baryon mass diff "<<_baryon->P4().M()<<" mass "<<_meson->PdgMass()<<std::endl;
    return kine::PDK(W,_meson->P4().M(),_baryon->P4().M());

  }
  inline double TwoBodyProduction::kinCM_MesonE(double W) const {
    auto m2_a =_meson->P4().M2();
    auto m2_b =_baryon->P4().M2();
    // std::cout<<"kinCM_MesonE "<< (W*W + m2_a - m2_b)/(2.0*W)<<std::endl;
    return (W*W + m2_a - m2_b)/(2.0*W);
  }
  inline double TwoBodyProduction::kin_tFromWCosTh(double W, double cosTh) const{
    if(Parent()->P4().P()==0){
      std::cerr<<"TwoBodyProduction::kin_tFromWCosTh, parent at rest, we require a valid parent particle to calculate the kinematics"<<std::endl;
      exit(0);
    }
    if(Parent()->P4().M()!=W){
      std::cerr<<"TwoBodyProduction::kin_tFromWCosTh, Parent mass "<<Parent()->P4().M()<<"  != W "<<W<<std::endl;
      exit(0);
	
    }
    auto cmBoost=Parent()->P4().BoostToCM();
    auto p1cm=boost(*_photon,cmBoost);
    // std::cout<<"kin_tFromWCosTh "<<p1cm.M2() + _meson->M2() - 2 * (p1cm.E()* kinCM_MesonE(W)-p1cm.P()* kinCM_MesonP(W)*cosTh)<<std::endl;

    return p1cm.M2() + _meson->M2() - 2 * (p1cm.E()* kinCM_MesonE(W)-p1cm.P()* kinCM_MesonP(W)*cosTh);
  }
  
  inline double TwoBodyProduction::kin_tFromWCosThatQ20(double W, double cosTh) const{
    return kine::tFromcosthW(cosTh,W,0.0,_target->M(),_meson->Mass(),_baryon->Mass());
  }

  inline double TwoBodyProduction::CalcCosThCM() const{
    auto cmBoost=Parent()->P4().BoostToCM();
    auto pm_cm=boost(_meson->P4(),cmBoost); //meson in CM
    auto pg_cm=boost(*_photon,cmBoost); //photon in CM
    // std::cout<<" TwoBodyProduction::CalcCosThCM() "<<Parent()->P4()<<" "<<*_photon<<" "<<_meson->P4()<<std::endl;
    //std::cout<<" TwoBodyProduction::CalcCosThCM() "<<pg_cm<<" "<<pm_cm<<std::endl;
    _cosThCM=TMath::Cos(ROOT::Math::VectorUtil::Angle( pg_cm,pm_cm));
    // in case we need t make sure up-to-date
    _t = kin_tFromWCosTh(_W,_cosThCM);
    //std::cout<<" TwoBodyProduction::CalcCosThCM() "<< _cosThCM<<" "<<_t<<std::endl;
    return  _cosThCM;
  }
  inline double TwoBodyProduction::PhaseSpaceFactor() const{
    return PhaseSpaceFactor_dCosTh();
  }
}
