//////////////////////////////////////////////////////////////
///
///Class:		DecayModelstM2
///Description:
///             Control behaviour of Particle decay to Particle products
///             Defined by
///             1) Two body phase space as a fucntion of s and t
///             
///
///            Note derived classes should include a constructor to initialise
///            DecayModelstM2( particle_ptrs , const std::vector<int> pdgs );
#pragma once

#include "DecayModelst.h"
#include "SDME.h"
#include "FunctionsForElectronScattering.h"
#include "DecayingParticle.h"
#include <TH1D.h>

namespace elSpectro{

 
  class DecayModelstM2 : public DecayModelst {

  public:
    
    DecayModelstM2()=delete;
    //constructor giving jpac amplitude pointer (which we will now own)
    //and decay particles 
    DecayModelstM2(  particle_ptrs parts,const std::vector<int> pdgs  );
    
    //
    double Intensity() const override; //this should perhaps be final so derived classes cannot overwrite...
    
    void PostInit(ReactionInfo* info) override;

    bool RegenerateOnFail() const noexcept override {return true;};
    bool HasAngularDistribution() override{return false; }

    const Particle* GetMeson() const noexcept{return _meson; }
    const Particle* GetBaryon() const noexcept{return _baryon; }

    

    /* double PgammaCMsq()const noexcept{*/
    // double PgammaCMsq() const noexcept{
    //   if(_photon->M()==0) return kine::PDK2(_W,0,_target->M());
    //   auto  pgammaCM= PgammaCM();
    //   return  pgammaCM* pgammaCM; //for dt phase space factor
    // }
    
    // double PgammaCM()const noexcept{
    //   //in case no photon 4-vector yet
    //   if(_photon->M()==0) return kine::PDK(_W,0,_target->M());
    //   //else PDK does not qork for virtual photons
    //   auto cmBoost=Parent()->P4().BoostToCM();
    //   auto p1cm=boost(*_photon,cmBoost);
      
    //   //     std::cout<<"PgammaCMsq M"<<_photon->M()<<" PLAB "<<_photon->P()<<" PCM "<<p1cm.P()<<" or "<<kine::PDK(_W,_photon->M(),_target->M())<<" or "<<kine::PDK(_W,_photon->M2(),_target->M())<<" or "<<kine::PDK(_W,0,_target->M())<<std::endl;
    //   return p1cm.P();
    // }
    
   
    const ReactionPhotoProd* ProductionInfo() const { return _prodInfo; }
    //    const ReactionElectroProd* ProductionInfo() const { return _prodInfo; }
    
    void HistIntegratedXSection_ds(TH1D& hist);
    void HistIntegratedXSection(TH1D& hist);
    void HistMaxXSection(TH1D& hist);
    
    // double PhaseSpaceFactor() const noexcept {
    //   /* auto fluxPhaseSpace = p1*_W;//eqn 47.28b https://pdg.lbl.gov/2019/reviews/rpp2019-rev-kinematics.pdf
    // 	 auto ans =  1./fluxPhaseSpace
    // 	 * kine::PhaseSpaceFactorDt(_W,p1,_meson->Mass(),_baryon->Mass())
    // 	 * kine::PDK(_W,_meson->Mass(),_baryon->Mass())/_W
    // 	 * PhaseSpaceNorm();//nbarn it
    // 	 */ //Note above full calculation simplifies to
    //   //return PhaseSpaceNorm()/_s/kine::PDK2(_W,_photon->M(),_target->M());
    //   //Please note kine::PDK2(_W,_photon->M(),_target->M()) does not give
    //   //correct momentum
    //   return PhaseSpaceNorm()/_s/PgammaCMsq();
    //   //this would not be the case if the differential was dcosth rather than t
    // }
    
    double PhaseSpaceCorrectForQ2() const noexcept {
      auto cmBoost=(_photon+_target).P4().BoostToCM();
      auto p1cm=boost(*_photon,cmBoost);
      return p1cm.P()/kine::PDK(get_W(),0,_target.M());
    }
 
  protected:
    
  
    double FindMaxOfIntensity();

  public:
    // double get_s() const noexcept{ return _s; }
    // double get_t() const noexcept { return _t; }
    // double get_W() const noexcept { return _W; }
    // double get_Q2() const noexcept { return -_photon->M2(); }
    double get_M2() const noexcept { return _baryon->M2(); }

  
    //    void HistIntegratedXSection(TH1D& hist) override=0;
    
    // based on jpac_photon :: inclusive_kinematics
    inline double M2fromXCOS(double x, double cos)
    {
      auto r = (x/(costh));
      return M2fromRCOS(r,costh);
    }

    // Energy of X in (x, y)
    inline double EfromXY(double x, double y)
    {
      //max momentum when minimum mass of baryon
      auto pMax =kine::pFromKallen(_s,MX2,_minMB*_minMB); 
      double e2 = _mX2 + pMax*pMax * (x*x + y*y);
      return sqrt(e2);
    };
    inline double M2fromXY(double x, double y)
    {
      auto MX2= _meson->M2();
      return _s + MX2 - 2. * _W * EfromXY(x, y);
    }; 
    // Similarly, momentum transfer from XY
    inline double TfromXY(double x, double y)
    {
      auto qGamma= _photon->P();
      auto pMax =kine::pFromKallen(_s,MX2,_minMB*_minMB); 
      return _mX2 - 2.*qGamma * (EfromXY(x, y) - pMax * x);
    };
       
    inline double M2fromRCOS(double r, double cos)
    {

      auto MT2 =  _target->M2();
      auto MX2= _meson->M2();
      //max momentum when minimum mass of baryon
      auto pMax =kine::pFromKallen(_s,MX2,_minMB*_minMB); 
      return _M2 = MX2 + _s - 2. * TMath::Sqrt(MX2 * _s + pMax * pMax * r * r * _s);
    };
    
    inline double TfromRCOS(double r, double costh)
    {
      double M2 = M2fromRCOS(r, cos);
      auto MT2 =  _target->M2();
      auto MX2= _meson->M2();
      auto pMax =kine::pFromKallen(_s,MX2,_minMB*_minMB); 
       double t = _mX2 - (_s - MT2) * (_s - M2 + MX2) / (2. * _s) + 2.* _photon->P() * pMax() * r * costh;
      return t;
    };
   double M2fromTX(double t, double x)
    {
      auto MT2 =  _target->M2();
      auto MX2= _meson->M2();
      double lami = kine::Kallen(_s, 0.,MT2);
      double lamf = kine::Kallen(_s, MX2, _minMB*_minMB);
      double num = MT2 * MX2 + MT2 * _s + MX2 * _s - _s*_s - 2.*_s*t + sqrt(lami * lamf) * x;
      return num / (MT2 - _s);
    };

  private:

    double DifferentialXSect() override const{//dont let others call this as need _s, _W and _t set
      //ignore possible longitudinal component
      return d3sigma();//differential stdM2
    }
       
    PhotonPolarisationVector* _photonPol={nullptr};
 

    ReactionPhotoProd* _prodInfo={nullptr};
 
    Particle* _baryon={nullptr};
    Particle* _meson={nullptr};
    LorentzVector* _photon={nullptr};
    LorentzVector* _target={nullptr};//{0,0,0,escat::M_pr()};
    const LorentzVector* _ebeam={nullptr};//{0,0,0,escat::M_pr()};
 
    mutable double _max={0};
    mutable double _s={0};
    mutable double _t={0};
    mutable double _W={0};
    mutable double _dt={0};
    mutable double _dsigma={0};
    double _Wmax={0};
    double _minMassB=0;
    
    bool _useSDME={false};
    bool _isElProd={true};
    
    ClassDefOverride(elSpectro::DecayModelstM2,1); //class DecayModelstM2
    
  };//class DecayModelstM2

}//namespace elSpectro
