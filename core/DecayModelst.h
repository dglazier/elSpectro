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

    const Particle* GetMeson() const noexcept{return _meson; }
    const Particle* GetBaryon() const noexcept{return _baryon; }

    void SetUseSDME(bool use=true){_useSDME=use;}
    

    double PgammaCMsq()const noexcept{
      //in case no photon 4-vector yet
      if(_photon->M()==0) return kine::PDK2(_W,0,_target->M());
      //else PDK does not qork for virtual photons
      auto cmBoost=Parent()->P4().BoostToCM();
      auto p1cm=boost(*_photon,cmBoost);
      return p1cm.P()*p1cm.P(); //for dt phase space factor
    }

    const ReactionElectroProd* ProductionInfo() const { return _prodInfo; }

    void HistIntegratedXSection(TH1D& hist);
  protected:
    
    virtual double MatrixElementsSquared_L() const {return 0; }
    virtual double MatrixElementsSquared_T() const {return 1; } //just real photon by default
    
    constexpr double  PhaseSpaceNorm() const {return 1./(2.56819E-6)/64/TMath::Pi();}
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
    
       
    SDME* const  GetMesonSDMEs() const {return _sdmeMeson;}
    SDME* const  GetBaryonSDMEs() const {return _sdmeBaryon;}

    virtual void CalcMesonSDMEs() const {};
    virtual void CalcBaryonSDMEs() const {};

    double FindMaxOfIntensity();

    double get_s() const noexcept{ return _s; }
    double get_t() const noexcept { return _t; }
    double get_W() const noexcept { return _W; }
 
  private:

    double DifferentialXSect() const{//dont let others call this as need _s, _W and _t set
      //Note if your derived model already gives differential cross section
      //you will need to divide by PhaseSpaceFactor to get MatrixElementSquared from it
      return PhaseSpaceFactor() * ( MatrixElementsSquared_T() +
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
 
    bool _useSDME={false};

    ClassDefOverride(elSpectro::DecayModelst,1); //class DecayModelst
    
  };//class DecayModelst

}//namespace elSpectro
