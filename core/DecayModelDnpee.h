//////////////////////////////////////////////////////////////
///
///Class:		DecayModelDnpee
///Description:
///             Control behaviour of npe+e- final states on deuteron
///             
///
///            Note derived classes should include a constructor to initialise
///            DecayModelDnpee( particle_ptrs , const std::vector<int> pdgs );
#pragma once

#include "DecayModel.h"
#include "PhaseSpaceDecay.h"
#include "FunctionsForElectronScattering.h"
#include "DecayingParticle.h"
#include <TH1D.h>

namespace elSpectro{

 
  class DecayModelDnpee : public DecayModel {

  public:
    
    DecayModelDnpee();
     
    //
    double Intensity() const override; //this should perhaps be final so derived classes cannot overwrite...
    
    void PostInit(ReactionInfo* info) override;

    bool RegenerateOnFail() const noexcept override {return true;};
    bool HasAngularDistribution() override{return false; }

    const Particle* GetNeutron() const noexcept{return _neutron; }
    const Particle* GetProton() const noexcept{return _proton; }
    const Particle* GetElectron() const noexcept{return _ele; }
    const Particle* GetPositron() const noexcept{return _pos; }

    void SetParent(DecayingParticle* pa) override ;
    double PhaseSpaceWeightSq(double W) override{
      return _phaseSpace.PhaseSpaceWeightSq(W);
    }

   double PhaseSpaceFactor() const noexcept {
     // auto fluxPhaseSpace = p1*_W;//eqn 47.28b https://pdg.lbl.gov/2019/reviews/rpp2019-rev-kinematics.pdf
     return PhaseSpaceNorm()/_s/PgammaCMsq();
     
    }
     double PgammaCMsq() const noexcept{
       // std::cout<<"PgammaCMsq() "<<kine::PDK2(_W,0,_target->M())<<" "<<_W<<" "<<_target->M()<<std::endl;
      //  if(_photon->M()==0) return kine::PDK2(_W,0,_target->M());
      auto  pgammaCM= PgammaCM();
      return  pgammaCM* pgammaCM; //for dt phase space factor
    }
    
    double PgammaCM()const noexcept{
      //in case no photon 4-vector yet
      if(_photon->M()==0) return kine::PDK(_W,0,_target->M());
      //else PDK does not qork for virtual photons
      auto cmBoost=Parent()->P4().BoostToCM();
      auto p1cm=boost(*_photon,cmBoost);
      
      //  std::cout<<"PgammaCMsq M "<<_photon->M()<<" PLAB "<<_photon->P()<<" PCM "<<p1cm.P()<<" or "<<kine::PDK(_W,_photon->M(),_target->M())<<" or "<<kine::PDK(_W,_photon->M2(),_target->M())<<" or "<<kine::PDK(_W,0,_target->M())<<std::endl;
      return p1cm.P();
    }
    virtual double MatrixElementsSquared_L() const {return 0; }
    virtual double MatrixElementsSquared_T() const {return 1; } //just real photon by default

    double FindMaxOfIntensity();

  private:
       constexpr double  PhaseSpaceNorm() const {return 1./(2.56819E-6)/64/TMath::Pi();}// Convert from GeV^-2 -> nb
 
    virtual double DifferentialXSect() const{//dont let others call this as need _s, _W and _t set
      //Note if your derived model already gives differential cross section
      //you will need to divide by PhaseSpaceFactor to get MatrixElementSquared from it
    //std::cout<<" DifferentialXSect() "<<_W<<" "<<PhaseSpaceFactor()<<"  MT "<< MatrixElementsSquared_T()<<" "<<PgammaCMsq()<<std::endl;
      return PhaseSpaceFactor() * ( MatrixElementsSquared_T());// + (_photonPol->Epsilon()+_photonPol->Delta())*MatrixElementsSquared_L()); 
    }
       
    Particle* _neutron={nullptr}; 
    Particle* _proton={nullptr};
    Particle* _ele={nullptr};
    Particle*  _pos={nullptr};
    
    LorentzVector* _photon={nullptr};
    LorentzVector* _target={nullptr};//{0,0,0,escat::M_pr()};
    const LorentzVector* _ebeam={nullptr};//{0,0,0,escat::M_pr()};

    PhaseSpaceDecay _phaseSpace;
    ReactionPhotoProd* _prodInfo={nullptr};

    double _Wmax={0};
    double _Wmin={0};
    mutable double _max={0};
    mutable double _dvolume={1};
    mutable double _s={0};
    mutable double _t={0};
    mutable double _W={0};

    bool _isElProd={true};
 
    ClassDefOverride(elSpectro::DecayModelDnpee,1); //class DecayModelst

  };

}
  
