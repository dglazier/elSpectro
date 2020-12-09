//////////////////////////////////////////////////////////////
///
///Class:		ElectronScattering
///Description:
///            Class to manage intial electron scatter
///             1) Initial state particles
///             
///             2) DecayVectors for scattered electron
///                and g*p (e.g. ScatterEelectron_xy)
///             3) DecayModelQ2W for e'(g*p) production
///                to some final state e.g. e'rho+p
///                (e' probably remains unchanged)           

#pragma once
#include "ProductionProcess.h"
#include "PhaseSpaceDecay.h"
#include "FunctionsForElectronScattering.h"

#include <TMath.h> //for Pi()

namespace elSpectro{
  

  class ElectronScattering : public ProductionProcess {

  public:
    
    //or if decayer already give required distribution
   
    ElectronScattering(double ep,double ionp,
		       double anglee,double anglep,
		       DecayVectors* decayer, DecayModel* model=new PhaseSpaceDecay{{},{11,-2211}},int ionpdg=2212);
    ElectronScattering(double ep,double ionp,
		       DecayVectors* decayer, DecayModel* model=new PhaseSpaceDecay{{},{11,-2211}},int ionpdg=2212);

    //Constructors using default ScatteredElectron_xy decayer
   ElectronScattering(double ep,double ionp,
		       double anglee,double anglep,DecayModel* model=new PhaseSpaceDecay{{},{11,-2211}},int ionpdg=2212);
    ElectronScattering(double ep,double ionp,DecayModel* model=new PhaseSpaceDecay{{},{11,-2211}},int ionpdg=2212);

    
    DecayStatus  GenerateProducts( ) override;

 
    void InitGen() override;

    double W2Max()const noexcept{
      return  sqrt(_massIon*_massIon + 2 * (_nuclRestElec.E() -escat::M_el())
		   *(_massIon - escat::M_el()));
    }


    void SetLimit_Q2min(double val){_Q2min=val;}
    void SetLimit_Q2max(double val){_Q2max=val;}
    void SetLimit_Xmin(double val){_Xmin=val;}
    void SetLimit_Xmax(double val){_Xmax=val;}

    //Note these limits apply to target rest frame!!
    void SetLimitTarRest_ePmin(double val){_ePmin=val;}
    void SetLimitTarRest_ePmax(double val){_ePmax=val;}
    void SetLimitTarRest_eThmin(double val){_eThmin=val;}
    void SetLimitTarRest_eThmax(double val){_eThmax=val;}


    double dsigma() const override {return _gStarN->Model()->dsigma()* Decayer()->dsigma();}
    //double dsigma() const override {return _gStarN->Model()->dsigma();}
    // double dsigma() const override {return Decayer()->dsigma();}
    
  private:
    
    ElectronScattering()=default;

    
    void SetBeamCondtion();
    
    double _pElectron={0}; //nominal e- beam energy
    double _pIon={0}; //nominal ion beam energy
    double _angleElectron={0}; //nominal electron crossing angle
    double _angleIon={0};  //nominal proton crossing angle
    double _massIon={0};  //nominal proton crossing angle
    int _pdgIon={2212}; //species of ion

    //user specified limits
    double _Q2min={0};
    double _Q2max={0};
    double _eThmin={0};
    double _eThmax={0};
    double _ePmin={0};
    double _ePmax={0};
    double _Xmin={0};
    double _Xmax={0};

    
    long _nsamples=0;

    Particle _beamElec;
    Particle _beamNucl;
    LorentzVector _nuclRestElec;
    LorentzVector _nuclRestNucl;

    DecayingParticle* _gStarN={nullptr}; 
    //PhotoProdInfo _info;
    ReactionElectroProd _reactionInfo;
    
    
    ClassDefOverride(elSpectro::ElectronScattering,1); //class ElectronScattering
 
  };


}//namespace elSpectro

