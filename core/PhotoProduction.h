//////////////////////////////////////////////////////////////
///
///Class:		PhotoProduction
///Description:
///            Class to manage intial real photoproduction
///             1) Initial state particles
///             
///             2) DecayVectors for bremsstrahlung electron
///                and gp
///             3) DecayModel0W for (gp) production
///                to some final state e.g. rho+p
///   

#pragma once
#include "ProductionProcess.h"
#include "PhaseSpaceDecay.h"

#include <TMath.h> //for Pi()

namespace elSpectro{
  
  
  class PhotoProduction : public ProductionProcess {

  public:
    static int NintegralsPhotoProduction;// keep track of number of integrations in same root sesssion

 
    PhotoProduction(CollidingParticle *electron,CollidingParticle* target,  DecayModel* model);
    
    DecayStatus  GenerateProducts( ) override;

 
    void InitGen() override;

  

    void SetLimit_Emin(double val){_Emin=val;}
    void SetLimit_Emax(double val){_Emax=val;}

   

    double IntegrateCrossSection() override;
    double IntegrateCrossSectionFast() override;
    LorentzVector MakeCollision();

    void SetCacheIntegrals(int doit=1){_cacheIntegrals=doit;}
    
  private:
    
    PhotoProduction()=delete;

    
    void SetBeamCondtion();
    void SetNominalBeamCondtion();

    Particle _beamPhot;
    Particle _beamNucl;
    LorentzVector _nuclRestPhot;
    LorentzVector _nuclRestNucl;
    ReactionPhotoProd _reactionInfo;

    double _pElectron={0}; //nominal e- beam energy
    double _pIon={0}; //nominal ion beam energy
    double _angleElectron={0}; //nominal electron crossing angle
    double _angleIon={0};  //nominal proton crossing angle
    double _massIon={0};  //nominal proton crossing angle
    double _Wmin={0};  //minumum CM mass
    
    double _Emin={0};
    double _Emax={0};

    
    long _nsamples=0;
    int _pdgIon={2212}; //species of ion

    short _cacheIntegrals={0};
    
    DecayingParticle* _gammaN={nullptr}; 
    CollidingParticle* _photonptr={nullptr};
    CollidingParticle* _targetptr={nullptr};
    
    ClassDefOverride(elSpectro::PhotoProduction,1); //class PhotoProduction
 
  };


}//namespace elSpectro

