//////////////////////////////////////////////////////////////
///
///Class:		FormationQ2W
///Description:
///             Control behaviour of Particle decay to Particle products
///             Defined by
///             1) list of Particle products
///             2) Intensity function dependent on Q2 and W
///
///            Note derived classes should include a constructor to initialise
///            FormationQ2W( particle_ptrs , const std::vector<int> pdgs );
#pragma once

#include "Formation.h"

namespace elSpectro{

  
  class FormationQ2W : public Formation {

  public:
    
    FormationQ2W()=delete;
    //delete default constructor so have to use threshold
    //so other 5 constructors also defaulted(rule of 5)

    //constructor giving W theshold, just produces scatted electron kinematics
    // FormationQ2W(  double thresh  );
    //constructor giving W theshold and subsequent primary decay of Nucl+gamma* system
    //FormationQ2W(  double thresh, DecayModel* gNmodel,DecayVectors* gNdecayer=nullptr);
    FormationQ2W(  double thresh, decaymodel_ptr gNmodel,decayer_ptr gNdecayer=CloneDecayer(TwoBodyFlat()));
    
    // Each model must define its intensity
    double Intensity() const override;

    void PostInit(ReactionInfo* info) override;

    bool RegenerateOnFail() const  noexcept override {return false;}


    const Particle& GetScatteredElectron() const noexcept{return *static_cast<const DecayingParticle*>(Product(_idxElectron)); }
  

    double getQ2() const noexcept{return -_gamma.M2();}

     //Q2 dependence from The H1 Collaboration: Elastic electroproduction of ρ mesons at HERA eqn 49 https://link.springer.com/content/pdf/10.1007/s100520000150.pdf
    constexpr double Q2H1RhoAt0() const  noexcept {return 3.0610097;}//1./TMath::Power((0.77549000*0.77549000),2.2); 
    double Q2H1Rho() const noexcept {return 1./TMath::Power((getQ2() + 0.77549000*0.77549000),2.2)/Q2H1RhoAt0(); }

    // void CreateExcitationSpectra();

  protected:
    
    //mutable PhotoProdInfo _myInfo;//!
    ReactionElectroProd* ElectroProdInfo() const noexcept {return  static_cast<ReactionElectroProd*> (_reactInfo);}
    
    uint _idxElectron=0;

  private:
    
    void Init();
    
   
    ClassDefOverride(elSpectro::FormationQ2W,1); //class FormationQ2W
    
  };//class FormationQ2W

}//namespace elSpectro
