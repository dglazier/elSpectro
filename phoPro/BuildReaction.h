#pragma once
#include "Interface.h"
#include "JpacTwoBody.h"
#include "CppHelperFuncs.h"
#include "ElectronScattering.h"
#include "DecayChannel.h"
#include "FormationQ2W.h"
#include "CollidingParticle.h"
#include "NuclearBreakup.h"
#include "QuasiFreeNucleon.h"

namespace phoPro{

  elSpectro::ElectronScattering*  Build_ep_Collision(double ebeamE,double pbeamE, elSpectro::JpacTwoBody& jpac){
    auto Nmass = (0.93827210);
    double pbeamP = sqrt(pbeamE*pbeamE - Nmass*Nmass);
   
    double ebeamP = sqrt(ebeamE*ebeamE - 0.00051099900*0.00051099900);

    auto elBeam =  elSpectro::CollidingParticle{11,ebeamP};
    elBeam.SetAngleThetaPhi(TMath::Pi(),0);

    //define pr beam, pdg =2212
    auto prBeam =  elSpectro::CollidingParticle{2212,pbeamP};

    //Decay g*p state
    //new FormationQ2W{0,&jpac};
    
    //combine beam, target and reaction products
    return new elSpectro::ElectronScattering(elBeam,prBeam, elSpectro::CloneModel(elSpectro::FormationQ2W{0,cpp::MakeBaseShared< elSpectro::DecayModel>(jpac)} ) );
 
  }
  
  elSpectro::ElectronScattering*  Build_eD_Collision(double ebeamE,double NbeamE, int tarPDG, elSpectro::JpacTwoBody& jpac){
    //convert energy to momenta
    auto Nmass = (0.93827210+0.93956540)/2;
    double NbeamP = sqrt(NbeamE*NbeamE - Nmass*Nmass);
    double DbeamP = 2*NbeamP;
    double ebeamP = sqrt(ebeamE*ebeamE - 0.00051099900*0.00051099900);
    
    auto elBeam =  elSpectro::CollidingParticle{11,ebeamP};
    elBeam.SetAngleThetaPhi(TMath::Pi(),0);

    // construct coliding particle with species=nucleonPDG from a nuclear breakup
    // to proton and neutron
    // mometum = 0 
    // Nucleus = 1000010020 = deuteron
    // use default initial nucleon momentum distribution
    //auto qfTarget = initial(nucleonPDG,0,1000010020,model(new NuclearBreakup(2212,2112)),new QuasiFreeNucleon()); //approximate paramterisation
    //Use CDBonn potential to give nucleon momentum
    auto qfTarget = elSpectro::CollidingParticle{tarPDG,DbeamP,1000010020,
						 elSpectro::CloneModel(elSpectro::NuclearBreakup{2212,2112}),
						 elSpectro::CloneDecayer(elSpectro::QuasiFreeNucleon(elSpectro::QuasiFree::CDBonnMomentum()))
    };//CDBonn potential

    //get 4-momentum of target nucleon
    //  auto nucleon=qfTarget->GetInteracting4Vector();
 
   //Decay g*p state
    //new FormationQ2W{0,&jpac};
    
    //combine beam, target and reaction products
    
    return new elSpectro::ElectronScattering(elBeam,qfTarget, elSpectro::CloneModel(elSpectro::FormationQ2W{0,cpp::MakeBaseShared< elSpectro::DecayModel>(jpac)} ) );
 
  }

}
