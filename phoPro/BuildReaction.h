#pragma once
#include "Interface.h"
#include "JpacTwoBody.h"
#include "CppHelperFuncs.h"
#include "ElectronScattering.h"
#include "DecayChannel.h"
#include "FormationQ2W.h"
#include "CollidingParticle.h"

namespace phoPro{

  elSpectro::ElectronScattering*  Build_ep_Collision(double ebeamE,double pbeamE, elSpectro::JpacTwoBody& jpac){
    auto elBeam =  elSpectro::CollidingParticle{11,ebeamE};
    auto elBeamP4=elBeam.GetInteracting4Vector();
    elBeam.SetAngleThetaPhi(TMath::Pi(),0);

    //define pr beam, pdg =2212
    auto prBeam =  elSpectro::CollidingParticle{2212,pbeamE};
    auto prBeamP4=prBeam.GetInteracting4Vector();

    //Decay g*p state
    //new FormationQ2W{0,&jpac};
    
    //combine beam, target and reaction products
    return new elSpectro::ElectronScattering(elBeam,prBeam, elSpectro::CloneModel(elSpectro::FormationQ2W{0,cpp::MakeBaseShared< elSpectro::DecayModel>(jpac)} ) );
 
  }
}
