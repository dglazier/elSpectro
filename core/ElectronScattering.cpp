#include "ElectronScattering.h"
#include "FunctionsForElectronScattering.h"
#include "FunctionsForGenvector.h"
#include "Manager.h"
#include <TDatabasePDG.h>


namespace elSpectro{

  // ElectronScattering::ElectronScattering(double ep,double ionp,
  // 					 double eangle,double ionangle,
  // 					 DecayModel* model, int ionpdg):
  //   _pElectron{ep},
  //   _pIon{ionp},
  //   _angleElectron{eangle},
  //   _angleIon{ionangle},
  //   _pdgIon{ionpdg},
  //   ProductionProcess{model}
  // {
  //   SetBeamCondtion();
  // }
  // //////////////////////////////////////////////////////////////////
  // ElectronScattering::ElectronScattering(double ep,double ionp,
  // 					 DecayModel* model, int ionpdg):
  //   _pElectron{ep},
  //   _pIon{ionp},
  //   _angleElectron{TMath::Pi()},
  //   _angleIon{0},
  //   _pdgIon{ionpdg},
  //   ProductionProcess{model}
  // {
  //   SetBeamCondtion();
  // }
  /////////////////////////////////////////////////////////////////////
  ElectronScattering::ElectronScattering(double ep,double ionp, DecayVectors* decayer,  DecayModel* model, int ionpdg):
    _pElectron{ep},
    _pIon{ionp},
    _angleElectron{TMath::Pi()},
    _angleIon{0},
    _pdgIon{ionpdg},
    ProductionProcess{0,decayer,model}{
      
      SetBeamCondtion();
    }
  /////////////////////////////////////////////////////////////////////
  ElectronScattering::ElectronScattering(double ep,double ionp,
		     double eangle,double ionangle,
		     DecayVectors* decayer,  DecayModel* model, int ionpdg):
    _pElectron{ep},
    _pIon{ionp},
    _angleElectron{eangle},
    _angleIon{ionangle},
    _pdgIon{ionpdg},
    ProductionProcess{0,decayer,model}
  {
  
      SetBeamCondtion();
  }

  /////////////////////////////////////////////////////////////////////
  void ElectronScattering::SetBeamCondtion(){
    
    _massIon=TDatabasePDG::Instance()->GetParticle(_pdgIon)->Mass();

    
    _beamElec.SetXYZT(0,0,_pElectron,
		      escat::E_el(_pElectron));
     genvector::LorentzRotateY(_beamElec,_angleElectron);
    std::cout<<"ElectronScattering::SetBeamCondtion() Electron "<< _beamElec<<std::endl;
    
    
    _beamNucl.SetXYZT(0,0,_pIon,
		      TMath::Sqrt(_pIon*_pIon + _massIon*_massIon));
    genvector::LorentzRotateY(_beamNucl,_angleIon);
    std::cout<<"ElectronScattering::SetBeamCondtion() Nucl "<< _beamNucl<<std::endl;
 
    //For decaying
    SetXYZT(_beamElec.X(),_beamElec.Y(),
	    _beamElec.Z(),_beamElec.T());

  }

  /////////////////////////////////////////////////////////////////////////
  double ElectronScattering::GenerateProducts(const CurrentEventInfo* parentInfo){
    //First, Eventually want to sample from beam divergence distributions
    LorentzVector collision = _beamElec + _beamNucl;

    //Boost into ion rest frame
    auto prBoost=_beamNucl.BoostToCM();
    collision=boost(collision,prBoost);
  
    //set decay parent for e -> e'g*
    SetXYZT(collision.X(),collision.Y(),collision.Z(),collision.T());

    //proceed through decay chain
    DecayingParticle::GenerateProducts(parentInfo);
    
    //Boost all stable particles back to lab
    Manager::Instance().Particles().BoostStable(-prBoost);

    return 1.0;
  }

}
