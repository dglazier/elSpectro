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


    //For nominal beam conditions set nucleon rest frame vectors
    //in case info is needed in PostInit stage
    //Boost into ion rest frame
    auto prBoost=_beamNucl.BoostToCM();
    _nuclRestNucl=LorentzVector(0,0,0,_beamNucl.M());
    _nuclRestElec= boost(_beamElec,prBoost);
    
  }
  /////////////////////////////////////////////////////////////////////////
  void ElectronScattering::InitGen(){
    //pass on lorentzvectors in nucleon rest frame
    //This is the internal frame for the generator
    _reactionInfo._target=&_nuclRestNucl;
    _reactionInfo._ebeam =&_nuclRestElec;
    DecayingParticle::PostInit(dynamic_cast<ReactionInfo*>(&_reactionInfo));
    auto& unproducts=Model()->UnstableProducts();
    if(unproducts.size()!=1) {
      std::cerr<<"ElectronScattering::PostInit need a Q2W model with just a gamma*N decay product"<<std::endl;
    }
    _gStarN = unproducts[0];//should only be gamma*N decaying product 
  }
  /////////////////////////////////////////////////////////////////////////
  DecayStatus  ElectronScattering::GenerateProducts(){
    //First, Eventually want to sample from beam divergence distributions
    LorentzVector collision = _beamElec + _beamNucl;

    //Boost into ion rest frame
    auto prBoost=_beamNucl.BoostToCM();
    collision=boost(collision,prBoost);
    _nuclRestNucl=LorentzVector(0,0,0,_beamNucl.M());
    _nuclRestElec= boost(_beamElec,prBoost);
    
    //set decay parent for e -> e'g*
    SetXYZT(collision.X(),collision.Y(),collision.Z(),collision.T());

    //proceed through decay chain
    
    //in case I have any useful info to pass to decay products
    //e.g. polarisations, SDMEs, moments,...
      //   const CurrentEventInfo* myInfo ={nullptr}; 

     
    long nProdsamples = 0 ;
    long nDecaysamples = 0 ;

    //generate scattered electron
    // DecayStatus status = DecayStatus::ReGenerate;
    //while(status == DecayStatus::ReGenerate){

    //find a W candidate
    while(DecayingParticle::GenerateProducts()!=DecayStatus::Decayed){
      nProdsamples++;
    }//DecayModelQ2W
    
    //  std::cout<< "ElectronScattering  " <<nProdsamples<<std::endl;
  
      //Decay gStarN->Products, may need to regenerate if this model intensity depends on W (in this case any DecayModelQ2W sWeight (if not =1)  must be divided out from this model's weight)
      //while((status=_gStarN->GenerateProducts()) == DecayStatus::TryAnother) {
    //	nDecaysamples++;
    // } 
    // if(status==DecayStatus::Decayed) break;
    // }
    
    // std::cout<<"Needed "<<nProdsamples<<" eletecrons for this event and "<<nDecaysamples<<" decays"<<std::endl;
    
    //Boost all stable particles back to lab
    Manager::Instance().Particles().BoostStable(-prBoost);

    return DecayStatus::Decayed;
  }

}
