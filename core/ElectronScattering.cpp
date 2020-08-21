#include "ElectronScattering.h"
#include "FunctionsForGenvector.h"
#include "Manager.h"
#include "Interface.h" //for generator
#include "ScatteredElectron_xy.h"
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
    _beamElec{11},
    _beamNucl{ionpdg},
    ProductionProcess{0,decayer,model}
  {
      
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
    _beamElec{11},
    _beamNucl{ionpdg},
    ProductionProcess{0,decayer,model}
  {
  
      SetBeamCondtion();
  }
  /////////////////////////////////////////////////////////////////////
  ElectronScattering::ElectronScattering(double ep,double ionp, DecayModel* model, int ionpdg):
    _pElectron{ep},
    _pIon{ionp},
    _angleElectron{TMath::Pi()},
    _angleIon{0},
    _pdgIon{ionpdg},
    _beamElec{11},
    _beamNucl{ionpdg},
    ProductionProcess{0,nullptr,model}
  {
      
      SetBeamCondtion();
  }
  /////////////////////////////////////////////////////////////////////
  ElectronScattering::ElectronScattering(double ep,double ionp,
		     double eangle,double ionangle,  DecayModel* model, int ionpdg):
    _pElectron{ep},
    _pIon{ionp},
    _angleElectron{eangle},
    _angleIon{ionangle},
    _pdgIon{ionpdg},
    _beamElec{11},
    _beamNucl{ionpdg},
    ProductionProcess{0,nullptr,model}
  {
  
      SetBeamCondtion();
  }

  /////////////////////////////////////////////////////////////////////
  void ElectronScattering::SetBeamCondtion(){
    
    _massIon=TDatabasePDG::Instance()->GetParticle(_pdgIon)->Mass();

    
    _beamElec.SetXYZT(0,0,_pElectron,
		      escat::E_el(_pElectron));
    auto p4=_beamElec.P4();
    genvector::LorentzRotateY(p4,_angleElectron);
    _beamElec.SetP4(p4);
    std::cout<<"ElectronScattering::SetBeamCondtion() Electron "<< _beamElec.P4()<<std::endl;
    
    
    _beamNucl.SetXYZT(0,0,_pIon,
		      TMath::Sqrt(_pIon*_pIon + _massIon*_massIon));
    p4=_beamNucl.P4();
    genvector::LorentzRotateY(p4,_angleIon);
    _beamNucl.SetP4(p4);
    std::cout<<"ElectronScattering::SetBeamCondtion() Nucl "<< _beamNucl.P4()<<std::endl;
 
    //For decaying
    SetXYZT(_beamElec.P4().X(),_beamElec.P4().Y(),
	    _beamElec.P4().Z(),_beamElec.P4().T());


    //For nominal beam conditions set nucleon rest frame vectors
    //in case info is needed in PostInit stage
    //Boost into ion rest frame
    auto prBoost=_beamNucl.P4().BoostToCM();
    _nuclRestNucl=LorentzVector(0,0,0,_beamNucl.Mass());
    _nuclRestElec= boost(_beamElec.P4(),prBoost);

    //set inital lab particles
    AddInitialParticlePtr(&_beamElec);
    AddInitialParticlePtr(&_beamNucl);


   }
  /////////////////////////////////////////////////////////////////////////
  void ElectronScattering::InitGen(){
    //pass on lorentzvectors in nucleon rest frame
    //This is the internal frame for the generator
    _reactionInfo._target=&_nuclRestNucl;
    _reactionInfo._ebeam =&_nuclRestElec;
    
    DecayingParticle::PostInit(dynamic_cast<ReactionInfo*>(&_reactionInfo));
    //now scattered electron should be set
    //_scattered= _reactionInfo._scattered;
    
    auto& unproducts=Model()->UnstableProducts();
    if(unproducts.empty()==true) return;
    
    if(unproducts.size()!=1) {
      std::cerr<<"ElectronScattering::PostInit need a Q2W model with just a gamma*N decay product"<<std::endl;
    }

    _gStarN = unproducts[0];//should only be gamma*N decaying product

    std::cout<<"Electron Scattering min mass "<<_gStarN->MinimumMassPossible()<<std::endl;
    //default scatteredelectron_xy, now have all parameters
    if(Decayer()==nullptr){
      //Need to give ebeam (in ion rest), mass of ion, W threshold
      auto tempDecayer=new ScatteredElectron_xy(_nuclRestElec.P(), _massIon, _gStarN->MinimumMassPossible());
      tempDecayer->SetModel(Model());
      SetDecayer(tempDecayer); //give it to a sink
      mutableDecayer()->PostInit(dynamic_cast<ReactionInfo*>(&_reactionInfo));
      generator().SetModelForMassPhaseSpace(_gStarN->Model());
     }
    
  }
  /////////////////////////////////////////////////////////////////////////
  DecayStatus  ElectronScattering::GenerateProducts(){
    //First, Eventually want to sample from beam divergence distributions
    LorentzVector collision = _beamElec.P4() + _beamNucl.P4();

    //Boost into ion rest frame
    auto prBoost=_beamNucl.P4().BoostToCM();
    collision=boost(collision,prBoost);
    _nuclRestNucl=LorentzVector(0,0,0,_beamNucl.Mass());
    _nuclRestElec= boost(_beamElec.P4(),prBoost);
    
    //set decay parent for e -> e'g*
    SetXYZT(collision.X(),collision.Y(),collision.Z(),collision.T());

    //proceed through decay chain
    //first get masses for all products
    // DetermineAllMasses(); //this will define W threshold...
    //make sure masses are below maximum kinematically possible
    // if(_gStarN->MinimumMassPossible() > W2Max() ){
      // DetermineAllMasses();
    //}
      
    //find a W candidate
    while(DecayingParticle::GenerateProducts()!=DecayStatus::Decayed){
      //std::cout<<"not decayed "<<_nsamples<<std::endl;
      _nsamples++;
    }//DecayModelQ2W
    
    
    //Boost all stable particles back to lab
    Manager::Instance().Particles().BoostStable(-prBoost);

    return DecayStatus::Decayed;
  }

}
