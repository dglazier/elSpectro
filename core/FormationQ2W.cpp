#include <TDatabasePDG.h>
#include "FormationQ2W.h"
#include "FunctionsForElectronScattering.h"


namespace elSpectro{
  //////////////////////////////////////////////////////
  ////Constructor for e- scattering kinematics only
  // FormationQ2W::FormationQ2W( double thresh) :
  //   _threshold{thresh},
  //   Formation{{},{-2211,11}}
  // {
  //   _name={"FormationQ2W"};

  //   Init();
  // }
  ///////////////////////////////////////////////////////
  //Create a decay model decaying to g*N(-2211) and e'
  FormationQ2W::FormationQ2W( double thresh) :
    Formation{thresh,{ DecayingParticle{-2211} },{11}}
  {
    _name={"FormationQ2W"};
    Init();
  }
  ///////////////////////////////////////////////////////
  //Create a decay model decaying to g*N(-2211) and e'
  FormationQ2W::FormationQ2W( double thresh,
			      decaymodel_ptr gNmodel,decayer_ptr gNdecayer) :
    Formation{thresh,{ DecayingParticle{-2211,gNmodel,gNdecayer} },{11}}
  {
    _name={"FormationQ2W"};
    Init();
  }
  ////////////////////////////////////////////////////////
  ///complete constructor
  void FormationQ2W::Init(){
    std::cout<<" FormationQ2W::Init() first "<<Products()[0]->Pdg()<<std::endl;
    if(Products()[0]->Pdg()==11){
      // _gstarNuc=Products()[1]; //-2211
      // _electron=Products()[0]; //11
      _idxGStarNuc=1;
      _idxElectron=0;
    }
    else{
      //_gstarNuc=Products()[0]; //-2211
      // _electron=Products()[1]; //11
      _idxGStarNuc=0;
      _idxElectron=1;
  
    }
 
    // GetScatteredElectron().Print();
    // GetGammaN().Print();
    
    /* auto gNprods=GetGammaN().Model()->Products();
    
    if( TString("Baryon")==TDatabasePDG::Instance()
	->GetParticle(gNprods[0]->Pdg())->ParticleClass() ){
      //Make sure meson is product 0 and baryon product 1
      dynamic_cast<DecayingParticle*>(_gstarNuc)->Model()->SwapProducts(0,1);
    }
    */
    
  }

  ////////////////////////////////////////////////////////
  void FormationQ2W::PostInit(ReactionInfo* info){
    
    _reactInfo = info;
    auto    prodInfo = ElectroProdInfo() ;
    //std::cout<<"FormationQ2W::PostInit "<<_electron->P4ptr()<<" "<<_gstarNuc->P4ptr()<<std::endl;
    if( prodInfo==nullptr){
      std::cerr<<"FormationQ2W PostInit not an ElectronScattering reaction"<<std::endl;
      exit(0);
      }

    /*
    if(prodInfo->Wmax()==0){
      //For backward compatability, should be done with
      //colliding particles now and set in ProductionProcess::PostInit().
      prodInfo->_Wmax=( (prodInfo->_target) + (prodInfo->_ebeam) ).M();
      std::cout<<"FormationQ2W wax "<< prodInfo ->_target<<" "<<prodInfo->_ebeam<<" "<<prodInfo->_Wmax<<std::endl;exit(0);
      }
    */
    // auto gNprods=dynamic_cast<DecayingParticle*>(_gstarNuc)->Model()->Products();

      // _prodInfo->_baryon=gNprods[1]->P4ptr();
      // _prodInfo->_meson=gNprods[0]->P4ptr();
      
      
      DecayModel::PostInit(prodInfo);
      
      //upate in case minumim mass changed..
      if(_threshold<MinimumMassPossible() )_threshold=MinimumMassPossible(); 
      
      CreateExcitationSpectra();
      
  }
  
  ////////////////////////////////////////////////////////
  double  FormationQ2W::Intensity() const{
    
    double W = getW();
    if(TMath::IsNaN(W)) return 0.0;
    if(W  < _threshold ) return 0.;
 
    auto prodInfo = ElectroProdInfo() ;
     
    //calculate virtual photon
    auto p4beam=(prodInfo->_ebeam);
    auto p4tar=(prodInfo->_target);
    auto p4scat=GetScatteredElectron().P4();
    
    _gamma = p4beam-p4scat;//can now use getQ2
    // std::cout<<"FormationQ2W::Intensity() "<<p4beam+p4tar<<" beam "<<p4beam<<" tar "<<p4tar<<" scat "<<p4scat<<_gamma<<" sum ebeam "<<p4scat+_gamma<<std::endl;//calculate photon polarisation
    auto epsilon = escat::virtualPhotonPolarisation(p4beam,p4tar,p4scat);
    //protect divide by 0
    auto delta = (epsilon==1)? 0: 2*escat::M2_el()/getQ2()*(1-epsilon);
    
    _photonPol.SetEpsilon(epsilon);
    _photonPol.SetDelta(delta);
    
    //   Get envelope weight from integrated cross section
    double weight=1.0;
    
    weight = CurrCrossSection().GetWeightFor( W  );
    
    
    weight*=Q2H1Rho();
    //std::cout<<"FormationQ2W "<<weight<<" "<<getQ2()<< " "<<W<<std::endl;
    
    //copy all currently known particle info
    prodInfo->_scattered=p4scat;
    //_prodInfo->_photoN=_gstarNuc->P4();
    prodInfo->_photon=_gamma;
    prodInfo->_photonPol=_photonPol;
    
     
    return weight;
  }
    
  
}
