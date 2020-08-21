#include "DecayModelQ2W.h"
#include "DecayingParticle.h"
#include "FunctionsForElectronScattering.h"
#include <TDatabasePDG.h>


namespace elSpectro{
  //////////////////////////////////////////////////////
  ////Constructor for e- scattering kinematics only
  DecayModelQ2W::DecayModelQ2W( double thresh) :
    _threshold{thresh},
    DecayModel{{},{-2211,11}}
  {
    _name={"DecayModelQ2W"};

    Init();
  }
  ///////////////////////////////////////////////////////
  ///constructor includes subseqent decay of Ngamma* system
  /* DecayModelQ2W::DecayModelQ2W( double thresh, DecayModel* gNmodel) :
    _threshold{thresh},
    DecayModel{{ new DecayingParticle{-2211,primary} },{11}}
  {
    _name={"DecayModelQ2W_with_primary_decay"};
    Init();
    }*/
  ///////////////////////////////////////////////////////
  ///constructor includes subseqent decay of Ngamma* system
  DecayModelQ2W::DecayModelQ2W( double thresh,
				DecayModel* gNmodel,DecayVectors* gNdecayer) :
    _threshold{thresh},
    DecayModel{{ new DecayingParticle{-2211,gNmodel,gNdecayer} },{11}}
  {
    _name={"DecayModelQ2W_with_primary_decay_and_decayer"};
    Init();
  }
  ////////////////////////////////////////////////////////
  ///complete constructor
  void DecayModelQ2W::Init(){
    if(Products()[0]->Pdg()==11){
      _gstarNuc=Products()[1]; //-2211
      _electron=Products()[0]; //11 
    }
    else{
      _gstarNuc=Products()[0]; //-2211
      _electron=Products()[1]; //11
    }
 
    _electron->Print();
    _gstarNuc->Print();
    
    auto gNprods=dynamic_cast<DecayingParticle*>(_gstarNuc)->Model()->Products();
    if( TString("Baryon")==TDatabasePDG::Instance()
	->GetParticle(gNprods[0]->Pdg())->ParticleClass() ){
      //Make sure meson is product 0 and baryon product 1
      dynamic_cast<DecayingParticle*>(_gstarNuc)->Model()->SwapProducts(0,1);
    }
 
    if(_threshold<MinimumMassPossible() )_threshold=MinimumMassPossible(); 
  }

  ////////////////////////////////////////////////////////
  void DecayModelQ2W::PostInit(ReactionInfo* info){
      _prodInfo = dynamic_cast<ReactionElectroProd*> (info);
      //std::cout<<"DecayModelQ2W::PostInit "<<_electron->P4ptr()<<" "<<_gstarNuc->P4ptr()<<std::endl;
      if( _prodInfo==nullptr) std::cerr<<"DecayModelQ2W PostInit not an ElectronScattering reaction"<<std::endl;
      _prodInfo->_scattered=_electron->P4ptr();
      _prodInfo->_photoN=_gstarNuc->P4ptr();
      _prodInfo->_photon=&_gamma;
      _prodInfo->_photonPol=&_photonPol;
      
      auto gNprods=dynamic_cast<DecayingParticle*>(_gstarNuc)->Model()->Products();
      // std::cout<<"DecayModelQ2W::PostInit "<<_prodInfo<<" "<<gNprods[0]->Pdg()<<" "<<gNprods[1]->Pdg()<<std::endl;
      
      _prodInfo->_baryon=gNprods[1]->P4ptr();
      _prodInfo->_meson=gNprods[0]->P4ptr();
      
      
      //std::cout<<"DecayModelQ2W::PostInit "<<std::endl;
      DecayModel::PostInit(_prodInfo);
     // std::cout<<"DecayModelQ2W::PostInit done"<<std::endl;
    
    }
  
  ////////////////////////////////////////////////////////
  double  DecayModelQ2W::Intensity() const{
    // std::cout<<"DecayModelQ2W::Intensity "<<MinimumMassPossible()<<" "<<ParentVector().M()<<std::endl;
    if(CheckThreshold()==false){
      return 0.;
    }
    if(GetGammaN()->P4().M()<GetGammaN()->MinimumMassPossible()) return 0;
    if(_gstarNuc->P4().M() < _threshold ) return 0.;
 
    //auto parent = ParentVector(); //e' + g*N
    //  auto parent = _electron->P4() + _gstarNuc->P4();

    //  std::cout<<"DecayModelQ2W::Intensity "<<parent<<" "<<_electron->P4()<<std::endl;
    //note we are in the nuc rest frame, so parent momentum = e- beam momentum
    // auto ebeam = parent;
    //ebeam.SetE(escat::E_el(parent.P()));

    //calculate virtual photon
    const auto& p4beam=*(_prodInfo->_ebeam);
    const auto& p4tar=*(_prodInfo->_target);
    const auto& p4scat=_electron->P4();
    
    //    _gamma = *(_prodInfo->_ebeam) - _electron->P4();
    _gamma = p4beam-p4scat;
    
    //calculate photon polarisation
    auto epsilon = escat::virtualPhotonPolarisation(p4beam,p4tar,p4scat);
    auto delta = 2*escat::M2_el()/(-_gamma.M2())*(1-epsilon);
    
    _photonPol.SetEpsilon(epsilon);
    _photonPol.SetDelta(delta);
    
    return 1.;
  
  }
  
}
