#include "DecayModelQ2W.h"
#include "DecayingParticle.h"
#include "FunctionsForElectronScattering.h"

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
  DecayModelQ2W::DecayModelQ2W( double thresh, DecayModel* primary) :
    _threshold{thresh},
    DecayModel{{ new DecayingParticle{-2211,primary} },{11}}
  {
    _name={"DecayModelQ2W_with_primary_decay"};
    Init();
  }
  void DecayModelQ2W::Init(){
    if(Products()[0]->Pdg()==11){
      _gstarNuc=Products()[1]; //-2211
      _electron=Products()[0]; //11 
    }
    else{
      _gstarNuc=Products()[0]; //-2211
      _electron=Products()[1]; //11
    }
    std::cout<<"DecayModelQ2W::DecayModelQ2W "<<std::endl;
    _electron->Print();
    _gstarNuc->Print();

    if(_threshold<MinimumMassPossible() )_threshold=MinimumMassPossible(); 
  }
  const CurrentEventInfo* DecayModelQ2W::Intensity(const CurrentEventInfo* info) const{

    if(CheckThreshold()==false){
      _myInfo._W=0;
      _myInfo._Q2=0;
      _myInfo._weight=0;
      return &_myInfo;
    }
    
    auto parent = ParentVector();
    //  auto parent = _electron->P4() + _gstarNuc->P4();

    //  std::cout<<"DecayModelQ2W::Intensity "<<parent<<" "<<_electron->P4()<<std::endl;
    //note we are in the nuc rest frame, so parent momentum = e- beam momentum
    auto ebeam = parent;
    ebeam.SetE(escat::E_el(parent.P()));
    
    auto gammastar = ebeam - _electron->P4();
    
    //set info for passing on
    _myInfo._W = _gstarNuc->P4().M();
    _myInfo._Q2= -1 * gammastar.M2();
    
    //evaluate model
    if(_myInfo._W > _threshold ) _myInfo._weight = 1;
    else  _myInfo._weight = 0;
   
    // std::cout<<"DecayModelQ2W::Intensity "<<_threshold<<" "<<ebeam<<" "<<gammastar<<" "<<_myInfo._W<<" "<<_myInfo._Q2<<"               "<<_myInfo._weight<<std::endl;

    return &_myInfo;
  }
  
}
