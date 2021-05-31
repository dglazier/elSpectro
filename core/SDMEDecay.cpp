#include "SDMEDecay.h"
#include <TDatabasePDG.h>
#include <array>

namespace elSpectro{

  SDMEDecay::SDMEDecay( particle_ptrs ps, const std::vector<int> pdgs):
     DecayModel{ps,pdgs}
  {

    _name={"SDMEDecay"};

  }

  ///////////////////////////////////////////////////////////////
  void SDMEDecay::PostInit(ReactionInfo* info){
    DecayModel::PostInit(info);
    auto prodInfo= dynamic_cast<ReactionPhotoProd*> (info); //I need Reaction info

    _rho=Parent()->GetSDME(); //get pointer to SDME values

    if(prodInfo->_meson!=Parent()->P4ptr() ){
      std::cerr<<"SDMEDecay::PostInit this is not the photoproduced meson! I have pdg # "<<Parent()->Pdg()<<" while the photoproduced meson was "<<std::endl;
      exit(0);
    }
    _meson = prodInfo->_meson;
    _baryon = prodInfo->_baryon;
    _photon = prodInfo->_photon;
    _photonPol = prodInfo->_photonPol;
    _child1 = Products()[0]->P4ptr();
  }
 
}
