#include "BremstrPhoton.h"
#include <TMath.h>

namespace elSpectro{

  void BremstrPhoton::PostInit(ReactionInfo* info){
   auto prodInfo = dynamic_cast<ReactionPhotoProd*> (info);
   _photonPol =  prodInfo->_photonPol;
   if(_photonPol==nullptr) return; 
   _photonPol->SetPolPlane(_polPlane);
   _photonPol->SetRealPhoto();
  }

  ////////////////////////////////////////////////////////////////////
  ///Caclulate two body decay from masses and random costh and phi
  ///Return a weight that gives phase-space distribution
  double BremstrPhoton::Generate(const LorentzVector& parent, const particle_ptrs& products)  {
    _weight=1;
    //Get brem e- energy
    auto gammaEnergy = (_bremDist->SampleSingle() )*_ebeam;

      //Get random flat cos theta
    auto costhtar  =  TMath::Cos(parent.Theta()); //assume along e- direction for now
    auto sinthtar      = 0.0;
    //calculate momentum components
    auto px      = gammaEnergy * sinthtar;
    auto py      = 0;
    auto pz      = gammaEnergy * costhtar;

    _photon.SetXYZT(px,py,pz,gammaEnergy);

    //get photon polarisation if required
    if(_linPolDist.get()) _photonPol->SetEpsilon(_linPolDist->GetValueFor(gammaEnergy));
    
    //Must make sure scattered products are in the same frame as the parent
    //theta=0 => moving along parent direction
    //i.e. boost vector should only have z component
    //Please note this needs checked for correct rotation
    //Also add the random phi angle...
    // BoostToParentWithRandPhi(parent,_gamma);
    products[0]->SetP4( _photon );
    // products[1]->SetP4(parent - _photon); //electron
    return _weight; 
  }


}
