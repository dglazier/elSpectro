#include "QuasiFreeNucleon.h"
#include <TMath.h>

namespace elSpectro{


  ////////////////////////////////////////////////////////////////////
  ///Caclulate two body decay from masses and random costh and phi
  ///Return a weight that gives phase-space distribution
  double QuasiFreeNucleon::Generate(const LorentzVector& parent, const particle_ptrs& products)  {
    _weight=1;
    //Get Fermi momentum sampled from distribution given in input file
    auto ptar = _fermiDist->SampleSingle();
    //Get random flat cos theta
    auto costhtar   = RandomCosTh();
    auto sinthtar      = TMath::Sqrt( 1-costhtar*costhtar );
    //Get random flat phi
    //  done in boost
    //  phtar      = RandomPhi();
    //calculate momentum components
    auto pxtar      = ptar * sinthtar;
    auto pytar      = 0;
    auto pztar      = ptar * costhtar;
  
    // Force spectator on mass shell
    //get its mass from input data and PDG database
    auto smass2=products[1]->M2();
    auto Espec = TMath::Sqrt(ptar*ptar + smass2);
    //set spectator 4 momentum (opposite p to quasi target)
    _spectator.SetXYZT(-pxtar,-pytar,-pztar, Espec);
  
    //calculate quasi target energy from energy conservation
    //This will give offshell mass
    auto Etar  = parent.M() - Espec;
    _nucleon.SetXYZT(pxtar,pytar,pztar, Etar);


    //Must make sure scattered products are in the same frame as the parent
    //theta=0 => moving along parent direction
    //i.e. boost vector should only have z component
    //Please note this needs checked for correct rotation
    //Also add the random phi angle...
    //  std::cout<<"spectator "<<_spectator<<" "<<_spectator.M()<<std::endl;
    //std::cout<<"nucleon "<<_nucleon<<" "<<_nucleon.M()<<std::endl;
    //_a.SetXYZT( x_a, y_a, z_a, e_a);
    BoostToParent(parent,_spectator);
    //std::cout<<"spectator boosted "<<_spectator<<" "<<_spectator.M()<<std::endl;
    products[1]->SetP4(_spectator);
    products[0]->SetP4( _nucleon= parent - _spectator );
    //std::cout<<"nucleon boosted"<<_nucleon<<" "<<_nucleon.M()<<std::endl;
    if(_nucleon.M()<0.01) _weight=0; //fix minimum offshell mass to 10MeV...
    return _weight; 
  }


}
