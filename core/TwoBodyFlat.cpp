#include "TwoBodyFlat.h"
#include "FunctionsForKinematics.h"
#include <TMath.h>

namespace elSpectro{


  ////////////////////////////////////////////////////////////////////
  ///Caclulate two body decay from masses and random costh and phi
  ///Return a weight that gives phase-space distribution
  double TwoBodyFlat::Generate(const LorentzVector& parent, const particle_ptrs& products)  {

    _W=parent.M(); //may also be needed by derived classes

    _weight=1;//reset weight
    
    //sample mass of products in case they are from a distribution
    auto m2_a =products[0]->M2();
    auto m2_b =products[1]->M2();
    
    
    if((_W - TMath::Sqrt(m2_a) - TMath::Sqrt(m2_b) ) < 0 ) return 0;//non physical
    auto e_a = (_W*_W + m2_a - m2_b)/(2.0*_W); // E decay product a
    auto p_a = TMath::Sqrt(e_a*e_a - m2_a); // p for both
    // auto e_b = TMath::Sqrt(p_a*p_a + m2_b); // E for decay product b

    auto costh = RandomCosTh();
    auto sinth=TMath::Sqrt(1-costh*costh);
    
    auto phi =  RandomPhi(); // sample phi when z-axis is rotated to simplify  i.e. in BoostToParent
    auto sinphi=TMath::Sin(phi);
    auto cosphi=TMath::Sqrt(1-sinphi*sinphi);
 
    //momentum components in CM frame
    auto x_a = p_a * sinth * cosphi;
    auto y_a = p_a * sinth * sinphi;
    auto z_a = p_a * costh;

 
    //Must make sure scattered products are in the same frame as the parent
    //theta=0 => moving along parent direction
    //i.e. boost vector should only have z component
    //Please note this needs checked for correct rotation
    //Also add the random phi angle...
  
    _a.SetXYZT( x_a, y_a, z_a, e_a);
    BoostToParent(parent,_a);
    products[0]->SetP4(_a);
    products[1]->SetP4( parent - _a );
    
    return _weight; 
  }


}
