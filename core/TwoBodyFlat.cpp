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
    
    //momentum components in CM frame
    //y applied in rotation function
    auto x_a = p_a* sinth;
    auto y_a = 0;
    auto z_a = p_a * costh;
    //auto x_a = p_a * sinth * cosphi;
    //auto y_a = p_a * sinth * sinphi;
    //auto z_a = p_a * costh;

 
    //Rotate out of parent rest frame and boost
    //to frame parent is in
    _a.SetXYZT( x_a, y_a, z_a, e_a);
    BoostToParentWithRandPhi(parent,_a);
    products[0]->SetP4(_a);
    products[1]->SetP4( parent - _a );
    //   std::cout<<"TwoBody "<<products[0]->P4()<<" "<<products[0]->P4().M()<<" "<<products[1]->P4()<<" "<<products[1]->P4().M()<<std::endl;
    return _weight; 
  }


}
