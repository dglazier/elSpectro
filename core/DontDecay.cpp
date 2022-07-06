#include "DontDecay.h"
#include "FunctionsForKinematics.h"
#include <TMath.h>

namespace elSpectro{


  ////////////////////////////////////////////////////////////////////
  ///Caclulate two body decay from masses and random costh and phi
  ///Return a weight that gives phase-space distribution
  double DontDecay::Generate(const LorentzVector& parent, const particle_ptrs& products)  {

    _W=parent.M(); //may also be needed by derived classes

    _weight=1;//reset weight
    
    _a = parent;
    products[0]->SetP4(_a);
    return _weight; 
  }


}
