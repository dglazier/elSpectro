//////////////////////////////////////////////////////////////
///
///Class:		GenericModelst
///Description:
///             Control behaviour of Particle decay to Particle products
///             Defined by
///             1) elSpectro Dist

#pragma once

#include "DecayModelst.h"
#include "Distribution.h"

namespace elSpectro{

  class GenericModelst : public DecayModelst {

  public:
    
    GenericModelst()=delete;
    //constructor giving jpac amplitude pointer (which we will now own)
    //and decay particles 
    GenericModelst( Distribution* dist, particle_ptrs parts,
		  const std::vector<int> pdgs  );
    
    bool HasAngularDistribution() override{return true; } //I have an angular distribution


    double MatrixElementsSquared_T() const override {
      
      return _dist->GetValueFor(get_W(),get_t());
    }
    
    
  private:

    std::unique_ptr<Distribution> _dist={nullptr};

    ClassDefOverride(elSpectro::GenericModelst,1); //class GenericModelst
    
  };//class GenericModelst

}//namespace elSpectro
