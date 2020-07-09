#pragma once

#include "LorentzVector.h"
#include <Math/LorentzRotation.h>
#include <Math/RotationX.h>
#include <Math/RotationY.h>
#include <Math/RotationZ.h>

namespace genvector{

  using elSpectro::LorentzVector;
  using elSpectro::BetaVector;


  using  ROOT::Math::RotationX;
  using  ROOT::Math::RotationY;
  using  ROOT::Math::RotationZ;
  using  ROOT::Math::LorentzRotation;
  //using  ROOT::Math::VectorUtil::ProjVector;
  
  
  inline void LorentzRotateX(LorentzVector &vec,double angle){
    LorentzRotation rot{RotationX{angle}};
    vec=rot*vec;
  }
  inline void LorentzRotateY(LorentzVector &vec,double angle){
    LorentzRotation rot{RotationY{angle}};
    vec=rot*vec;
  }
  inline void LorentzRotateZ(LorentzVector &vec,double angle){
    LorentzRotation rot{RotationZ{angle}};
    vec=rot*vec;
  }

  inline void SetZAxis(const LorentzVector &zaxis,LorentzVector &vec){
    LorentzRotateY(vec,zaxis.Theta());
  }

}
