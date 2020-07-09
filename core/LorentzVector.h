/////////////////////////////////////////////////////////
///
///Class:		LorentzVector
///Description:
///             1) Typedef to ROOT::Math::LorentzVector


#pragma once

#include <Math/Vector4D.h> //For XyZTVector
#include <Math/Rotation3D.h> //For XyZTVector
#include <Math/VectorUtil.h> //For when LorentzVector is used in other classes

namespace elSpectro{

  using LorentzVector= ROOT::Math::XYZTVector;
  using BetaVector=ROOT::Math::DisplacementVector3D< ROOT::Math::Cartesian3D<ROOT::Math::Rotation3D::Scalar> >;

}
