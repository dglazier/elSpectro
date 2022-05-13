//////////////////////////////////////////////////////////////
///
///Class:		PhotonPolarisationVector
///Description:
///             Representation of photon polarisation
///             Taken from eqn(78)
///             Schilling, Wolf
///             Nuclear Physics B61 (1973) 381-413.
///             First 4 components give real photon case
///             As from Schilling, Seyboth, Wolf
///             Nuclear Physics BI5 (1970) 397-412.


#pragma once

#include <array>
#include <TMath.h>

namespace elSpectro{

  class PhotonPolarisationVector {

    ///For now take  R->defintion of SDMEs as cannot seperate L and T
    //see eqn(91)
  public:

    void SetRealPol(double eps,double phi){
      SetEpsilonPhi(eps,phi);//[1] and [2]
    }
    void SetEpsilonDeltaPhi(double eps,double delta,double phi){
      SetEpsilonPhi(eps,phi);//[1] and [2]
      auto sinphi=TMath::Sin(phi);
      auto cosphi=TMath::Cos(phi);

      _elements[4] = eps+delta;
      auto factor= TMath::Sqrt(2*eps*(1+eps+2*delta));
      _elements[5] = factor*cosphi;
      _elements[6] = factor*sinphi;
      
    }
    void SetEpsilonPhi(double eps,double phi){
      //auto sin2phi=TMath::Sin(2*phi);
      //auto cos2phi=TMath::Sqrt(1- sin2phi*sin2phi);
      _elements[1] = - eps*TMath::Cos(2*phi); 
      _elements[2] = - eps*TMath::Sin(2*phi); 

    }
    
    void SetLongLepton_Epsilon_h(double eps,double h){
      _elements[3] = TMath::Sqrt(1-eps*eps) * h;
    }
    //need to add SetTransLepton_Epsilon_Phi_h for elements [7] and [8]

    void SetEpsilon(double eps){_epsilon =eps;}
    void SetDelta(double del){_delta =del;}
    void SetPhi(double phi){_phi =phi;}
    void SetPolPlane(double phi){_phiPol =phi;}
    void SetRealPhoto(){_realPhoto=true;}

    void Calc(){
      if(_realPhoto)return SetRealPol(_epsilon,_phi-_phiPol);
      //else electron scattering
      SetEpsilonDeltaPhi(_epsilon,_delta,_phi);
    }
    
    double& operator[](int i)  {return _elements[i];}

    double Epsilon()const noexcept{return _epsilon;}
    double Delta()const noexcept{return _delta;}
    double Phi()const noexcept{return _phi;}
  private:

    std::array<double,8> _elements={1,0,0,0,0,0,0,0};
    double _epsilon={0};
    double _delta = {0};
    double _phi={0};
    double _phiPol={0};//linear polarisation plane
    bool _realPhoto={false};
  };
}
