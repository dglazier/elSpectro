//////////////////////////////////////////////////////////////
///
///Class:		VectorSDMEDecay
///Description:
///             Calculate intensity based on vector SDME
///             See Schilling and Wolf formalism

#pragma once

#include "SDMEDecay.h"

namespace elSpectro{

 
  class VectorSDMEDecay : public SDMEDecay {

  public:
    
    VectorSDMEDecay()=default;
    //only declaring default constructor
    //so other 5 constructors also defaulted(rule of 5)
    //constructor to decay into particles

    //decay type = 0 for decay to spin 0; =1 for decay to leptons
    VectorSDMEDecay( particle_ptrs , const std::vector<int> pdgs, int decayType );

    // Each model must define its intensity
    double Intensity() const final;
    
    //void PostInit(ReactionInfo* info) final;
    
    short Spin() const final{ return 1;}

    double W0spin0(double cosSqTh,double sin2Th,double cosPh,double sinSqTh,double cos2Ph) const noexcept;
    double W1spin0(double sinSqTh,double cosSqTh,double sin2Th, double cosPh, double cos2Ph) const noexcept;
    double W2spin0(double sin2Th, double sinPh,double sinSqTh, double sin2Ph) const noexcept;
    double W3spin0(double sin2Th, double sinPh,double sinSqTh, double sin2Ph) const noexcept;

    
    double W0lepto(double cosSqTh,double sin2Th,double cosPh,double sinSqTh,double cos2Ph) const noexcept;
    double W1lepto(double sinSqTh,double cosSqTh,double sin2Th, double cosPh, double cos2Ph) const noexcept;
    double W2lepto(double sin2Th, double sinPh,double sinSqTh, double sin2Ph) const noexcept;
    double W3lepto(double sin2Th, double sinPh,double sinSqTh, double sin2Ph) const noexcept;
  
  private:
    int _type=0;
    
   ClassDefOverride(elSpectro::VectorSDMEDecay,1); //class VectorSDMEDecay
    
  };//class VectorSDMEDecay
  inline   double VectorSDMEDecay::W0spin0(double cosSqTh,double sin2Th,double cosPh,double sinSqTh,double cos2Ph) const noexcept{
 
    return   0.5 * (1 -_rho->Re(0,0,0) )
      + 0.5 * (3*_rho->Re(0,0,0) - 1) *cosSqTh
      - TMath::Sqrt2() * _rho->Re(0,1,0)*sin2Th*cosPh
      - _rho->Re(0,1,-1)*sinSqTh*cos2Ph;
  }

  inline   double VectorSDMEDecay::W1spin0(double sinSqTh,double cosSqTh,double sin2Th, double cosPh, double cos2Ph) const noexcept{

    return   _rho->Re(1,1,1) * sinSqTh
      + _rho->Re(1,0,0) * cosSqTh
      - TMath::Sqrt2() * _rho->Re(1,1,0)*sin2Th*cosPh
      - _rho->Re(1,1,-1)*sinSqTh*cos2Ph;
  }
  
  inline   double VectorSDMEDecay::W2spin0(double sin2Th, double sinPh,double sinSqTh, double sin2Ph) const noexcept{

    return TMath::Sqrt2() * _rho->Im(2,1,0)*sin2Th*sinPh
      + _rho->Im(2,1,-1)*sinSqTh*sin2Ph;
  }
  inline   double VectorSDMEDecay::W3spin0(double sin2Th, double sinPh,double sinSqTh, double sin2Ph) const noexcept{

    return TMath::Sqrt2() * _rho->Im(3,1,0)*sin2Th*sinPh
      + _rho->Im(3,1,-1)*sinSqTh*sin2Ph;
  }

  inline   double VectorSDMEDecay::W0lepto(double cosSqTh,double sin2Th,double cosPh,double sinSqTh,double cos2Ph) const noexcept{
 
    return   0.5 * (1 + _rho->Re(0,0,0) )
      - 0.5 * (3*_rho->Re(0,0,0) - 1) *cosSqTh
      + TMath::Sqrt2() * _rho->Re(0,1,0)*sin2Th*cosPh
      + _rho->Re(0,1,-1)*sinSqTh*cos2Ph;
  }
  inline   double VectorSDMEDecay::W1lepto(double sinSqTh,double cosSqTh,double sin2Th, double cosPh, double cos2Ph) const noexcept{

    return   _rho->Re(1,1,1) * (1+cosSqTh)
      + _rho->Re(1,0,0) * sinSqTh
      + TMath::Sqrt2() * _rho->Re(1,1,0)*sin2Th*cosPh
      + _rho->Re(1,1,-1)*sinSqTh*cos2Ph;
  }
  inline   double VectorSDMEDecay::W2lepto(double sin2Th, double sinPh,double sinSqTh, double sin2Ph) const noexcept{

    return -TMath::Sqrt2() * _rho->Im(2,1,0)*sin2Th*sinPh
      - _rho->Im(2,1,-1)*sinSqTh*sin2Ph;
  }
  inline   double VectorSDMEDecay::W3lepto(double sin2Th, double sinPh,double sinSqTh, double sin2Ph) const noexcept{

    return -TMath::Sqrt2() * _rho->Im(3,1,0)*sin2Th*sinPh
      - _rho->Im(3,1,-1)*sinSqTh*sin2Ph;
  }

}//namespace elSpectro
