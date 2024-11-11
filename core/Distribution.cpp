#include "Distribution.h"

namespace elSpectro{

  double Distribution::Integrate1DX(double xlow,double xhigh) const{
    if((_cacheIntegralLow==xlow) && (_cacheIntegralHigh==xhigh) )
      return _cacheIntegral;

    
    if(xlow==xhigh) {
      xlow = GetMinX();
      xhigh = GetMaxX();
    }
    else{
      xlow = xlow < GetMinX() ? GetMinX() : xlow; //enforce true min
      xhigh = xhigh > GetMaxX() ? GetMaxX() : xhigh; //enforce true max
    }
    
    ROOT::Math::GSLIntegrator ig(ROOT::Math::IntegrationOneDim::kADAPTIVE,
				 ROOT::Math::Integration::kGAUSS61);
    auto F = [this](double x){
      return GetValueFor(x);
    };
    ROOT::Math::Functor1D wF(F);
    ig.SetFunction(wF);

    _cacheIntegralLow =xlow; 
    _cacheIntegralHigh =xhigh;
    auto result = ig.Integral(xlow,xhigh);
    _cacheIntegral= result;
    return result;
  }
  
  double Distribution::Mean1DX(double xlow,double xhigh) const{
    if((_cacheMeanLow==xlow) && (_cacheMeanHigh==xhigh) )
      return _cacheMean;

    if(xlow==xhigh) {
      xlow = GetMinX();
      xhigh = GetMaxX();
    }
    else{
      xlow = xlow < GetMinX() ? GetMinX() : xlow; //enforce true min
      xhigh = xhigh > GetMaxX() ? GetMaxX() : xhigh; //enforce true max
    }
    
    ROOT::Math::GSLIntegrator ig(ROOT::Math::IntegrationOneDim::kADAPTIVE,
				 ROOT::Math::Integration::kGAUSS61);
    auto F = [this](double x){
      return GetValueFor(x) * x;
    };
    ROOT::Math::Functor1D wF(F);
    ig.SetFunction(wF);

    _cacheMeanLow =xlow; 
    _cacheMeanHigh =xhigh;
    auto result = ig.Integral(xlow,xhigh)/Integrate1DX(xlow,xhigh);
    _cacheMean= result;
    return result;
  }
  
}
