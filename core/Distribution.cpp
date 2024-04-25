#include "Distribution.h"

namespace elSpectro{

  double Distribution::Integrate1DX(double xlow,double xhigh) const{
    if(xlow==xhigh) {
      xlow = GetMinX();
      xhigh = GetMaxX();
    }
    
    ROOT::Math::GSLIntegrator ig(ROOT::Math::IntegrationOneDim::kADAPTIVE,
				 ROOT::Math::Integration::kGAUSS61);
    //ROOT::Math::Functor1D wF(*this,&Distribution::GetValueFor);
    auto F = [this](double x){
      return GetValueFor(x);
    };
    ROOT::Math::Functor1D wF(F);
    ig.SetFunction(wF);

    return ig.Integral(xlow,xhigh);
  }
  
}
