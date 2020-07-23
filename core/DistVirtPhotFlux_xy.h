//////////////////////////////////////////////////////////////
///
///Class:		DistVirtPhotFlux_xy
///Description:
///             Virtual photon flux distribution as a function of x and y
///             Uses accept or reject from log(x,y) flux function

#pragma once

#include "Distribution.h"
#include <TMath.h>
#include <utility>

namespace elSpectro{

  class DistVirtPhotFlux_xy : public Distribution {

    
  public :

    DistVirtPhotFlux_xy(float ebeam,float xmin,float xmax,float ymin,float ymax);
 
    double SampleSingle()   noexcept final {
      return 0;
    }
    
    dist_pair SamplePair()   noexcept final {
      FindWithAcceptReject();
      return _xy;
    }

    double CurrentValue() const noexcept final {return _val;}
    double MaxValue() const noexcept final {return _max_val;}
    double MinValue() const noexcept final {return 0;}

    double GetMinX() const noexcept final{return TMath::Exp( _lnxmin);}
    double GetMaxX() const noexcept final{return TMath::Exp( _lnxmax);}
    double GetMinY() const noexcept final{return TMath::Exp( _lnymin);}
    double GetMaxY() const noexcept final{return TMath::Exp( _lnymax);}

    void FindWithAcceptReject();

  private:
    //no one should use default constructor
    DistVirtPhotFlux_xy()=default;

    dist_pair _xy{0,0};
    double _val{0};

    //limits
    float _ebeam{0};
    float _lnxmin{0};
    float _lnxmax{0};
    float _lnymin{0};
    float _lnymax{0};

    double _max_val{0};

    
    // ClassDef(elSpectro::DistVirtPhotFlux_xy,1); //class DistVirtPhotFlux_xy
 

  };

}
