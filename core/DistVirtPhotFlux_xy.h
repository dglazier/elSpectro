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

    //DistVirtPhotFlux_xy(double ebeam,float xmin,float xmax,float ymin,float ymax);
    DistVirtPhotFlux_xy(double eb, double mion, double Wmin);
 
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

    double GetMinY() const noexcept final{return TMath::Exp( _lnymin);}
    double GetMaxY() const noexcept final{return TMath::Exp( _lnymax);}
    double GetMinX() const noexcept final{return 0;}
    double GetMaxX() const noexcept final{return 1;}

    void FindWithAcceptReject();
    
    void SetElecE(double ee){_ebeam=ee;}
    void SetM(double m){_mTar=m;}
    void SetWThreshold(double Wmin);
    
  private:
    //no one should use default constructor
    DistVirtPhotFlux_xy()=default;

    dist_pair _xy{0,0};
    double _val{0};

    //limits
    double _ebeam{0};
    double _mTar{0};
    double _lnymin{0};
    double _lnymax{0};
    double _Wthresh2{0};
    double _max_val{0};
    double _requestQ2min{-1};
    double _requestQ2max{1E10};

    
    ClassDef(elSpectro::DistVirtPhotFlux_xy,1); //class DistVirtPhotFlux_xy
 

  };

}
