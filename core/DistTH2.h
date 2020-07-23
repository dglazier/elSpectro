//////////////////////////////////////////////////////////////
///
///Class:		DistTF1
///Description:
///             wrapper for TH1 distributions
/// 

#pragma once

#include "Distribution.h"
#include <TH2.h>

namespace elSpectro{

  class DistTH2 : public Distribution {

    
  public :

    DistTH2(const TH2D& ff);
 
    double SampleSingle()   noexcept final {
      return 0;
    }
    
    dist_pair SamplePair()   noexcept final {
      _th2.GetRandom2(_x,_y);
      _val = _th2.Interpolate(_x,_y);
      return dist_pair{_x,_y};
    } 

    double CurrentValue() const noexcept final {return _val;}
    double MaxValue() const noexcept final {return _max_val;}
    double MinValue() const noexcept final {return _min_val;}

    double GetX() const noexcept { return _x;}
    double GetY() const noexcept { return _x;}

    //void SetX(double v) {_x=v;}
    //void SetY(double v) {_y=v;}
    
    double GetMinX() const noexcept final{return _th2.GetXaxis()->GetXmin();}
    double GetMaxX() const noexcept final{return _th2.GetXaxis()->GetXmax();}
    double GetMinY() const noexcept final{return _th2.GetYaxis()->GetXmin();}
    double GetMaxY() const noexcept final{return _th2.GetYaxis()->GetXmax();}

    double GetWeightForXY(double valX,double valY) const {return _th2.Interpolate(valX,valY)/_max_val;}
    
    const TH2& GetTH2() const noexcept {return _th2;}
    
  private:
    //no one should use default constructor
    DistTH2()=default;

    TH2D _th2;
    double _val{0};
    double _x{0};
    double _y{0};
    double _max_val{0};
    double _min_val{0};
    
    ClassDef(elSpectro::DistTH2,1); //class Distribution
 

  };

}
