//////////////////////////////////////////////////////////////
///
///Class:		DistUniform
///Description:
///             wrapper for uniform distribution
/// 

#pragma once

#include "Distribution.h"
#include <TRandom.h>

namespace elSpectro{

  class DistUniform : public Distribution {

    
  public :

  DistUniform(double xmin,double xmax , double val=1):_xmin{xmin},_xmax{xmax},_val{val}{};
  DistUniform(double xmin,double xmax , double ymin,double ymax , double val=1):
    _xmin{xmin},_xmax{xmax},_ymin{xmin},_ymax{xmax},_val{val}{};
 
    double SampleSingle()  noexcept final {
      return _x=gRandom->Uniform(_xmin,_xmax);
    }
    
    dist_pair SamplePair()   noexcept final {
      _x=gRandom->Uniform(_xmin,_xmax);
      _y=gRandom->Uniform(_ymin,_ymax);
      return dist_pair{_x,_y};
    }

    double CurrentValue() const noexcept final {return _val;}
    double MaxValue() const noexcept final {return _val;}
    double MinValue() const noexcept final {return _val;}

  
    double GetMinX() const noexcept final{return _xmin;}
    double GetMaxX() const noexcept final{return _xmax;}
    double GetMinY() const noexcept final{return _ymin;}
    double GetMaxY() const noexcept final{return _ymax;}

    double GetX() const noexcept { return _x;}
    double GetY() const noexcept { return _y;}
 
    double GetWeightFor(double valX)  {return 1;}
    
  private:
    //no one should use default constructor
    DistUniform()=default;
    
    double _val{0};
    double _x{0};
    double _xmax{0};
    double _xmin{0};
    double _y{0};
    double _ymax{0};
    double _ymin{0};
 
    ClassDef(elSpectro::DistUniform,1); //class Distribution
 

  };

}
