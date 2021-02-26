//////////////////////////////////////////////////////////////
///
///Class:		DistTF1
///Description:
///             wrapper for TF1 distributions
/// 

#pragma once

#include "Distribution.h"
#include <TF1.h>
#include <string>

namespace elSpectro{

  class DistTF1 : public Distribution {

    
  public :

    DistTF1(const TF1& ff);
 
    double SampleSingle()   noexcept final {
      _x=_tf1.GetRandom();
      _val=_tf1.Eval(_x);
      return _x;
    }
    double SampleSingle(double xmin,double xmax)   noexcept final {
      _x=_tf1.GetRandom(xmin,xmax);
      _val=_tf1.Eval(_x);
      return _x;
    }
    
    dist_pair SamplePair()   noexcept final {return dist_pair{0,0};} ;

    double CurrentValue() const noexcept final {return _val;}
    double MaxValue() const noexcept final {return _max_val;}
    double MinValue() const noexcept final {return _min_val;}

    double GetX() const noexcept { return _x;}
    
    double GetMinX() const noexcept final{return _tf1.GetXmin();}
    double GetMaxX() const noexcept final{return _tf1.GetXmax();}

    double GetWeightFor(double valX) const {return _tf1.Eval(valX)/_max_val;}
    double GetValueFor(double valX,double valY=0) final {return _tf1.Eval(valX);}
    
    const TF1& GetTF1() const noexcept {return _tf1;}
    
  private:
    //no one should use default constructor
    DistTF1()=default;

    TF1 _tf1;
    double _val{0};
    double _x{0};
    double _max_val{0};
    double _min_val{0};
    
    ClassDef(elSpectro::DistTF1,1); //class Distribution
 

  };

}
