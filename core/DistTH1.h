//////////////////////////////////////////////////////////////
///
///Class:		DistTF1
///Description:
///             wrapper for TH1 distributions
/// 

#pragma once

#include "Distribution.h"
#include <TH1.h>

namespace elSpectro{

  class DistTH1 : public Distribution {

    
  public :

    DistTH1(const TH1D& ff);
 
    double SampleSingle()  noexcept final {
      _x=_th1.GetRandom();
      _val=_th1.Interpolate(_x);
      return _x;
    }
    
    dist_pair SamplePair()   noexcept final {return dist_pair{0,0};} ;

    double CurrentValue() const noexcept final {return _val;}
    double MaxValue() const noexcept final {return _max_val;}
    double MinValue() const noexcept final {return _min_val;}

    double GetX() const noexcept { return _x;}
    
    double GetMinX() const noexcept final{return _th1.GetXaxis()->GetXmin();}
    double GetMaxX() const noexcept final{return _th1.GetXaxis()->GetXmax();}

    //  double GetWeightFor(double valX)  {return  (static_cast<TH1D*>(&_th1))->GetBinContent((static_cast<TH1D*>(&_th1))->FindBin(valX))/_max_val;}
    double GetWeightFor(double valX)  {return (static_cast<TH1D*>(&_th1))->Interpolate(valX)/_max_val;}
    double GetValueFor(double valX,double valY=0) final  {return (static_cast<TH1D*>(&_th1))->Interpolate(valX);}
    
    const TH1& GetTH1() const noexcept {return _th1;}
    
  private:
    //no one should use default constructor
    DistTH1()=default;

    TH1D _th1;
    double _val{0};
    double _x{0};
    double _max_val{0};
    double _min_val{0};
    
    ClassDef(elSpectro::DistTH1,1); //class Distribution
 

  };

}
