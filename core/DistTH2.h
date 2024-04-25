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
 
    double SampleSingle()   noexcept override {
      return 0;
    }
    
    dist_pair SamplePair()   noexcept override {
      RandomXY();
      FindVal();
      return dist_pair{_x,_y};
    } 
    
    double CurrentValue() const noexcept override {return _val;}
    double MaxValue() const noexcept override {return _max_val;}
    double MinValue() const noexcept override {return _min_val;}

    double GetX() const noexcept { return _x;}
    double GetY() const noexcept { return _y;}

    void SetX(double v) const {_x=v;}
    void SetY(double v) const {_y=v;}
    
    double GetMinX() const noexcept override{return _th2.GetXaxis()->GetXmin();}
    double GetMaxX() const noexcept override{return _th2.GetXaxis()->GetXmax();}
    double GetMinY() const noexcept override{return _th2.GetYaxis()->GetXmin();}
    double GetMaxY() const noexcept override{return _th2.GetYaxis()->GetXmax();}

    double GetWeightForXY(double valX,double valY) const {return ((TH2D*)(&_th2))->Interpolate(valX,valY)/_max_val;}
    double GetValueFor(double valX,double valY=0) const override {return GetValueForXY(valX,valY);}
    double GetValueForXY(double valX,double valY) const {return (static_cast<const TH2D*>(&_th2))->Interpolate(valX,valY);}
    
    const TH2& GetTH2() const noexcept {return _th2;}

  protected:
    void SetVal(double val){_val=val;}
    void FindVal(){_val = _th2.GetBinContent(_th2.FindFixBin(_x,_y));}
    virtual void RandomXY(){_th2.GetRandom2(_x,_y);}

  private:
    //no one should use default constructor
    DistTH2()=default;

    TH2D _th2;
    double _val{0};
    mutable double _x{0};
    mutable double _y{0};
    double _max_val{0};
    double _min_val{0};
    
    ClassDef(elSpectro::DistTH2,1); //class Distribution
 

  };

}
