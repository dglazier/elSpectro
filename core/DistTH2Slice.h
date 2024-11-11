//////////////////////////////////////////////////////////////
///
///Class:		DistTH2Slice
///Description:
///             Fix x sample single random number from y
/// 

#pragma once

#include "DistTH2.h"
#include "DistTH1.h"

namespace elSpectro{

  class DistTH2Slice : public DistTH2 {

    
  public :

    DistTH2Slice(const TH2D& ff);
 
    double SampleSingle()   noexcept override {
      //x is fixed then sample y
      RandomXY();
      return GetY();
    }
    double MaxValue() const noexcept override {return _current_max;}

    double LinearInterpolateAlongX(Double_t valx,UInt_t biny,UInt_t binx1,UInt_t binx2){
      // std::cout<<" LinearInterpolateAlongX "<<valx<<" "<<biny<<" "<<binx1<<" "<<binx2<<std::endl;      //Takes the current value of x and returns the value at biny
      //interpolated between 2 x bins, one of which must contain valx
      auto& th2 = GetTH2();
      auto x1 = th2.GetXaxis()->GetBinCenter(binx1);
      auto x2 = th2.GetXaxis()->GetBinCenter(binx2);
      //fractional x
      auto frac_x = (valx-x1)/(x2-x1);
      //value at x,y
      auto val1 = th2.GetBinContent(binx1,biny);
      auto val2 = th2.GetBinContent(binx2,biny);
      // std::cout<<" LinearInterpolateAlongX "<<x1<<" "<<x2<<" "<<frac_x<<" "<<val1<<" "<<val2<<std::endl;
      return val1+frac_x*(val2-val1);
      
    }
  protected:
    
    void RandomXY() noexcept override;
    
  private:

    DistTH1 _distHighX;
    std::vector<std::vector<double >> _binIntegrals;
    std::vector<double > _binMax;
    std::vector<double > _binTotal;
    mutable double _current_max=0;
    double _maxIntegral=0;
    double _maxBinVal=0;
    
    ClassDefOverride(elSpectro::DistTH2Slice,1); //class Distribution
 

  };

}
