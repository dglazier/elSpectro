//////////////////////////////////////////////////////////////
///
///Class:		DistTH2Slice
///Description:
///             Fix x sample single random number from y
/// 

#pragma once

#include "DistTH2.h"

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
  
  protected:
    
    void RandomXY() noexcept override;
    
  private:

    std::vector<std::vector<double >> _binIntegrals;
    std::vector<double > _binMax;
    std::vector<double > _binTotal;
    mutable double _current_max=0;
    double _maxIntegral=0;
    double _maxBinVal=0;
    
    ClassDefOverride(elSpectro::DistTH2Slice,1); //class Distribution
 

  };

}
