//////////////////////////////////////////////////////////////
///
///Class:		DistConst
///Description:
///             wrapper for Constant, i.e. no actual distribution
/// 

#pragma once

#include "Distribution.h"

namespace elSpectro{

  class DistConst : public Distribution {

    
  public :

  DistConst(double var, double val=1):_var1{var},_val{val}{};
  DistConst(double var1,double var2, double val):_var1{var1},_var2{var2},_val{val}{};
 
    double SampleSingle()  noexcept final {
      return _var1;
    }
    
    dist_pair SamplePair()   noexcept final {return dist_pair{_var1,_var2};} ;

    double CurrentValue() const noexcept final {return _val;}
    double MaxValue() const noexcept final {return _val;}
    double MinValue() const noexcept final {return _val;}

    double GetVar1() const noexcept { return _var1;}
    double GetVar2() const noexcept { return _var2;}

    double GetMinX() const noexcept final{return _var1;}
    double GetMaxX() const noexcept final{return _var1;}

 
    double GetWeightFor(double valX)  {return 1;}
    double GetValueFor(double valX,double valY=0)final  {return _val;}
    
  private:
    //no one should use default constructor
    DistConst()=default;
    
    double _val{0};
    double _var1{0};
    double _var2{0};
     
    ClassDef(elSpectro::DistConst,1); //class Distribution
 

  };

}
