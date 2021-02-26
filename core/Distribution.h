//////////////////////////////////////////////////////////////
///
///Class:		Distribution
///Description:
///             interface to  sampling distributions wrappers

#pragma once

#include<utility> //for pair
#include<memory> //for unique_ptr

namespace elSpectro{

  using dist_pair = std::pair<double,double>;

  
  class Distribution {

  public:
    
    Distribution()=default;
    
  public :

    virtual double SampleSingle()  noexcept = 0 ;
    virtual double SampleSingle(double xmin,double xmax)  noexcept { return SampleSingle();};//default ignores limits in case not applicable for derived class
    virtual dist_pair SamplePair()  noexcept = 0 ;

    virtual double CurrentValue() const noexcept=0;
    virtual double MaxValue() const noexcept=0;
    virtual double MinValue() const noexcept=0;

    virtual double GetMinX() const noexcept = 0 ;
    virtual double GetMaxX() const noexcept = 0 ;
    virtual double GetMinY() const noexcept {return 0;}
    virtual double GetMaxY() const noexcept {return 0;}

    double GetCurrentWeight() const noexcept { return CurrentValue()/MaxValue();}

    virtual double GetValueFor(double valX,double valY=0)= 0 ;
    
  protected :

 
  private:
    
    //   ClassDef(elSpectro::Distribution,1); //class Distribution
 
  };
  
  using dist_uptr=std::unique_ptr<Distribution>;
}
