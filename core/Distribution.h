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
    virtual dist_pair SamplePair()  noexcept = 0 ;

    virtual double CurrentValue() const noexcept=0;

    virtual double GetMinX() const noexcept = 0 ;
    virtual double GetMaxX() const noexcept = 0 ;
    virtual double GetMinY() const noexcept {return 0;}
    virtual double GetMaxY() const noexcept {return 0;}
    
  private:
    
    //   ClassDef(elSpectro::Distribution,1); //class Distribution
 
  };
  
  using dist_uptr=std::unique_ptr<Distribution>;
}
