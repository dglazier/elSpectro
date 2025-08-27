//////////////////////////////////////////////////////////////
///
///Class:		Distribution
///Description:
///             interface to  sampling distributions wrappers

#pragma once

#include<utility> //for pair
#include<memory> //for unique_ptr
#include <Math/GSLIntegrator.h>
#include <Math/IntegrationTypes.h>
#include <Math/Functor.h>

namespace elSpectro{

  using dist_pair = std::pair<double,double>;

  
  class Distribution {

  public:
    
    Distribution()=default;
    virtual ~Distribution()=default;
    
  public :

    virtual double SampleSingle()  noexcept = 0 ;
    
    virtual double SampleSingle(double xmin,double xmax)  noexcept {std::cerr<<"Distribution::SampleSingle(double xmin,double xmax) not implmented for this derived distribution"<<std::endl;exit(0); return SampleSingle();};//default ignores limits in case not applicable for derived class

    virtual dist_pair SamplePair()  noexcept = 0 ;

    virtual double CurrentValue() const noexcept=0;
    virtual double MaxValue() const noexcept=0;
    virtual double MinValue() const noexcept=0;

    virtual double GetMinX() const noexcept = 0 ;
    virtual double GetMaxX() const noexcept = 0 ;
    virtual double GetMinY() const noexcept {return 0;}
    virtual double GetMaxY() const noexcept {return 0;}

    double GetCurrentWeight() const noexcept { return CurrentValue()/MaxValue();}
    double GetWeightFor(double valX,double valY=0) const {return GetValueFor(valX,valY)/MaxValue();}

    virtual double GetValueFor(double valX,double valY=0) const = 0 ;

    virtual double Integrate1DX(double xlow=0,double xhigh=0) const;
    virtual double Mean1DX(double xlow=0,double xhigh=0) const;

    virtual void Print(){};
  protected :

 
  private:

    mutable double _cacheIntegralLow=-1;
    mutable double _cacheIntegralHigh=-1;
    mutable double _cacheIntegral=0;
    mutable double _cacheMeanLow=-1;
    mutable double _cacheMeanHigh=-1;
    mutable double _cacheMean=0;

    //ClassDef(elSpectro::Distribution,1); //class Distribution
 
  };
  
  using dist_uptr=std::unique_ptr<Distribution>;
}
