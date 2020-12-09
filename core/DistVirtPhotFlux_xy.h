//////////////////////////////////////////////////////////////
///
///Class:		DistVirtPhotFlux_xy
///Description:
///             Virtual photon flux distribution as a function of x and y
///             Uses accept or reject from log(x,y) flux function

#pragma once

#include "Distribution.h"
#include "FunctionsForElectronScattering.h"
#include <TMath.h>
#include <utility>

namespace elSpectro{

  class DistVirtPhotFlux_xy : public Distribution {

    
  public :

    //DistVirtPhotFlux_xy(double ebeam,float xmin,float xmax,float ymin,float ymax);
    DistVirtPhotFlux_xy(double eb, double mion, double Wmin);
 
    double SampleSingle()   noexcept final {
      return 0;
    }
    
    dist_pair SamplePair()   noexcept final {
      _forIntegral==false ? FindWithAcceptReject() : FindFlat();
      return _xy;
    }

    double CurrentValue() const noexcept final {return _val;}
    double MaxValue() const noexcept final {return _max_val;}
    double MinValue() const noexcept final {return 0;}

    double GetMinY() const noexcept final{return TMath::Exp( _lnymin);}
    double GetMaxY() const noexcept final{return TMath::Exp( _lnymax);}
    double GetMinX() const noexcept final{return 0;}
    double GetMaxX() const noexcept final{return 1;}

    void FindWithAcceptReject();
    void FindFlat();
    
    void SetElecE(double ee){_ebeam=ee;}
    void SetM(double m){_mTar=m;}
    void SetWThreshold(double Wmin);

    void SetQ2min(double val){_requestQ2min=val;}
    void SetQ2max(double val){_requestQ2max=val;}
    void SetThmin(double val){_requestThmin=val;_requestCosThmin=TMath::Cos(val);}
    void SetThmax(double val){_requestThmax=val;_requestCosThmax=TMath::Cos(val);}
    void SetXmin(double val){_requestXmin=val;}
    void SetXmax(double val){_requestXmax=val;}

    void SetYmax(double val){_requestYmax=val;}
    void SetYmin(double val){_requestYmin=val;}


    void ForIntegrate(bool integ){_forIntegral=integ;}
  protected:
    double XMin(double y);
    double XMax(double y);
    
  private:
    //no one should use default constructor
    DistVirtPhotFlux_xy()=default;

    dist_pair _xy{0,0};
    double _val{0};

    //limits
    double _ebeam={0};
    double _mTar={0};
    double _lnymin={0};
    double _lnymax={0};
    double _Wthresh2={0};
    double _max_val={0};
    double _requestQ2min={0};
    double _requestQ2max={0};
    double _requestCosThmin={-1};
    double _requestCosThmax={-1};
    double _requestThmin={0};
    double _requestThmax={0};
    double _maxPossiblexRange={1};
    
    //no Ymin as set by threshold
    double _requestYmax={0};
    double _requestYmin={0};
    double _requestXmin={0};
    double _requestXmax={1};

    bool _forIntegral=false;
    
    ClassDef(elSpectro::DistVirtPhotFlux_xy,1); //class DistVirtPhotFlux_xy
 

  };
  
  inline  double DistVirtPhotFlux_xy::XMin(double y){
      double r = 2*_mTar*_ebeam*y;
      double Q2min = escat::M2_el()*y*y/(1-y);
      // double Q2max = r + _mTar*_mTar - _Wthresh2;
      double avail_xmin =escat::M2_el()*y/(2*_mTar*_ebeam)/(1-y);
       if( Q2min<_requestQ2min)
	avail_xmin = _requestQ2min/r;
       if(_requestThmin>0){
	 auto Q2fromTh=escat::Q2_cosThy(_ebeam,_requestCosThmin,y);
	 //if limit on P, thre will be a limit on Q2 for a given y
	 if(Q2fromTh>Q2min)
	   avail_xmin = Q2fromTh/r;
       }
      if(_requestXmin> avail_xmin) avail_xmin=_requestXmin;
      return avail_xmin;
    }
   inline  double DistVirtPhotFlux_xy::XMax(double y){
      double r = 2*_mTar*_ebeam*y;
      double Q2max = r + _mTar*_mTar - _Wthresh2;
      double avail_xmax = 1 + (_mTar*_mTar - _Wthresh2 )/r;

      if( _requestQ2max>0)
	if( Q2max>_requestQ2max)//take the smallest "max"
	  avail_xmax = _requestQ2max/r;
      if(_requestThmax>0){
	auto Q2fromTh=escat::Q2_cosThy(_ebeam,_requestCosThmax,y);
	if(Q2fromTh<Q2max)
	  avail_xmax =Q2fromTh/r;
      }
      if(_requestXmax < avail_xmax) avail_xmax=_requestXmax;
      return avail_xmax;
    }

}
