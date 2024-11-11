//////////////////////////////////////////////////////////////
///
///Class:		DistYGivenX
///Description:
///             Fix x sample single random number from y
/// 

#pragma once

#include "Distribution.h"
#include <Rtypes.h>
#include <TMath.h>
#include <TRandom.h>
#include <TH1D.h>

namespace elSpectro{

  class DistYGivenX : public Distribution {

    using bins = std::vector<double>;
    
  public :

    DistYGivenX() = default;
    DistYGivenX(const bins& x,const std::vector<bins>& y,const std::vector<bins>& vals);
 
    double SampleSingle()   noexcept override {
      //x is fixed then sample y
      RandomXY();
      return GetY();
    }
    dist_pair SamplePair()   noexcept override {
      //!!!!!!This is not a correct algorithm
      //Really need to integrate over y bins
      //then make CDF for x to select a random X
      SetX(gRandom->Uniform(GetMinX(),GetMaxX()));
      RandomXY();
      return dist_pair{_x,_y};
    } 

    double MaxValue() const noexcept override {return _current_max;}
    //double MaxValue() const noexcept override {return 1.;}
    double MinValue() const noexcept override {return 0;}
    double CurrentValue() const noexcept override {return _val;}

    double GetX() const noexcept { return _x;}
    double GetY() const noexcept { return _y;}

    void SetX(double v) const {
      _x=v;
      _current_xbin = TMath::BinarySearch(_nBinsX,_xbins_lowedges.data(),v);
    }
    void SetY(double v) const {_y=v;}
 
    double GetMinX() const noexcept override{return _xbins_lowedges.front();}
    double GetMaxX() const noexcept override{return _xbins_lowedges.back();}
    double GetMinY() const noexcept override{return _ybins_lowedges[_current_xbin].front();}
    double GetMaxY() const noexcept override{return _ybins_lowedges[_current_xbin].back();}

    double GetValueFor(double valX,double valY=0) const override {
      auto xb = TMath::BinarySearch(_nBinsX,_xbins_lowedges.data(),valX);
      auto yb = TMath::BinarySearch(_ybins_lowedges[xb].size(),_ybins_lowedges[xb].data(),valY);
      return _vals[xb][yb];
    }
    TH1D GetMaxValVersusX(){
      TH1D his{"max_versus_x","max_versus_x",static_cast<Int_t>(_nBinsX),_xbins_lowedges.data()};
      for(auto i=0;i<_nBinsX;++i){
	his.SetBinContent(i+1,_xbins_max[i]);
      }
      return his;
    }
    // double LinearInterpolateAlongX(Double_t valx,UInt_t biny,UInt_t binx1,UInt_t binx2){
    //   // std::cout<<" LinearInterpolateAlongX "<<valx<<" "<<biny<<" "<<binx1<<" "<<binx2<<std::endl;      //Takes the current value of x and returns the value at biny
    //   //interpolated between 2 x bins, one of which must contain valx
    //   auto& th2 = GetTH2();
    //   auto x1 = th2.GetXaxis()->GetBinCenter(binx1);
    //   auto x2 = th2.GetXaxis()->GetBinCenter(binx2);
    //   //fractional x
    //   auto frac_x = (valx-x1)/(x2-x1);
    //   //value at x,y
    //   auto val1 = th2.GetBinContent(binx1,biny);
    //   auto val2 = th2.GetBinContent(binx2,biny);
    //   // std::cout<<" LinearInterpolateAlongX "<<x1<<" "<<x2<<" "<<frac_x<<" "<<val1<<" "<<val2<<std::endl;
    //   return val1+frac_x*(val2-val1);
      
    // }
    double LinearInterpolate(double pos,double pos1,double pos2,double val1,double val2){
      auto frac = (pos-pos1)/(pos2-pos1);
      return pos1+frac*(pos2-pos1);
    }
  protected:
    
    virtual void RandomXY() noexcept;
    void SetVal(double val){_val=val;}

  private:
    TH1D _maxXHist;
    
    mutable double _x{0};
    mutable double _y{0};
    mutable double _val{0};

    bins _xbins; //x bins are fixed for all y
    std::vector<bins> _ybins; //y bins can be different for all x
    std::vector<bins> _vals;  //distribution value, mapped to ybins
    
    std::vector<bins> _xbins_integrals; //CDF
    bins _xbins_lowedges;
    bins _xbins_max; //maximum value of distribution in each x bin
    std::vector<bins> _ybins_lowedges;
    std::vector<bins> _ybins_widths;
    
    mutable double _current_max=0;
    mutable uint _current_xbin=0;
    size_t _nBinsX=0;
    
    ClassDef(elSpectro::DistYGivenX,1); //class Distribution
 

  };

 
}
