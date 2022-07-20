#include "DistTH1.h"
namespace elSpectro{

  DistTH1::DistTH1(const TH1D& ff):
    _th1{ff}{
    _max_val = ff.GetMaximum();
    _min_val = ff.GetMinimum();
 
  }
  double DistTH1::GetMinX() const noexcept{

    Int_t ib=0;
    for(ib=1;ib<_th1.GetNbinsX();ib++){
     auto fval=_th1.GetBinContent(ib);
     if(fval>0)break;
    }
    
    return _th1.GetBinLowEdge(ib);
  }
  double DistTH1::GetMaxX() const noexcept{
    Int_t ib=0;
    for(ib=_th1.GetNbinsX(); ib>0 ;ib--){
     auto fval=_th1.GetBinContent(ib);
     if(fval>0)break;
    }
    
    return _th1.GetBinLowEdge(ib)+_th1.GetBinWidth(ib);

  }

}
