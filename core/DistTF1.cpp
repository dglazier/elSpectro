#include "DistTF1.h"
#include "TH1.h"
namespace elSpectro{

  DistTF1::DistTF1(const TF1& ff):
    _tf1{ff}{
    _max_val = _tf1.GetMaximum();
    _min_val = _tf1.GetMinimum();
    _tf1.SetNpx(1E4);
  }
  double DistTF1::GetMinX() const noexcept{

    Int_t ib=0;
    for(ib=1;ib<_tf1.GetHistogram()->GetNbinsX();ib++){
     auto xval= _tf1.GetHistogram()->GetBinCenter(ib);
     auto fval=_tf1.Eval(xval);
     if(fval>0)break;
    }
    
    return _tf1.GetHistogram()->GetBinLowEdge(ib);
  }
  double DistTF1::GetMaxX() const noexcept{
    Int_t ib=0;
    for(ib=_tf1.GetHistogram()->GetNbinsX(); ib>0 ;ib--){
     auto xval= _tf1.GetHistogram()->GetBinCenter(ib);
     auto fval=_tf1.Eval(xval);
     if(fval>0)break;
    }
    double max=_tf1.GetHistogram()->GetBinLowEdge(ib)+_tf1.GetHistogram()->GetBinWidth(ib);
    return max;

  }

}
