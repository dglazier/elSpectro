#include "DistTH1.h"
namespace elSpectro{

  DistTH1::DistTH1(const TH1D& ff):
    _th1{ff}{
    _max_val = ff.GetMaximum();
    _min_val = ff.GetMinimum();
 
  }
 
}
