#include "DistTH2.h"

namespace elSpectro{

  DistTH2::DistTH2(const TH2D& ff):
    _th2{ff}{
    _max_val = _th2.GetMaximum();
    _min_val = _th2.GetMinimum();
 
  }
 
}
