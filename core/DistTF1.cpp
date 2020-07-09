#include "DistTF1.h"
namespace elSpectro{

  DistTF1::DistTF1(const TF1& ff):
    _tf1{ff}{
    _tf1.SetNpx(1E4);
  }
 
}
