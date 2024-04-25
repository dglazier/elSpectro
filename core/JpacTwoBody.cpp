#include "JpacTwoBody.h"
#include <TDatabasePDG.h>

namespace elSpectro{
  ///////////////////////////////////////////////////////
  ///constructor includes subseqent decay of Ngamma* system
  JpacTwoBody::JpacTwoBody( jpacPhoto::amplitude* amp ,
			      particle_ptrs parts, const std::vector<int> pdgs) :
    _amp{amp},
    TwoBodyProduction{ parts, pdgs }
  {
    _name={"JpacTwoBody"};

  }
 
 
}
