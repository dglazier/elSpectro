#include "GenericModelst.h"
#include "FunctionsForJpac.h"
#include "FunctionsForGenvector.h"
#include <TDatabasePDG.h>

namespace elSpectro{
  ///////////////////////////////////////////////////////
  ///constructor includes subseqent decay of Ngamma* system
  GenericModelst::GenericModelst( Distribution* dist ,
			      particle_ptrs parts, const std::vector<int> pdgs) :
    _dist{dist},
    DecayModelst{ parts, pdgs }
  {
    _name={"GenericModelst"};
  }
}
