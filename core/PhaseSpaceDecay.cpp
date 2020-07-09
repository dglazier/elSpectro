#include "PhaseSpaceDecay.h"


namespace elSpectro{

  PhaseSpaceDecay::PhaseSpaceDecay( particle_ptrs ps, const std::vector<int> pdgs):
     DecayModel{ps,pdgs}
  {

    _name={"PhaseSpaceDecay"};

  }

}
