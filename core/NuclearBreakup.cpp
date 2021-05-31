#include "NuclearBreakup.h"

namespace elSpectro{

  ///////////////////////////////////////////////////
  ///For nuclei breakup into 2 parts, 1 of which is used
  /// as initial particle in collision
  NuclearBreakup::NuclearBreakup(int pdg1, int pdg2 ):
    DecayModel{{},{pdg1,pdg2}}
  {
    

  }

}
