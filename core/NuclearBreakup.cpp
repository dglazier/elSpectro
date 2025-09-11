#include "NuclearBreakup.h"

namespace elSpectro{

  ///////////////////////////////////////////////////
  ///For nuclei breakup into 2 parts, 1 of which is used
  /// as initial particle in collision
  NuclearBreakup::NuclearBreakup(int tar_pdg, int spec_pdg ):
    DecayModel{ {},{Particle{tar_pdg},Particle{spec_pdg}}}
   {
    

  }

}
