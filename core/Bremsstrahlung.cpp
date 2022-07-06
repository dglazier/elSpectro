#include "Bremsstrahlung.h"

namespace elSpectro{

  ///////////////////////////////////////////////////
  ///For nuclei breakup into 2 parts, 1 of which is used
  /// as initial particle in collision
  Bremsstrahlung::Bremsstrahlung():
    DecayModel{{},{22}} //just produce real photon
  {
    

  }

}
