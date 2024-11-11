#include "JpacTwoBody.h"
#include <TDatabasePDG.h>

namespace elSpectro{
  ///////////////////////////////////////////////////////
  ///constructor includes subseqent decay of Ngamma* system
  JpacTwoBody::JpacTwoBody( jpacAmp_ptr amp ,const decaying_objs& decs, const particle_objs& stables):TwoBodyProduction( decs,stables ),_amp(amp)
  {
    _name={"JpacTwoBody"};

  }
 
 
}
