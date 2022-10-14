// Class to sum up any number of generic amplitudes and build observables.
// Amplitudes are loaded up in a vector and summed incoherently
//
// Author:       Daniel Winney (2020)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#include "amplitude_blend.hpp"

// Evaluate the sum for given set of helicites, and mandelstam invariant s and t
std::complex<double> jpacPhoto::amplitude_blend::helicity_amplitude(std::array<int, 4> helicities, double s, double t)
{
   int index = find_helicity(helicities, _kinematics->_jp[0], _kinematics->_mB);
 
    if(s<_low_max){
      _amp_low->check_cache(s, t);
      return _amp_low->_cached_helicity_amplitude[index];
    }
    else if(s>_high_min){
      _amp_high->check_cache(s, t);
      return _amp_high->_cached_helicity_amplitude[index];
    }
    else{
    
    _amp_low->check_cache(s, t);
    _amp_high->check_cache(s, t);

    double blend_range= _high_min - _low_max;
    double fraction_low = 1 - (s-_low_max)/blend_range;

    auto low_contrib = _amp_low->_cached_helicity_amplitude[index] * fraction_low;
    auto high_contrib = _amp_high->_cached_helicity_amplitude[index] * (1 - fraction_low);
    
    return   low_contrib + high_contrib ;
    }
};
