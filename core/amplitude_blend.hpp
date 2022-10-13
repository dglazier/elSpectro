// Class to sum up any number of generic amplitudes and build observables.
// Amplitudes are loaded up in a vector and summed incoherently
//
// Author:       Daniel Winney (2020)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------
#pragma once

#include "core/amplitude.hpp"

// ---------------------------------------------------------------------------
// The amplitude_blend class can take a vector of the above amplitude objects
// and add them together to get observables!
// ---------------------------------------------------------------------------

namespace jpacPhoto
{
  class amplitude_blend : public amplitude
  {
  private:
    // Store a vector of all the amplitudes you want to sum incoherently
    amplitude* _amp_low={nullptr};
    amplitude* _amp_high={nullptr};

    double _low_max = {0};
    double _high_min = {0};
    
  public:
    // Empty constructor
    amplitude_blend(reaction_kinematics * xkinem, amplitude* amp_low, double lowTo, amplitude* amp_high,double highFrom,std::string identifer = "amplitude_blend")
      : amplitude(xkinem, identifer),
	_amp_low{amp_low},_amp_high{amp_high},
	_low_max{lowTo},_high_min{highFrom}
    {
      if(lowTo>highFrom){
	std::cerr<<"Error : amplitude_blend low limit "<<lowTo<<" greater then high"<<highFrom<<std::endl;
	exit(0);
      }
      _isSum = true; //kind of. Actually for cache purposes.
    };

    //assume low amp allowedJP, 
    // inline std::vector<std::array<int,2>> allowedJP()
    // {
    //   //return _amp_low->allowedJP();
    //     return _amp_low->allowed_meson_JP();
    // };

    //take low amp allowed JP
    std::vector<std::array<int,2>> allowed_meson_JP() override{ return _amp_low->allowed_meson_JP();}
    std::vector<std::array<int,2>> allowed_baryon_JP() override{
          return _amp_low->allowed_baryon_JP();
    }

    helicity_channel helicity_CM_frame() override {return _amp_low->helicity_CM_frame();}
    // TODO: Add a set_params which timesi in one vector and allocates approriaten number of
    // params to each sub amplitude

    // Evaluate the sum for given set of helicites, energy, and cos
    std::complex<double> helicity_amplitude(std::array<int, 4> helicities, double s, double t) override;
    
   };
}

