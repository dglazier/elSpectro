
#include "reaction_kinematics.hpp"
#include "regge_trajectory.hpp"
#include "amplitudes/pseudoscalar_exchange.hpp"
#include "amplitudes/pomeron_exchange.hpp"
#include "amplitudes/amplitude_sum.hpp"

#include <TMath.h>
#include <TH1D.h>

namespace elSpectro {

  namespace jpacFun{

    void HistProbabilityDistribution_s(jpacPhoto::amplitude* amp, TH1D&  hist);


    double FindMaxOfProbabilityDistribution(jpacPhoto::amplitude* amp,double Wmax);

    
  }

}
