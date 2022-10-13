
#include "core/reaction_kinematics.hpp"
#include "core/regge_trajectory.hpp"
#include "core/pseudoscalar_exchange.hpp"
#include "core/pomeron_exchange.hpp"
#include "core/amplitude_sum.hpp"

#include <TMath.h>
#include <TH1D.h>

namespace elSpectro {

  namespace jpacFun{

    void HistProbabilityDistribution_s(jpacPhoto::amplitude* amp, TH1D&  hist);


    double FindMaxOfProbabilityDistribution(jpacPhoto::amplitude* amp,double Wmax);

    TH1D HistFromLargestBins(const TH1D& h1,const TH1D& h2);
  }

}
