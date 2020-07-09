#pragma once

#include "LorentzVector.h"
#include <TMath.h>
#include <numeric> //for accumulate

namespace elSpectro {
  namespace kine {

    using elSpectro::LorentzVector;

    inline double PDK2(double a, double b, double c)
    {
      return (a-b-c)*(a+b+c)*(a-b+c)*(a+b-c)/(4*a*a);
    }
    
    inline double PDK(double a, double b, double c){
      return TMath::Sqrt( kine::PDK2(a,b,c) );
    }

    inline double PhaseSpaceWeightMax(double W, const std::vector<double>& masses){

       //subtract masses from W to get TCM
      double TCM= std::accumulate(masses.begin(),masses.end(),W,  std::minus<double>());
      
      double emmax = TCM + masses[0];
      double emmin = 0.;
      double wtmax = 1.;
      const auto nStable = masses.size();
      
      for (uint n=1; n<nStable; n++) {
         emmin += masses[n-1];
         emmax += masses[n];
	 wtmax *= kine::PDK(emmax, emmin, masses[n]);
      }
      return wtmax;

    }

  }
}
