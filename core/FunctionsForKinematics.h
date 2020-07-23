#pragma once

#include "LorentzVector.h"
#include <TMath.h>
#include <numeric> //for accumulate

namespace elSpectro {
  
  namespace kine {

    using elSpectro::LorentzVector;

    constexpr double M_pr(){return 0.93827208816;}

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

    inline double tmin(double W,double Mx,double Mg,double Mt,double Mr){

      return ( (Mr*Mr-Mx*Mx-Mt*Mt)/2/W )* ( (Mr*Mr-Mx*Mx-Mt*Mt)/2/W)
	- ( PDK(W,Mg,Mt)-PDK(W,Mx,Mr) )*( PDK(W,Mg,Mt)-PDK(W,Mx,Mr) );
    }
    
    inline double t0(double W,double M1,double M2,double M3,double M4){
      double p1 = PDK(W,M1,M2);
      double p3 = PDK(W,M3,M4);
      
      double E1 = sqrt(M1*M1 + p1*p1);
      double E3 = sqrt(M3*M3 + p3*p3);
      return M1*M1 + M3*M3  - 2 * ( E1*E3 -p1*p3 ); 

  }
    inline double tmax(double W,double M1,double M2,double M3,double M4){
      
      double p1 = PDK(W,M1,M2);
      double p3 = PDK(W,M3,M4);
      double tmax = t0(W,M1,M2,M3,M4) - 4*p1*p3;
      return  t0(W,M1,M2,M3,M4) - 4*p1*p3 ; 
      
    }
    inline double costhFromt(double t, double W,double M1,double M2,double M3,double M4){
      double p1 = PDK(W,M1,M2);
      double p3 = PDK(W,M3,M4);
      
      double E1 = sqrt(M1*M1 + p1*p1);
      double E3 = sqrt(M3*M3 + p3*p3);

      double t0 =  M1*M1 + M3*M3  - 2 * ( E1*E3 -p1*p3 ); 

      return 1 - (t0-t)/2/p1/p3;
    }
   inline double tFromcosthW(double costh, double W,double M1,double M2,double M3,double M4){
      double p1 = PDK(W,M1,M2);
      double p3 = PDK(W,M3,M4);
      
      double E1 = sqrt(M1*M1 + p1*p1);
      double E3 = sqrt(M3*M3 + p3*p3);

  
      return  M1*M1 + M3*M3  - 2 * ( E1*E3 -p1*p3*costh ); ;
      
    }
    
  }
}
