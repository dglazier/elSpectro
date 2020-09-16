#pragma once

#include "LorentzVector.h"
#include <TMath.h>
#include <numeric> //for accumulate

namespace elSpectro {
  
  namespace kine {

    using elSpectro::LorentzVector;
    using elSpectro::BetaVector;
    using elSpectro::MomentumVector;
    using ROOT::Math::VectorUtil::boost;
    
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
      //std::cout<<"max cac "<<W<<" "<<TCM<<" "<<std::endl;
      //for (uint n=0; n<nStable; n++) std::cout<< masses[n]<<" ";
      //std::cout<<std::endl;

      for (uint n=1; n<nStable; n++) {
	emmin += masses[n-1];
	emmax += masses[n];
	wtmax *= kine::PDK(emmax, emmin, masses[n]);
      }

      auto massSum= W-TCM;
      auto factor=1.;

      //new and above threshold behaviour
      if(massSum/(nStable-1)<TCM)factor=(TMath::Factorial(nStable-1));
      else factor=(TMath::Factorial(nStable-2));
   
      return wtmax/factor;

    }
    inline double PhaseSpaceWeightMaxFromEquDist(double W, const std::vector<double>& masses){
      //A good approximation of the max weight can be found from
      //assuming each two particle mass combination increases with
      //an equal fraction of TCM
      double TCM= std::accumulate(masses.begin(),masses.end(), W,  std::minus<double>());

      auto const Nt= masses.size();
      std::vector<double> invMass(Nt);
      
      std::vector<Double_t> rno(Nt); //for fractional increases
      rno[0] = 0;
      
      for (uint n=1; n<Nt-1; ++n)  rno[n]=(n)*(1./(Nt-1)); 
      rno[Nt-1] = 1;
      
      double sum=0;
      for (uint n=0; n<Nt; ++n) {
	sum      += masses[n];
	invMass[n] = rno[n]*TCM + sum;
      }
      //now calculate weight
      double wt=1;
   
      for (uint n=0; n<Nt-1; ++n) {
	wt*= PDK(invMass[n+1],invMass[n],masses[n+1]);
      }
      
      return wt;
    } 
    inline double PhaseSpaceFactorDt(double W,double P1,double M3,double M4){
      //this version takes initial momentum P1, allows high Q2
      //if xsect is dt, need change of variables, cosTh -> t and dt/dcosTh = 2*p1*p3 where p1 and p3 in CM frame
      return 1./2/P1/TMath::Sqrt(PDK2(W,M3,M4)); //1 sqrt call
    }
    inline double PhaseSpaceFactorDt(double W,double M1,double M2,double M3,double M4){
      //if xsect is dt, need change of variables, cosTh -> t and dt/dcosTh = 2*p1*p3 where p1 and p3 in CM frame
      return 1./2/TMath::Sqrt(PDK2(W,M1,M2)*PDK2(W,M3,M4)); //1 sqrt call
    }
    inline double FluxPhaseSpaceFactor(const LorentzVector v1,const LorentzVector v2){
      //eqn 47.27 https://pdg.lbl.gov/2019/reviews/rpp2019-rev-kinematics.pdf
      //v1 and v2 are 2 incident particles
      double pdot = v1.T()*v2.T() - (v1.X()*v2.X() + v1.Y()*v2.Y() + v1.Z()*v2.Z());
      return 4*TMath::Sqrt( pdot*pdot  - v1.M2()*v2.M2() );
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
    inline double tFromcosthWP1P3(double costh, double W,double p1,double p3,double M1,double M2,double M3,double M4){
      //with M1<0 PDK can sometimes not return finite number
      double E1 = sqrt(M1*M1 + p1*p1);
      double E3 = sqrt(M3*M3 + p3*p3);
  
      return  M1*M1 + M3*M3  - 2 * ( E1*E3 -p1*p3*costh ); ;
      
    }
    inline double colliderMomentumInRestFrameOf(double m1, double p1, double mrest, double prest){
      //calculate the momentum of 1 with m1 and p1
      //in rest from of rest with mrest and prest
      LorentzVector p(0,0,-p1,sqrt( m1*m1 + p1*p1 ));
      LorentzVector rest(0,0,prest,sqrt( mrest*mrest + prest*prest ));
      auto restBoost=rest.BoostToCM();
      auto p1rest=ROOT::Math::VectorUtil::boost(p,restBoost);
      return p1rest.P();
    }
    
    ////////////////////////////////////////////////////////
    ///z-axis along gamma direction in meson rest frame
    inline void mesonDecayGJ(const LorentzVector* gamma,const LorentzVector* meson,const LorentzVector* baryon,const LorentzVector* d1,MomentumVector* angles){
      auto decBoost=meson->BoostToCM();
      auto decBar=boost(*baryon,decBoost);
      auto decGamma=boost(*gamma,decBoost);
      auto zV=decGamma.Vect().Unit();
      //auto zV=-decBar.Vect().Unit();
      auto yV=decBar.Vect().Cross(decGamma.Vect()).Unit();
      auto xV=yV.Cross(zV).Unit();
    
      LorentzVector decD1=boost(*d1,decBoost);
    
      angles->SetXYZ(decD1.Vect().Dot(xV),decD1.Vect().Dot(yV),decD1.Vect().Dot(zV));
    }
    inline void electroCMDecay(const LorentzVector* CM,const LorentzVector* ein,const LorentzVector* gamma,const LorentzVector* meson,MomentumVector* angles){
    //CM frame defined by e-scattering
       auto CMBoost=CM->BoostToCM();

       //LorentzVectors in CM
       auto CMBeam=boost(*ein,CMBoost);
       auto CMMes=boost(*meson,CMBoost);
       auto CMGamma=boost(*gamma,CMBoost);

       //definition of e- scattering frame
       auto zV=CMGamma.Vect().Unit();
       auto yV=CMGamma.Vect().Cross(CMBeam.Vect()).Unit();
       auto xV=yV.Cross(zV).Unit();
    
       angles->SetXYZ(CMMes.Vect().Dot(xV),CMMes.Vect().Dot(yV),CMMes.Vect().Dot(zV));
    }
    
  }
}
