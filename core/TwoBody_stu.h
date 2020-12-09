//////////////////////////////////////////////////////////////
///
///Class:		TwoBody_stu
///Description:
///            Derive methods for 
///            Generating LorentzVectors for parent->1,2
#pragma once

#include "TwoBodyFlat.h"
#include "FunctionsForKinematics.h"
#include <Math/VectorUtil.h> //for boosts etc.

namespace elSpectro{


  class TwoBody_stu : public TwoBodyFlat {

  public:

    TwoBody_stu(double s,double t,double t_slope,double u,double u_slope);

    double MyRandomCosTh() const noexcept final{
      _weight=1;
      auto randChannel = gRandom->Uniform();
      //  std::cout<<" TwoBody_stu RandomCosTh() "<<randChannel<<std::endl;
      //select s,t or u channel
      double costh=0;
      if(randChannel<=_t_strength)
	costh= CosThFrom_t();
      else if(randChannel<=(_t_strength+_s_strength))
	costh= CosThFrom_s();
      else
	costh= CosThFrom_u();
      
      _weight = CalcWeight();
      return costh;
    }
   
    double CosThFrom_t() const noexcept {
       
      double W = _CM->M();
      double M1=_p1->M();
      // double M1=0;
      double M2=_p2->M();
      double M3=_p3->M();
      double M4=_p4->M();
      
      double P3 = kine::PDK(W,M3,M4);

      //PDK does not always have a valid solution for g*
      //Direct boost of g* into cm rest frame
      auto cmBoost=_CM->BoostToCM();
      auto p1cm=boost(*_p1,cmBoost);
      auto P1=p1cm.P();
    
  
      double E1 = sqrt(M1*M1 + P1*P1);
      double E3 = sqrt(M3*M3 + P3*P3);
     
      //std::cout<<W-M1-M2<<" "<<W-TMath::Abs(M1)-M2<<" "<<kine::PDK2(W,M1,M2)<<" "<<sqrt(-kine::PDK2(W,M1,M2))<<" "<< kine::PDK(W,TMath::Abs(M1),M2)<<std::endl;
        
      /*auto prBoost=_p2->BoostToCM();
      auto p2cm=boost(*_p2,prBoost);
      p2cm=boost(p2cm,cmBoost);
      std::cout<<"P1CM "<<p1cm.M2()<<" "<<p1cm.P2()<<" "<<p1cm.P()<<" "<<p1cm.E()<<" "<<std::endl;
      std::cout<<"P2CM "<<p2cm.M2()<<" "<<p2cm.P2()<<" "<<p2cm.P()<<" "<<p2cm.E()<<" "<<std::endl;
      std::cout<<" check "<<2*_p1->E()*M2<<" "<<_p1->M2()<<std::endl;
      */
      //limits on t
      double tmin =  M1*M1 + M3*M3  - 2 * ( E1*E3 -P1*P3 ); 
      double tmax = tmin - 4*P1*P3 ;
      _t = tmax*2; //start off >tmax
      while( (_t=tmin - gRandom->Exp(1./_t_slope)) < tmax ){}; //tau=1/b0
      
      // while( t < tmax || t > tmin ){t= - gRandom->Exp(1./_t_slope);}; //tau=1/b0
      // t = -0.1+tmin;
      //weight is value/max_value as for Distribution
      //here max value =1 @ 0
      //*= in case any other weighting used
      //_weight *= TMath::Exp( (t) * _t_slope) * _t_strength;
      //_weight *= TMath::Exp( (_t-tmin) * _t_slope) * _t_strength + _strength;

         // if(TMath::IsNaN(tmin)) exit(0);
      return (1 - (tmin - _t)/2/P1/P3); //cos(theta) from t
      
    }
    double CalcWeight() const{
      double W = _CM->M();
      double M1=_p1->M();
      // double M1=0;
      double M2=_p2->M();
      double M3=_p3->M();
      double M4=_p4->M();
      
      double P3 = kine::PDK(W,M3,M4);

      //PDK does not always have a valid solution for g*
      //Direct boost of g* into cm rest frame
      auto cmBoost=_CM->BoostToCM();
      auto p1cm=boost(*_p1,cmBoost);
      auto P1=p1cm.P();
      
  
      double E1 = sqrt(M1*M1 + P1*P1);
      double E3 = sqrt(M3*M3 + P3*P3);
      double tmin =  M1*M1 + M3*M3  - 2 * ( E1*E3 -P1*P3 ); 
      return TMath::Exp( (_t-tmin) * _t_slope) * _t_strength + _s_strength;
    }
    double CosThFrom_u() const noexcept {
       //Needs fixed to return correct cos theta and applu u distribution
       return 0.;
       /*
      double W = _CM->M();
      double M1=_p1->M();
      double M2=_p2->M();
      double M3=_p3->M();
      double M4=_p4->M();
      
      double P4 = kine::PDK(W,M3,M4);
      //PDK does not always have a valid solution for g*
      //Direct boost of g* into cm rest frame
      auto cmBoost=_CM->BoostToCM();
      auto p1cm=boost(*_p1,cmBoost);
      auto P1=p1cm.P();
 
      double E1 = sqrt(M1*M1 + P1*P1);
      double E4 = sqrt(M4*M4 + P4*P4);
     std::cout<<W<<" "<<M1<<" "<<" "<<M2<<" "<<" "<<M3<<" "<<" "<<M4<<" p1 "<<P1<<" "<<P4<<" "<<E1<<" "<<E4<<std::endl;

      //limits on t
      double umax =  M1*M1 + M4*M4  - 2 * ( E1*E4 + P1*P4 ); 
      double umin = umax + 4*P1*P4 ;
      std::cout<<"umin "<<umin<<" "<<umax<<std::endl;
      double u = umax*2; //start off >tmax
      while( (u=umin - gRandom->Exp(1./_u_slope)) < umax ){}; //tau=1/b0
      //u = -0.1;
 
      //weight is value/max_value as for Distribution
      //here max value =1
      _weight *= TMath::Exp( (u-umin) * _u_slope);// * _u_strength;
      
     std::cout<<u<<" u "<<umin<<" "<<umax<<" "<<1 - (u-umax)/2/P1/P4<<" "<<_weight<<std::endl;
      if(TMath::IsNaN(umin)) exit(0);
      return (1 - (u-umax)/2/P1/P4); //cos(theta) from u
       */
    }

    double CosThFrom_s() const noexcept {
      //weight is value/max_value as for Distribution
      //here max value =1
      // _weight *=_s_strength;
      double costh= gRandom->Uniform(-1,1);
      double W = _CM->M();
      double M1=_p1->M();
      double M2=_p2->M();
      double M3=_p3->M();
      double M4=_p4->M();
      
   //PDK does not always have a valid solution for g*
      //Direct boost of g* into cm rest frame
      auto cmBoost=_CM->BoostToCM();
      auto p1cm=boost(*_p1,cmBoost);
      auto P1=p1cm.P();
    
  
      _t=kine::tFromcosthWP1P3(costh, P1, kine::PDK(W,M3,M4), W,M1,M2,M3,M4);
    
      return  costh;
    }

    void PostInit(ReactionInfo* info);

  private:

    double _s_strength={0};
    double _t_strength={0};
    double _t_slope={0};
    double _u_strength={0};
    double _u_slope={0};

    mutable double _t;
  
    LorentzVector* _p1={nullptr};
    LorentzVector* _p2={nullptr};
    LorentzVector* _p3={nullptr};
    LorentzVector* _p4={nullptr};
    LorentzVector* _CM={nullptr};
    
    ClassDef(elSpectro::TwoBody_stu,1); //class DecayVectors
 

  };

}
  
