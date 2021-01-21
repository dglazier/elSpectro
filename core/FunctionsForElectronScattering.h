//////////////////////////////////////////////////////////////
/// utiltiy functions for electron scattering kinematics

#pragma once

#include "LorentzVector.h"
#include <Math/VectorUtil.h> //for boosts etc.

#include "TMath.h"
//#include "Math/RotationY.h"
#include <iostream>

namespace elSpectro{

  namespace escat{
    
    using ROOT::Math::VectorUtil::boost;
    using ROOT::Math::VectorUtil::Angle;
    using ROOT::Math::VectorUtil::CosTheta;
    using ROOT::Math::VectorUtil::InvariantMass2;
  
 

    constexpr double M_el(){return 0.00051099895000;}
    constexpr double M2_el(){return M_el()*M_el();}
    constexpr double M_pr(){return 0.93827208816;}
    constexpr double M2_pr(){return M_pr()*M_pr();}
    
    constexpr double Alpha_by2Pi(){return 1./(137*2*TMath::Pi());}
    
    inline double E2_el(double p){return (p*p+M2_el());}
    inline double E2_pr(double p){return (p*p+M2_pr());}
    inline double E_el(double p){return TMath::Sqrt(E2_el(p));}
    inline double E_pr(double p){return TMath::Sqrt(E2_pr(p));}

    inline double P2_el(double en){return (en*en-M2_el());}
    inline double P2_pr(double en){return (en*en-M2_pr());}
    inline double P_el(double en){return TMath::Sqrt(P2_el(en));}
    inline double P_pr(double en){return TMath::Sqrt(P2_pr(en));}

    inline double y(const LorentzVector& e_in,const LorentzVector& p_in,const LorentzVector& e_sc){
      return  1 - e_sc.Dot(p_in)/e_in.Dot(p_in);
    }
    inline double y_fast(const LorentzVector& e_in,const LorentzVector& e_sc){
      return  1-e_sc.T()/e_in.T();
    }
    
    inline  double SinSqTheta(const LorentzVector& v1,const LorentzVector& v2){
      double cosTh=CosTheta(v1,v2);
      return 1-cosTh*cosTh;
    }
    inline  double SinTheta(const LorentzVector& v1,const LorentzVector& v2){
      return TMath::Sqrt(SinTheta(v1,v2));
    }
    inline  double SinSqHalfTheta(double theta){
      double sinhalf= TMath::Sin(theta/2);
      return sinhalf*sinhalf;
    }
    
    inline double Q2_b(const LorentzVector& e_in,const LorentzVector& e_sc){
      return 4*e_in.T()*e_sc.T()*SinSqHalfTheta(Angle(e_in,e_sc));
    }

    inline double CosTheta_Q2EinEsc(double Q2,double Ein,double Esc ){
      //from Q2 = -(in-out)^2 = 2(Ein*Esc- pin*psc*cos(th) -m^2)
      return  (-0.5*Q2 - M2_el() + Ein*Esc)/P_el(Ein)/P_el(Esc);
     
    }
    inline double virtualPhotonPolarisation(const LorentzVector& e_in,const LorentzVector& p_in,const LorentzVector& e_sc){
      auto yy  = y(e_in,p_in,e_sc);
      return 2*(1-yy)/(1+(1-yy)*(1-yy));
    }
  
    inline double Q2(const LorentzVector& e_in,const LorentzVector& e_sc){
      return InvariantMass2(e_in,e_sc);
    }
    
    inline double Q2min_fast(const LorentzVector& e_in,const LorentzVector& e_sc){
      double y = y_fast(e_in,e_sc);
      double m1y=(1-y);
      return M2_el()*y*y/m1y;
      
    }
    inline double Q2min_y(double y){
      double m1y=(1-y);
      return M2_el()*y*y/m1y;
    }
    inline double  PhotonE(const LorentzVector& e_in,const LorentzVector& e_sc){
      return e_in.T()-e_sc.T();
    }
    
    inline double x(const LorentzVector& e_in,const LorentzVector& e_sc){
      return Q2(e_in,e_sc)/(2*M_pr()*PhotonE(e_in,e_sc));
    }

    inline double ymin_Th(double e_in,double Th,double mTar){
      auto cosTh=TMath::Cos(Th);
      auto cosSqTh=cosTh*cosTh;
      auto sinSqTh = 1 - cosSqTh;
      auto Pel=P_el(e_in);
      auto sinSqThby2 = TMath::Sin(Th/2)*TMath::Sin(Th/2);
      auto n = (1+cosTh)*2*Pel*Pel/mTar/e_in* sinSqThby2;
      auto d = 1 + cosTh*TMath::Sqrt(1 - M2_el()/mTar*mTar*sinSqTh)
	+ (1+cosTh)*2*e_in/mTar*sinSqThby2;
      return n/d;
	
    }
    //Virtual photon flux from x and y
    //egamma=e_sc = e_in * y
    inline double Q2_xy(double e_in,double xx, double yy){
      return 2 * escat::M_pr() * e_in * yy * xx;
    }
    inline double Q2_cosThy(double e_in,double cosTh, double yy){
      double e_sc=e_in*(1-yy);
      return 2*(e_in*e_sc - P_el(e_in)*P_el(e_sc)*cosTh - M2_el());
    }
    inline double y_cosThQ2(double e_in,double cosTh, double Q2){
      double e_sc = (P_el(e_in)*P_el(e_sc)*cosTh + M2_el() + 0.5*Q2)/e_in;
      double y = 1 - e_sc/e_in ;
      return y;
    }
    inline double CosTh_xy(double e_in,double xx, double yy){
      double e_sc=e_in*(1-yy);     
      return  (e_in*e_sc - 0.5*Q2_xy(e_in,xx,yy) - M2_el() )/P_el(e_in)/P_el(e_sc);
    }
    inline double K_xy(double e_in,double xx, double yy){
      double K= (1-xx) * e_in*yy;
      return K<0 ? 0 : K; //protect -ve
    }
    inline double L_xy(double e_in,double xx, double yy){
      double m1y=1-yy;
      double L= ( (1 + m1y*m1y) / yy)  - ( 2 * M2_el() * yy ) / Q2_xy(e_in,yy,xx);
       return L<0 ? 0 : L; //protect -ve
    }
    inline double KLbyE_xy(double e_in,double xx, double yy){
      double m1y=1-yy;
      return (1 + m1y*m1y)  * (1-xx)  -  ( 2 * M2_el() * yy ) / Q2_xy(e_in,yy,xx) *yy; 
    }
    
    inline double K_W2(double W2){
      double K= W2-M2_pr();
      K/=M_pr();
      K/=2;
      return K<0 ? 0 : K; //protect -ve
    }
    /*
      inline double KLbyE_Q2y(double e_in,double Q2, double yy){
      //      double y = (W2 - M2_pr() + Q2)/2/M_pr()/e_in;
      if(y<=0) return 0;
      double m1y=1-yy;
      double L= ( (1 + m1y*m1y) / yy)  - ( 2 * M2_el() * yy ) / Q2;
      
      return L<0 ? 0 : K_W2(ww)*L/e_in; //protect -ve

      }
    */
    inline double flux_dxdy(double e_in,double xx, double yy){
      //protect zero values on denominator
      //mimimum escatter = M_el
      auto m1y=(1-yy);
      // std::cout<<"Q2 "<< Q2_xy(e_in,xx, yy)<<" "<<M2_el()*yy*yy/m1y<<" "<< K_xy(e_in,xx,yy)<<" "<<L_xy(e_in,xx,yy)<< " "<<1./yy<<" "<<1./xx<<std::endl;
      return ((1-yy)*e_in>=M_el()) * xx  ==0 ? 0 :
	Alpha_by2Pi() * K_xy(e_in,xx,yy) * L_xy(e_in,xx,yy) / e_in / yy / xx; 
    }
    
    inline double flux_dlnxdlny(double e_in,double ln_xx, double ln_yy){
      if(ln_xx>0||ln_yy>0) return 0;
      //Range of ln_xx [-infinity , 0 ]
      //protect zero values on denominator
      //mimimum escatter = M_el
      // auto m1y=(1-yy);
      double xx = TMath::Exp(ln_xx);
      double yy = TMath::Exp(ln_yy);
      auto m1y=(1-yy);
      return ((1-yy)*e_in>=M_el()) * xx  ==0 ? 0 :
	Alpha_by2Pi() * K_xy(e_in,xx,yy) * L_xy(e_in,xx,yy) / e_in ; 
    }
    /*
      inline double flux_dlnQ2dlnrW(double e_in,double Q2, double yy){
      //Range of ln_xx [-infinity , 0 ]
      //protect zero values on denominator
      //mimimum escatter = M_el
      // auto m1y=(1-yy);
      
 
      
      auto m1y=(1-yy);
         
      return ((1-yy)*e_in>=M_el()) * xx  ==0 ? 0 :
      Alpha_by2Pi() * KLbyE_Q2y(e_in, Q2, yy); 
      }
    */
  
    inline double Frixione(double e0, double y){

    
      double q2_max = -1 * M2_el()*y*y/(1 - y) ;
      //if(q2_max>-0.1)q2_max=-0.1 ;
      double eg=e0*y;
      double q2_min = - 2*eg*M_pr();//add this
    
      auto flux = Alpha_by2Pi() * (2*M2_el()*y*(1/q2_max - 1/q2_min) + (1 + (1 - y)*(1-y))/y * log(q2_min/q2_max));
    
      if(flux==TMath::Infinity()) return 0.;
      if(TMath::IsNaN(flux)) return 0.;
    
      if(flux<0) return 0.;
    
      return flux;
    }

  }
}
