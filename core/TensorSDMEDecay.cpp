#include "TensorSDMEDecay.h"
#include <TDatabasePDG.h>
#include <array>

namespace elSpectro{

  TensorSDMEDecay::TensorSDMEDecay( particle_ptrs ps, const std::vector<int> pdgs):
     SDMEDecay{ps,pdgs}
  {

    _name={"TensorSDMEDecay"};

  }


  //////////////////////////////////////////////////////////////
  double  TensorSDMEDecay::Intensity() const{
    // std::cout<<" TensorSDMEDecay::Intensity() "<<_photon<<" "<<_meson<<" "<<_baryon<<" "<<_child1<<std::endl;
    //std::cout<<"SDME2 "<< _rho->Val(1,1,1)<<" "<<std::endl;

    //get decay angles in GJ frame
    MomentumVector decayAngles={0,0,1};
    kine::mesonDecayGJ(_photon,_meson,_baryon,_child1,&decayAngles);
    // std::cout<<"VEctor "<<decayAngles.Theta()<<" "<<decayAngles.Phi()<<std::endl;
    //avoid calling trig functions where possible....
    auto theta= decayAngles.Theta();
    auto phi= decayAngles.Phi();
    auto cosTh=TMath::Cos(theta);
    auto cosSqTh = cosTh*cosTh;
    auto sinSqTh= 1 - cosSqTh;
    auto cos2Th = TMath::Cos( 2* theta);
    auto sin2Th = TMath::Sin( 2* theta);
    auto cosPh=TMath::Cos(phi);
    auto sinPh=TMath::Sin(phi);
    auto cos2Ph=TMath::Cos(2*phi);
    auto cos3Ph=TMath::Cos(3*phi);
    auto cos4Ph=TMath::Cos(4*phi);
    auto sin2Ph=TMath::Sin(2*phi);
    auto sin3Ph=TMath::Sin(3*phi);
    auto sin4Ph=TMath::Sin(4*phi);
    auto sinCubeTh=sinSqTh*TMath::Sqrt(sinSqTh);
    auto sinSq2Th =sin2Th*sin2Th;
    /*
    auto cos2Th=cosSqTh - sinSqTh;
    auto sinSq2Th = ( 1 - cos2Th*cos2Th );
    auto sin2Th = TMath::Sin();;
    auto cosPh=TMath::Cos(decayAngles.Phi());
    auto sinPh=TMath::Sqrt(1-cosPh*cosPh);
    auto cos2Ph=TMath::Cos(2*decayAngles.Phi());
    auto cos3Ph=TMath::Cos(3*decayAngles.Phi());
    auto cos4Ph=TMath::Cos(4*decayAngles.Phi());
    auto sin2Ph=TMath::Sqrt(1-cos2Ph*cos2Ph);
    auto sin3Ph=TMath::Sqrt(1-cos3Ph*cos3Ph);
    auto sin4Ph=TMath::Sqrt(1-cos4Ph*cos4Ph);
    */

    std::array<double,8> W={0,0,0,0,0,0,0,0};
    //  Eqn E11 in https://arxiv.org/pdf/2005.01617.pdf "Exclusive tensor meson photoproduction, V. Mathieu"


    W[0] = (
	    1./16*_rho->Re(0,0,0)*(1+3*cos2Th)*(1+3*cos2Th)
	    - 0.75*_rho->Re(0,1,-1)*sinSq2Th*cos2Ph
	    -TMath::Sqrt(3./8)* _rho->Re(0,1,0)*sin2Th*cosPh*(1+3*cos2Th)
	    +0.75*_rho->Re(0,1,1)*sinSq2Th
	    +3*_rho->Re(0,2,-1)*cosTh*sinCubeTh*cos3Ph
	    +3./4*_rho->Re(0,2,-2)*sinSqTh*sinSqTh*cos4Ph
	    +TMath::Sqrt(3./8)*_rho->Re(0,2,0)*sinSqTh*cos2Ph*(1+3*cos2Th)
	    -3*_rho->Re(0,2,1)*cosTh*sinCubeTh*cosPh
	    +0.75*_rho->Re(0,2,2)*sinSqTh*sinSqTh
	    );

    W[1]=(
	    1./16*_rho->Re(1,0,0)*(1+3*cos2Th)*(1+3*cos2Th)
	    - 0.75*_rho->Re(1,1,-1)*sinSq2Th*cos2Ph
	    -TMath::Sqrt(3./8)* _rho->Re(1,1,0)*sin2Th*cosPh*(1+3*cos2Th)
	    +0.75*_rho->Re(1,1,1)*sinSq2Th
	    +3*_rho->Re(1,2,-1)*cosTh*sinCubeTh*cos3Ph
	    +0.75*_rho->Re(1,2,-2)*sinSqTh*sinSqTh*cos4Ph
	    +TMath::Sqrt(3./8)*_rho->Re(1,2,0)*sinSqTh*cos2Ph*(1+3*cos2Th)
	    -3*_rho->Re(1,2,1)*cosTh*sinCubeTh*cosPh
	    +0.75*_rho->Re(1,2,2)*sinSqTh*sinSqTh
	  );
    //Assume 3/4i * _rho_(2,1,-1) = 3/4 *Im{rho2_1-1} as Re=0;
    //i.e. i*_rho->Im(2,1,-1) cancels i in  /4i
    W[2]= (

	   TMath::Sqrt(3./8)* _rho->Im(2,1,0)*sin2Th*sinPh*(1+3*cos2Th)
	   + 0.75*_rho->Im(2,1,-1)*sinSq2Th*sin2Ph
	   - TMath::Sqrt(3./8)*_rho->Im(2,2,0)*sinSqTh*sin2Ph*(1+3*cos2Th)
	   +3*_rho->Im(2,2,1)*cosTh*sinCubeTh*sinPh
	   -3*_rho->Im(2,2,-1)*cosTh*sinCubeTh*sin3Ph
	   -0.75*_rho->Im(2,2,-2)*sinSqTh*sinSqTh*sin4Ph
	   );
    //Assume same form as W[2], not given in paper
    W[3]=  (

	   TMath::Sqrt(3./8)* _rho->Im(3,1,0)*sin2Th*sinPh*(1+3*cos2Th)
	   + 0.75*_rho->Im(3,1,-1)*sinSq2Th*sin2Ph
	   - TMath::Sqrt(3./8)*_rho->Im(3,2,0)*sinSqTh*sin2Ph*(1+3*cos2Th)
	   +3*_rho->Im(3,2,1)*cosTh*sinCubeTh*sinPh
	   -3*_rho->Im(3,2,-1)*cosTh*sinCubeTh*sin3Ph
	   -0.75*_rho->Im(3,2,-2)*sinSqTh*sinSqTh*sin4Ph
  
 	    );
    
    //+ other elctroproduced see eqn(83-85) Schilling and Wolf for Vector

    //equation (82) Schilling and Wolf, cf eqn (29) Schilling Seyboth Wolf
    double result=0;
    _photonPol->Calc(); //epsilon,delta and phi should all now be set
    for(uint alpha=0; alpha < 8; alpha++ ){
      result += ( W[alpha] * (*_photonPol)[alpha] ) ;
    }
    result/=1.05; //make sure < 1
    
    if(result>1.0 || result<0){
      std::cout<<"TensorSDMEDecay invalid result "<< result<<std::endl;
      for(uint alpha=0; alpha < 8; alpha++ ){
	std::cout<<"      "<< W[alpha]<<" "<<(*_photonPol)[alpha]<<" "<<( W[alpha] * (*_photonPol)[alpha])<<std::endl;
	//  result += ( W[alpha] * (*_photonPol)[alpha]) ;
      }
      std::cout<<_rho->Re(0,0,0)<<" 010 "<<_rho->Re(0,1,0)<<" 01-1 "<<_rho->Re(0,1,-1)<<" 111 "<<_rho->Re(1,1,1)<<" 100 "<<_rho->Re(1,0,0)<<" 110 "<<_rho->Re(1,1,0)<<" 11-1 "<<_rho->Re(1,1,-1)<<" 210 "<<_rho->Re(2,1,0)<<" 21-1 "<<_rho->Re(2,1,-1)<<" "<<std::endl;
      std::cout<<"epsilon " <<_photonPol->Epsilon()<<" delta " <<_photonPol->Delta()*TMath::RadToDeg()<<" phi " <<_photonPol->Phi()<<std::endl;
      std::cout<<"CONDITION 1 "<<(_rho->Re(0,0,0)<=1) <<" "<< (_rho->Re(0,0,0)>=0) <<std::endl;
      std::cout<<"CONDITION 2 "<<(TMath::Abs(_rho->Re(0,1,-1))<=0.5*(1-_rho->Re(0,0,0)) )<<std::endl;
      std::cout<<"CONDITION 3 "<<(_rho->Re(0,1,0)*_rho->Re(0,1,0) <=0.25*_rho->Re(0,0,0)*(2-_rho->Re(0,0,0)-_rho->Re(0,1,-1)) )<<std::endl;
      std::cout<<"CONDITION 6 "<<(TMath::Abs(_rho->Re(1,0,0))<=_rho->Re(0,0,0)) <<std::endl;
      std::cout<<"CONDITION 7 "<<(TMath::Abs(_rho->Re(1,1,1))<=0.5*(1-_rho->Re(0,0,0)) )<<std::endl;
      std::cout<<"CONDITION 8 "<<(TMath::Abs(_rho->Re(1,1,-1))<=0.5*(1-_rho->Re(0,0,0))) <<std::endl;
      std::cout<<"CONDITION 9 "<<(TMath::Abs(_rho->Re(1,1,0))<=TMath::Sqrt( 0.5*_rho->Re(0,0,0)*(1-_rho->Re(0,0,0)) )) <<std::endl;
  
    }
    //std::cout<<"Vector weight "<<result<<std::endl;
   
    return result; //max value is already 1, but should check
  }
}
