#include "VectorSDMEDecay.h"
#include <TDatabasePDG.h>
#include <array>

namespace elSpectro{

  VectorSDMEDecay::VectorSDMEDecay( particle_ptrs ps, const std::vector<int> pdgs):
     DecayModel{ps,pdgs}
  {

    _name={"VectorSDMEDecay"};

  }

  ///////////////////////////////////////////////////////////////
  void VectorSDMEDecay::PostInit(ReactionInfo* info){
    DecayModel::PostInit(info);
    auto prodInfo= dynamic_cast<ReactionPhotoProd*> (info); //I need Reaction info

    _rho=Parent()->GetSDME(); //get pointer to SDME values

    if(prodInfo->_meson!=Parent()->P4ptr() ){
      std::cerr<<"VectorSDMEDecay::PostInit this is not the photoproduced meson! I have pdg # "<<Parent()->Pdg()<<" while the photoproduced meson was "<<std::endl;
      exit(0);
    }
    _meson = prodInfo->_meson;
    _baryon = prodInfo->_baryon;
    _photon = prodInfo->_photon;
    _photonPol = prodInfo->_photonPol;
    _child1 = Products()[0]->P4ptr();
  }
  //////////////////////////////////////////////////////////////
  double  VectorSDMEDecay::Intensity() const{
    // std::cout<<" VectorSDMEDecay::Intensity() "<<_photon<<" "<<_meson<<" "<<_baryon<<" "<<_child1<<std::endl;
    //std::cout<<"SDME2 "<< _rho->Val(1,1,1)<<" "<<std::endl;

    //get decay angles in GJ frame
    MomentumVector decayAngles={0,0,1};
    kine::mesonDecayGJ(_photon,_meson,_baryon,_child1,&decayAngles);
    // std::cout<<"VEctor "<<decayAngles.Theta()<<" "<<decayAngles.Phi()<<std::endl;
    //avoid calling trig functions where possible....
    auto cosTh=TMath::Cos(decayAngles.Theta());
    auto cosSqTh = cosTh*cosTh;
    auto sinSqTh= 1 - cosSqTh;
    auto cos2Th=cosSqTh - sinSqTh;
    auto sin2Th = TMath::Sqrt( 1 - cos2Th*cos2Th );
    auto cosPh=TMath::Cos(decayAngles.Phi());
    auto sinPh=TMath::Sqrt(1-cosPh*cosPh);
    auto cos2Ph=TMath::Cos(2*decayAngles.Phi());
    auto sin2Ph=TMath::Sqrt(1-cosPh*cosPh);
    
    //perhaps faster (Cache friendly) to get a copy of rho
    //auto rho = *_rho; and _photonPol ?

    std::array<double,8> W={0,0,0,0,0,0,0,0};
    //  eqn (31) Schilling,Seyboth and Wolf + factor const_3by4pi()* 


    W[0] = (
	       0.5 * (1 -_rho->Re(0,0,0) )
	       + 0.5 * (3*_rho->Re(0,0,0) - 1) *cosSqTh
	       - TMath::Sqrt2() * _rho->Re(0,1,0)*sin2Th*cosPh
	       - _rho->Re(0,1,-1)*sinSqTh*cos2Ph);

    W[1]=(
	     _rho->Re(1,1,1) * sinSqTh
	     + _rho->Re(1,0,0) * cosSqTh
	     - TMath::Sqrt2() * _rho->Re(1,1,0)*sin2Th*cosPh
	     - _rho->Re(1,1,-1)*sinSqTh*cos2Ph
	     );
    W[2]= (
	      TMath::Sqrt2() * _rho->Im(2,1,0)*sin2Th*sinPh
	      - _rho->Im(2,1,-1)*sinSqTh*sin2Ph
	      );
    
    W[3]=  (
	    TMath::Sqrt2() * _rho->Im(3,1,0)*sin2Th*sinPh
	    - _rho->Im(3,1,-1)*sinSqTh*sin2Ph
	    );
    
    //+ other elctroproduced see eqn(83-85) Schilling and Wolf

    
    //equation (82) Schilling and Wolf, cf eqn (29) Schilling Seyboth Wolf
    double result=0;
    _photonPol->Calc(); //epsilon,delta and phi should all now be set
    for(uint alpha=0; alpha < 8; alpha++ ){
      result += ( W[alpha] * (*_photonPol)[alpha] ) ;
    }
    result/=1.05; //max seems to be slightly>1 should check this
    
    if(result>1.0 || result<0){
      std::cout<<"VectorSDMEDecay invalid result "<< result<<std::endl;
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
