#include "VectorSDMEDecay.h"
#include <TDatabasePDG.h>
#include <array>

namespace elSpectro{

  VectorSDMEDecay::VectorSDMEDecay( particle_ptrs ps, const std::vector<int> pdgs, int decayType):
    SDMEDecay{ps,pdgs},_type{decayType}
  {

    _name={"VectorSDMEDecay"};
    if(pdgs[0]==11||pdgs[0]==13||pdgs[1]==11||pdgs[1]==13){
      if(_type!=1)
	Fatal("VectorSDMEDecay"," using Spin 0 decay equations, but with leptons! If you want to continue to do this , change decaytype to 100, but better to change decay type to 1");
      else if(_type==100) //lets use of spin 0 equations for leptons, Should not really be used
	_type = 0;
    }

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
    auto sin2Ph=TMath::Sin(2*phi);
    
    //perhaps faster (Cache friendly) to get a copy of rho
    //auto rho = *_rho; and _photonPol ?

    std::array<double,8> W={0,0,0,0,0,0,0,0};
    //  eqn (31) Schilling,Seyboth and Wolf + factor const_3by4pi()* 

    if(_type==1){//decay to lepton e.g. e+ e-
      W[0]= W0lepto(cosSqTh,sin2Th,cosPh,sinSqTh,cos2Ph)*3./8/TMath::Pi();
      W[1]= W1lepto(sinSqTh,cosSqTh,sin2Th,cosPh,cos2Ph)*3./8/TMath::Pi();
      W[2]= W2lepto(sin2Th, sinPh,sinSqTh, sin2Ph)*3./8/TMath::Pi();
      W[3]= W3lepto(sin2Th,sinPh,sinSqTh,sin2Ph)*3./8/TMath::Pi();
    }
    else{ //decay to spin0 e.g. pi+ pi-
      W[0]= W0spin0(cosSqTh,sin2Th,cosPh,sinSqTh,cos2Ph)*3./4/TMath::Pi();
      W[1]= W1spin0(sinSqTh,cosSqTh,sin2Th,cosPh,cos2Ph)*3./4/TMath::Pi();
      W[2]= W2spin0(sin2Th, sinPh,sinSqTh, sin2Ph)*3./4/TMath::Pi();
      W[3]= W3spin0(sin2Th,sinPh,sinSqTh,sin2Ph)*3./4/TMath::Pi();
 
    }
    //+ other elctroproduced see eqn(83-85) Schilling and Wolf

    //equation (82) Schilling and Wolf, cf eqn (29) Schilling Seyboth Wolf
    double result=0.;
    _photonPol->Calc(); //epsilon,delta and phi should all now be set

    for(uint alpha=0; alpha < 8; alpha++ ){
      result += ( W[alpha] * (*_photonPol)[alpha] ) ;
    }

    //  std::cout<<(*_photonPol)[1]<<" "<<TMath::Cos(2*_photonPol->Phi())<<std::endl;
    // return 0.5+0.5*(*_photonPol)[1];
    
    //result/=2*W[0]; //Divide by max value to get weight
    result/=1.1; //max seems to be slightly>1 should check this
    //std::cout<<"      "<< W[0]<<" "<<(*_photonPol)[0]<<" "<<W[1]<<" "<<(*_photonPol)[1]<<" phi "<<TMath::RadToDeg()*_photonPol->Phi()<<" "<<0.5 + 0.5*TMath::Cos(2*_photonPol->Phi())<<std::endl;
    // std::cout<<"    Sigma  0 : "<<_rho->Re(0,1,1)<<" "<<_rho->Re(0,0,0)<<" "<<_rho->Re(0,1,-1)<<std::endl;
    // std::cout<<"    Sigma  1 :  "<<_rho->Re(1,1,1)<<" "<<_rho->Re(1,0,0)<<" "<<_rho->Re(1,1,-1)<<std::endl;
    // std::cout<<"    Sigma  2 :  "<<_rho->Im(2,1,1)<<" "<<_rho->Im(2,1,-1)<<std::endl;
    // _countRho000+=_rho->Re(0,0,0);
    // _countRho011+=_rho->Re(0,1,1);
    // _countRho01m1+=_rho->Re(0,1,-1);
    // _countRho100+=_rho->Re(1,0,0);
    // _countRho111+=_rho->Re(1,1,1);
    // _countRho11m1+=_rho->Re(1,1,-1);
    // _countRho210+=_rho->Im(2,1,0);
    // _countRho211+=_rho->Im(2,1,1);
    // _countRho21m1+=_rho->Im(2,1,-1);
    // ++_countN;
    // std::cout<<"Mean SDMEs "<<" rho000 "<<_countRho000/_countN<<" rho011 "<<_countRho011/_countN<<" rho01m1 "<<_countRho01m1/_countN<<" rho100 "<<_countRho100/_countN<<" rho111 "<<_countRho111/_countN<<" rho11m1 "<<_countRho11m1/_countN<<" rho210 "<<_countRho210/_countN<<" rho211 "<<_countRho211/_countN<<" rho21m1 "<<_countRho21m1/_countN<<(*_photonPol)[1]<<" "<<(*_photonPol)[2]<<std::endl;
    //result = W[0] + 0.5*W[1]*TMath::Cos(2*_photonPol->Phi());
    //return result = 0.5 + 0.5*TMath::Cos(2*_photonPol->Phi());
    //std::cout<<"result "<<result<<" "<<W[0]<<" "<<W[1]<<" "<<W[2]<<" "<<W[3]<<" "<<W[4]<<" "<<W[5]<<std::endl;
    if(result>1.0 || result<0){
      std::cout<<"VectorSDMEDecay invalid result "<< result<<std::endl;
      double result2=0.;
      for(uint alpha=0; alpha < 8; alpha++ ){
	std::cout<<alpha<<"      "<< W[alpha]<<" "<<(*_photonPol)[alpha]<<" "<<( W[alpha] * (*_photonPol)[alpha])<<std::endl;
	result2 += ( W[alpha] * (*_photonPol)[alpha]) ;
	std::cout<<result2<<std::endl;
      }
      std::cout<<_rho->Re(0,0,0)<<" 010 "<<_rho->Re(0,1,0)<<" 01-1 "<<_rho->Re(0,1,-1)<<" 111 "<<_rho->Re(1,1,1)<<" 100 "<<_rho->Re(1,0,0)<<" 110 "<<_rho->Re(1,1,0)<<" 11-1 "<<_rho->Re(1,1,-1)<<" 210 "<<_rho->Re(2,1,0)<<" 21-1 "<<_rho->Re(2,1,-1)<<" "<<std::endl;
      std::cout<<"epsilon " <<_photonPol->Epsilon()<<" delta " <<_photonPol->Delta()<<" phi " <<_photonPol->Phi()*TMath::RadToDeg()<<std::endl;
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
