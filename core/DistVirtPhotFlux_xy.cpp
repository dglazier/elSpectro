#include "DistVirtPhotFlux_xy.h"
#include "FunctionsForElectronScattering.h"

#include <TRandom.h>

namespace elSpectro{

  DistVirtPhotFlux_xy::DistVirtPhotFlux_xy(double eb, double mion, double Wmin):
    _ebeam(eb),
    _mTar(mion)
   {
     SetWThreshold(Wmin);
   }
  
  void DistVirtPhotFlux_xy::SetWThreshold(double Wmin){//actually want to set mass nucleon
    if(_ebeam==0) {
      std::cerr<<"DistVirtPhotFlux_xy::SetWThreshold ebeam not set"<<std::endl;
      exit(0);
    }
    //r = W^2 - M^2 + Q2 = 2M(Eg) = 2M*ebeam*y
    //r_min = Wthreshold^2 - M^2 + 0 (assume Q2min =0)  i.e.-> 0 for elastic
    //r_max = 2M(Egmax) = 2M(Ebeam-m)
    //y = r/(2*M*ebeam)
    //Q2min = m^2*y^2/(1-y)
    //Q2max(r) = r + M^2 - Wthresh^2 = W^2 + Q2 - Wthres^2
    
    _Wthresh2=Wmin*Wmin;
    double rmin=_Wthresh2 - _mTar*_mTar;
    double ymin = rmin/(2*_mTar*_ebeam);
    if(ymin==0) ymin=1E-16;
    _lnymin=TMath::Log(ymin);

    //    double ymax = _ebeam - escat::M_el();
    _lnymax=TMath::Log(1.);

    double xmin =escat::M2_el()*ymin/(1-ymin)/(2*_mTar*_ebeam);
    _max_val=escat::flux_dlnxdlny(_ebeam,TMath::Log(xmin),_lnymin);
    }

  void DistVirtPhotFlux_xy::FindWithAcceptReject(){

    double lny = gRandom->Uniform(_lnymin,_lnymax);
    double y = TMath::Exp(lny);

    //now calculate x limits
    //y = r/2ME and x = Q2/r 
    double r = 2*_mTar*_ebeam*y;
    // double Q2min = M2_el()*y*y/(1-y);
    double xmin =escat::M2_el()*y/(1-y)/(2*_mTar*_ebeam); //Q2min/r;
    
    double Q2max = r + _mTar*_mTar - _Wthresh2;  
    double xmax = Q2max/r;  //= 1 when _Wthres2=_mTar*_mTar

    double lnx = gRandom->Uniform(TMath::Log(xmin),TMath::Log(xmax)); 
    
    while(  gRandom->Uniform()*_max_val >
	    (_val=escat::flux_dlnxdlny(_ebeam,lnx,lny)) ){

      if(_val>_max_val){ std::cerr<<"DistVirtPhotFlux_xy "<<_val<<" < " <<_max_val<<std::endl; exit(0);}

      lny = gRandom->Uniform(_lnymin,_lnymax);
      y = TMath::Exp(lny);
      r = 2*_mTar*_ebeam*y;
      xmin =escat::M2_el()*y/(1-y)/(2*_mTar*_ebeam); 
      
      Q2max = r + _mTar*_mTar - _Wthresh2;  
      xmax = Q2max/r;
   
      lnx = gRandom->Uniform(TMath::Log(xmin),TMath::Log(xmax)); 
    }
   
     
    _xy=std::make_pair(TMath::Exp(lnx),TMath::Exp(lny));
   
  }
}
