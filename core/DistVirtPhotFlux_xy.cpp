#include "DistVirtPhotFlux_xy.h"
#include "FunctionsForElectronScattering.h"

#include <TRandom.h>

namespace elSpectro{

  DistVirtPhotFlux_xy::DistVirtPhotFlux_xy(float ebeam,float xmin,float xmax,float ymin,float ymax):
    _ebeam(ebeam),
    _lnxmin(TMath::Log(xmin)),
    _lnxmax(TMath::Log(xmax)),
    _lnymin(TMath::Log(ymin)),
    _lnymax(TMath::Log(ymax))
  {
    _max_val=escat::flux_dlnxdlny(_ebeam,_lnxmin,_lnymin);
  }
  
  void DistVirtPhotFlux_xy::FindWithAcceptReject(){
    //   std::cout<<"DistVirtPhotFlux_xy::FindWithAcceptReject() "<<_lnxmin<<" "<<_lnxmax<<" "<<_lnymin<<" "<<_lnymax<<" "<<std::endl;
    //std::cout<<"DistVirtPhotFlux_xy::FindWithAcceptReject() "<<TMath::Exp(_lnxmin)<<" "<<TMath::Exp(_lnxmax)<<" "<<TMath::Exp(_lnymin)<<" "<<TMath::Exp(_lnymax)<<" "<<std::endl;
    double lnx = gRandom->Uniform(_lnxmin,_lnxmax); 
    double lny = gRandom->Uniform(_lnymin,_lnymax);
    
   
    while(  gRandom->Uniform()*_max_val >
	    (_val=escat::flux_dlnxdlny(_ebeam,lnx,lny)) ){
      if(_val>_max_val){ std::cerr<<"DistVirtPhotFlux_xy "<<_val<<" < " <<_max_val<<std::endl; exit(0);}
      lnx = gRandom->Uniform(_lnxmin,_lnxmax); 
      lny = gRandom->Uniform(_lnymin,_lnymax);
    }
    //delog
    // double x=TMath::Exp(lnx);
    //double y=TMath::Exp(lny);

    //    std::cout<<"   result "<<lnx<<" "<<lny<<" "<<TMath::Exp(lnx)<<" "<<TMath::Exp(lny)<<std::endl;
    
    _xy=std::make_pair(TMath::Exp(lnx),TMath::Exp(lny));
   
  }
}
