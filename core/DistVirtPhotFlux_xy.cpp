#include "DistVirtPhotFlux_xy.h"

#include <TRandom.h>

namespace elSpectro{

  DistVirtPhotFlux_xy::DistVirtPhotFlux_xy(double eb, double mion, double Wmin):
    _ebeam(eb),
    _mTar(mion)
   {
     SetWThreshold(Wmin);
   }
  
  void DistVirtPhotFlux_xy::SetWThreshold(double Wmin){
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

    //Find the minimum possible y given Q2 or e-' theta limits
    double rmin=_Wthresh2 - _mTar*_mTar ;
    
    if(rmin<0) rmin=0; //can be very small -ve
    
    if( _requestQ2min>0)
      rmin=rmin + _requestQ2min;
    
    if(_requestThmin>0){
      auto rminTh = escat::ymin_Th(_ebeam,_requestThmin,_mTar)*(2*_mTar*_ebeam);
      if(rminTh>rmin)rmin=rminTh;
      _requestCosThmin=TMath::Cos(_requestThmin);
    }
    
    double ymin = rmin/(2*_mTar*_ebeam);
    if(ymin==0) ymin=1E-30;
       //now find the y limits with this threshold
    std::cout<<"ymin "<<ymin<<" "<<rmin<<" "<<_Wthresh2<<" "<<_mTar<<" "<<_ebeam<<" "<<_Wthresh2-_mTar*_mTar<<" "<<_requestQ2min<<std::endl;
 
    //if(_requestYmin!=0){
    if(_requestYmin>ymin)
      _lnymin=TMath::Log(_requestYmin);
    else
      _lnymin=TMath::Log(ymin);
    //}
    
    if(_requestYmax!=0)
      _lnymax=TMath::Log(_requestYmax);
    else
      _lnymax=TMath::Log(1.);
    
    double xmin =XMin(ymin);						 
    
    //random search for the maximum value
    _max_val=0;
    for(int i=0;i<1E5;i++){
      //note r=Q2/A2 when W=Mtar so xmin==1, so just scan from low x instead
      double val  = escat::flux_dlnxdlny(_ebeam,gRandom->Uniform(TMath::Log(1E-40),TMath::Log(1)),gRandom->Uniform(_lnymin,_lnymax));
      if(val>_max_val)_max_val=val;
    }
    _max_val*=1.01; //to be sure got max

    
    std::cout<<"DistVirtPhotFlux_xy max "<<_max_val<<" within y limits "<<TMath::Exp(_lnymin)<<" "<<TMath::Exp(_lnymax)<<std::endl;
    std::cout<<"  Other limits : "<<std::endl;
    
    if(_requestQ2min!=0)std::cout<<"\tQ2min "<<_requestQ2min<<std::endl;
    if(_requestQ2max!=0)std::cout<<"\tQ2max "<<_requestQ2max<<std::endl;
    if(_requestThmin!=0)std::cout<<"\tThmin "<<_requestThmin*TMath::RadToDeg()<<std::endl;
    if(_requestThmax!=0)std::cout<<"\tThmax "<<_requestThmax*TMath::RadToDeg()<<std::endl;
    if(_requestXmin!=0)std::cout<<"\tXmin "<<_requestXmin<<std::endl;
    if(_requestXmax!=1)std::cout<<"\tXmax "<<_requestXmax<<std::endl;
    
  }
  

  
  void DistVirtPhotFlux_xy::FindWithAcceptReject(){
    
    double lny=0;
    double lnx=0;
    
    auto getRandomXY = [&lny,&lnx,this](){
      
      lny = gRandom->Uniform(_lnymin,_lnymax);
      double y = TMath::Exp(lny);

      //calculate the fraction of x-space available
      //now calculate x limits
      //y = r/2ME and x = Q2/r 
      double avail_xmax = XMax(y);
      double avail_xmin = XMin(y);
  
      auto xrange= avail_xmax-avail_xmin;
      if(xrange<1) //correct integral over x for x-range that is not 0-1, i.e.e to get the correct y dependendence
	if( gRandom->Uniform()> (xrange) ){ lnx=1;return; }
      
    
      lnx = gRandom->Uniform(TMath::Log(avail_xmin),TMath::Log(avail_xmax));
      
    };

    getRandomXY();//intial sample

   
    while(  gRandom->Uniform()*_max_val >
	    (_val=escat::flux_dlnxdlny(_ebeam,lnx,lny)) ) {
      if(_val>_max_val){
	_max_val=_val;
	std::cout<<" MAX REACHED "<<_val<<" "<<_max_val<<std::endl;
	exit(0);
      }
    
      getRandomXY(); //sample anothe pair
    
    }
    
  
  _xy=std::make_pair(TMath::Exp(lnx),TMath::Exp(lny));
  
  }
 
}//namespace
