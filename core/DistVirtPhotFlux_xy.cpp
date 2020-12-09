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
    if(ymin==0) ymin=1E-50;
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

    _maxPossiblexRange=1;
    for(int i=0;i<1E5;i++){
      double avail_xmin = XMin(static_cast<double>(i)/1E5);
      if(avail_xmin<_maxPossiblexRange)
	_maxPossiblexRange=avail_xmin;
    }
    
    std::cout<<"DistVirtPhotFlux_xy max "<<_max_val<<" within y limits "<<TMath::Exp(_lnymin)<<" "<<TMath::Exp(_lnymax)<<" minimum possible x "<<_maxPossiblexRange<<std::endl;
    std::cout<<"  Other limits : "<<std::endl;

    // if(_maxPossiblexRange==0)
    // _maxPossiblexRange=200;
    //else
    //  _maxPossiblexRange=-TMath::Log(_maxPossiblexRange);
    
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
  
      if(avail_xmin>=1){ lnx=1;return; }

      if(_maxPossiblexRange)
	//for efficiency we need not sample below lowest possible x value
	lnx = gRandom->Uniform(TMath::Log(_maxPossiblexRange),TMath::Log(1));
      else 
 	lnx = gRandom->Uniform(TMath::Log(1E-50),TMath::Log(1));
  
      //check if we are within allowed x-range
      //if not return and throw another y value
      auto currx=TMath::Exp(lnx);
      if(currx>avail_xmax){ lnx=1;return; }
      if(currx<avail_xmin){ lnx=1;return; }
    };

    getRandomXY();//intial sample

   
    while(  gRandom->Uniform()*_max_val >
	    (_val=escat::flux_dlnxdlny(_ebeam,lnx,lny)) ) {
      if(_val>_max_val){
	_max_val=_val;
	std::cout<<" MAX REACHED "<<_val<<" "<<_max_val<<std::endl;
	exit(0);
      }
      getRandomXY(); //sample another pair
    }
    

    //now we want the value of ths function to be Photon flux as function of x and y
    auto x=TMath::Exp(lnx);
    auto y=TMath::Exp(lny);
    
    _val=escat::flux_dxdy(_ebeam,x,y);
     
    //return x and y values
    _xy=std::make_pair(x,y);
  
  }
   void DistVirtPhotFlux_xy::FindFlat(){

     _val=0;
     _xy=std::make_pair(-1,-1); //unphysical should never be used
    
     double y = gRandom->Uniform( TMath::Exp(_lnymin),TMath::Exp(_lnymax) );
     
     //calculate the fraction of x-space available
     //now calculate x limits
     //y = r/2ME and x = Q2/r 
     double avail_xmax = XMax(y);
     double avail_xmin = XMin(y);
  
     if(avail_xmin>=1){ return; }

     
     double x=-1;
     if(_maxPossiblexRange)
       //for efficiency we need not sample below lowest possible x value
       x = gRandom->Uniform(_maxPossiblexRange,1);
     else 
       x = gRandom->Uniform(0,1);
     
     //check if we are within allowed x-range
     if( x<avail_xmax && x>avail_xmin )
       _val=escat::flux_dxdy(_ebeam,x,y);
     
     //return x and y values
     _xy=std::make_pair(x,y);
     
     
      
     return;
     
   }
  /*
   void DistVirtPhotFlux_xy::FindFlat(){

     TF1 sampleX("sampleX","TMath::Exp(-x[0])",0,1);
     sampleX.SetRange(0,1);
     sampleX.SetNpx(1E4);
     TF1 sampleY("sampleY","TMath::Exp(-x[0])",0,1);
     sampleY.SetRange(TMath::Exp(_lnymin),TMath::Exp(_lnymax));
     sampleY.SetNpx(1E4);
    
     _val=0;
     _xy=std::make_pair(-1,-1); //unphysical should never be used
    
      double y = sampleY.GetRandom();
      double sampleYVal=sampleY.Eval(y)/sampleYIntegral;
     
     //calculate the fraction of x-space available
     //now calculate x limits
     //y = r/2ME and x = Q2/r 
     double avail_xmax = XMax(y);
     double avail_xmin = XMin(y);
  
     if(avail_xmin>=1){ return; }

     double x=-1;
     while(x<_maxPossiblexRange)
        x = sampleX.GetRandom();
     
     double sampleXVal=sampleX.Eval(x)/sampleXIntegral;
    
     //check if we are within allowed x-range
     if( x<avail_xmax && x>avail_xmin )
       _val=escat::flux_dxdy(_ebeam,x,y)/sampleYVal/sampleXVal;
     
     //return x and y values
     _xy=std::make_pair(x,y);
     
     
      
     return;
     
   }

  */
}//namespace
