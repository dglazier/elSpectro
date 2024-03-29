#include "DistVirtPhotFlux_xy.h"

#include <TRandom.h>

//For pdf integration
#include <Math/Functor.h>
#include <RooFunctorBinding.h>
#include <RooRealVar.h>
#include <RooArgList.h>
#include <utility>
#include <functional>

namespace elSpectro{

  DistVirtPhotFlux_xy::DistVirtPhotFlux_xy(double eb, double mion, double Wmin):
    _ebeam(eb),
    _mTar(mion)
   {

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
    std::cout<<Wmin<<" "<<_mTar<<"  ymin "<<ymin<<" "<<rmin<<" "<<_Wthresh2<<" "<<_mTar<<" "<<_ebeam<<" "<<_Wthresh2-_mTar*_mTar<<" "<<_requestQ2min<<std::endl;
 
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
    


     
    //Seacrh for lowest possible x...
    _maxPossiblexRange=1;
    ymin = TMath::Exp(_lnymin);
    double ymax = TMath::Exp(_lnymax);
    double xmin =XMin(ymin);						 


    for(int i=0;i<1E5;i++){
      //     double avail_xmin = XMin(static_cast<double>(i*(ymax-ymin) + ymin)/1E5);
      //Need xmin with no experiment based limits on Q2 or theta
      double yformin=static_cast<double>(i*(ymax-ymin) + ymin)/1E5;
      double avail_xmin =escat::M2_el()*yformin/(2*_mTar*_ebeam)/(1-yformin);
       if(avail_xmin<_maxPossiblexRange)
	_maxPossiblexRange=avail_xmin;
    }

    FindMaxVal(); //For sampling

    std::cout<<"DistVirtPhotFlux_xy::SetWThreshold new Wmin "<<GetWMin() <<" "<<TMath::Sqrt(2*_mTar*_ebeam*ymin+_mTar*_mTar)<<" "<<TMath::Sqrt(rmin + _mTar*_mTar)-_mTar*_mTar*ymin*ymin/(1-ymin)<<" "<<sqrt( _mTar*(_mTar + 2*ymin*_ebeam ) -  escat::Q2_xy( _ebeam,_maxPossiblexRange,ymin))<<" xmin "<<xmin<<" "<<XMin(ymin)<<" ymin "<<ymin<<std::endl;
    _Wthresh2 =_mTar*(_mTar + 2*ymin*_ebeam )-  escat::Q2_xy( _ebeam,_maxPossiblexRange,ymin);
  
    std::cout<<"DistVirtPhotFlux_xy max "<<_max_val<<" within y limits "<<TMath::Exp(_lnymin)<<" "<<TMath::Exp(_lnymax)<<" minimum possible x "<<_maxPossiblexRange<<std::endl;
    std::cout<<"  Other limits : "<<std::endl;
    
    if(_requestQ2min!=0)std::cout<<"\tQ2min "<<_requestQ2min<<std::endl;
    if(_requestQ2max!=0)std::cout<<"\tQ2max "<<_requestQ2max<<std::endl;
    if(_requestThmin!=0)std::cout<<"\tThmin "<<_requestThmin*TMath::RadToDeg()<<std::endl;
    if(_requestThmax!=0)std::cout<<"\tThmax "<<_requestThmax*TMath::RadToDeg()<<std::endl;
    if(_requestXmin!=0)std::cout<<"\tXmin "<<_requestXmin<<std::endl;
    if(_requestXmax!=1)std::cout<<"\tXmax "<<_requestXmax<<std::endl;
    if(_requestYmin!=0)std::cout<<"\tYmin "<<_requestYmin<<std::endl;
    if(_requestYmax!=1)std::cout<<"\tYmax "<<_requestYmax<<std::endl;


   
    _lnxmin=TMath::Log(_maxPossiblexRange);
    _lnxmax=0;

    //Finally, Integrate over photon flux
    auto xvar = RooRealVar("x","x",-(_lnxmax-_lnxmin)/2,_lnxmin,_lnxmax,"");
    auto yvar = RooRealVar("y","y",-(_lnymax-_lnymin)/2,_lnymin,_lnymax,"");
    xvar.Print();
    yvar.Print();
    
    auto flambda = [this](const double *x)
      {
	if(x[0]==0) return 0.;
	if(x[1]==0) return 0.;
	auto val = Eval(x);
	return val;
      };

    
    auto wrapPdf=ROOT::Math::Functor( flambda , 2);
    auto pdf = RooFunctorPdfBinding("PdfDistVirtPhotFlux_xy", "PdfDistVirtPhotFlux_xy", wrapPdf, RooArgList(xvar,yvar));
    auto roovars= RooArgSet(xvar,yvar);
 
    _integral=pdf.getNorm(roovars);

    std::cout<<"DistVirtPhotFlux_xy INTEGRAL "<<pdf.getNorm(roovars)<<" at proton rest frame e- energy "<<_ebeam<<" and W threshold "<<TMath::Sqrt(_Wthresh2)<< std::endl;
  
  }

  void DistVirtPhotFlux_xy::FindMaxVal(){
    //random search for the maximum value
    _max_val=0;
    
    for(int i=0;i<1E7;i++){
      //note r=Q2/A2 when W=Mtar so xmin==1, so just scan from low x instead
      double ranX=gRandom->Uniform(TMath::Log(_maxPossiblexRange),TMath::Log(1));
      double ranY=gRandom->Uniform(_lnymin,_lnymax);
      double val  = escat::flux_dlnxdlny(_ebeam,ranX,ranY);
      //double val  = escat::flux_dlnxdlny(_ebeam,ranX,ranY)*WeightForW(TMath::Exp(ranX),TMath::Exp(ranY));
      if(val>_max_val){_max_val=val;}
    }
    _max_val*=1.02; //to be sure got max
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
    //  std::cout<<"DistVirtPhotFlux_xy::max "<<_max_val<<std::endl;
    while(  gRandom->Uniform()*_max_val >
	    (_val=escat::flux_dlnxdlny(_ebeam,lnx,lny)) ) {
      if(_val>_max_val){
	_max_val=_val;
	std::cout<<"DistVirtPhotFlux_xy::FindWithAcceptReject() MAX REACHED "<<_val<<" "<<_max_val<<std::endl;
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

 /*
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

    auto x=TMath::Exp(lnx);
    auto y=TMath::Exp(lny);
    while(  gRandom->Uniform()*_max_val >
    	    (_val=escat::flux_dlnxdlny(_ebeam,lnx,lny)*WeightForW(x,y) )) {
      //while(  gRandom->Uniform()*_max_val >
      //    (_val=escat::flux_dlnxdlny(_ebeam,lnx,lny)) ) {
       //if(gRandom->Uniform()>WeightForW(x,y))
       // continue;
       
      if(_val>_max_val){
    std::cout<< "check "<<escat::flux_dlnxdlny(_ebeam,TMath::Log(6.70592e-12), TMath::Log( 0.00908638))*WeightForW(6.70592e-12,0.00908638)<<std::endl;
	std::cout<<"DistVirtPhotFlux_xy::FindWithAcceptReject() MAX REACHED "<<_val<<" "<<_max_val<<" x "<<x <<" y "<<y<<" eb "<<_ebeam<<std::endl;
	exit(0);
      }
      getRandomXY(); //sample another pair
      x=TMath::Exp(lnx);
      y=TMath::Exp(lny);
    }
    

    //now we want the value of ths function to be Photon flux as function of x and y
    
    //  _val=escat::flux_dxdy(_ebeam,x,y)*WeightForW(x,y);
 
    // std::cout<<"Dist weight    "<<WeightForW(x,y)<<std::endl;
    //return x and y values
    _xy=std::make_pair(x,y);
  
  }
  */  
}//namespace
