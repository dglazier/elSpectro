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
    //Wmin=1;
    _Wthresh2=Wmin*Wmin;
    double rmin=_Wthresh2 - _mTar*_mTar;
    double ymin = rmin/(2*_mTar*_ebeam);
    if(ymin==0) ymin=1E-30;
    _lnymin=TMath::Log(ymin);

    //    double ymax = _ebeam - escat::M_el();
    _lnymax=TMath::Log(1.);

    double xmin =escat::M2_el()*ymin/(1-ymin)/(2*_mTar*_ebeam);


    //random search for the maximum value
    _max_val=0;
    for(int i=0;i<1E4;i++){
      double val  = escat::flux_dlnxdlny(_ebeam,gRandom->Uniform(TMath::Log(xmin),TMath::Log(1)),gRandom->Uniform(_lnymin,_lnymax));
    if(val>_max_val)_max_val=val;
    }
    std::cout<<"DistVirtPhotFlux_xy max "<<_max_val<<" x "<<xmin<<" y "<<ymin<<std::endl;
  }

  void DistVirtPhotFlux_xy::FindWithAcceptReject(){

    double lny=0;
    double lnx=0;
    
    auto getRandomXY = [&lny,&lnx,this](){
      
      lny = gRandom->Uniform(_lnymin,_lnymax);
      double y = TMath::Exp(lny);

      //now calculate x limits
      //y = r/2ME and x = Q2/r 
      double r = 2*_mTar*_ebeam*y;
    
      //double xmin =escat::M2_el()*y/(2*_mTar*_ebeam)/(1-y); //Q2min/r;
      double xmin = 1E-25;
        
      // double xmax = Q2max/r;  //= 1 when _Wthres2=_mTar*_mTar
      //Note the calculation below is more numerically stable with very small numbers
      double xmax = 1 + (_mTar*_mTar - _Wthresh2)/r;
      
      lnx = gRandom->Uniform(TMath::Log(xmin),TMath::Log(xmax));
      //lnx = gRandom->Uniform(TMath::Log(xmin),TMath::Log(1));
      auto x= TMath::Exp(lnx);
      // if(x>xmax){lnx=1;return;} //unphysical will resample
      auto Q2=escat::Q2_xy(_ebeam,x,y );

      double Q2max = r + _mTar*_mTar - _Wthresh2;  
      Q2max = Q2max>_requestQ2max ?_requestQ2max : Q2max ; //override with given max if exists
      double Q2min = escat::M2_el()*y*y/(1-y);
      Q2min = Q2min<_requestQ2min ?_requestQ2min : Q2min ; //override with given min if exists
 
      if(Q2 < Q2min ) lnx=1; //unphysical will resample
      if(Q2 > Q2max ) lnx=1; //unphysical will resample
      // if(xmax!=1) std::cout<<"x "<<TMath::Exp(lnx)<<" "<<xmin<<" \t"<<xmax<<" "<<r<<" y \t"<<y<<" "<<"\t"<<escat::flux_dlnxdlny(_ebeam,lnx,lny)  <<"\t"<<escat::flux_dlnxdlny(_ebeam,lnx,lny)/_max_val  <<" max is " <<_max_val<<std::endl;
    };

    getRandomXY();
    
    while(  gRandom->Uniform()*_max_val >
	    (_val=escat::flux_dlnxdlny(_ebeam,lnx,lny)) ) {
      //|| escat::Q2_xy( _ebeam,TMath::Exp(lnx),TMath::Exp(lny) )< _requestQ2min
      //if(_val>_max_val){ std::cerr<<"DistVirtPhotFlux_xy "<<_val<<" < " <<_max_val<<std::endl; exit(0);}
      if(_val>_max_val){ _max_val=_val;}
      
      getRandomXY(); //sample anothe pair

      /*
      lny = gRandom->Uniform(_lnymin,_lnymax);
      y = TMath::Exp(lny);
      r = 2*_mTar*_ebeam*y;
      xmin =escat::M2_el()*y/(1-y)/(2*_mTar*_ebeam); 
      
      Q2max = r + _mTar*_mTar - _Wthresh2;  
      xmax = Q2max/r;
   
      lnx = gRandom->Uniform(TMath::Log(xmin),TMath::Log(xmax));
      */
    }
    
     
    _xy=std::make_pair(TMath::Exp(lnx),TMath::Exp(lny));
   
  }
}
  /*

   double lny = gRandom->Uniform(_lnymin,_lnymax);
    double y = TMath::Exp(lny);

    //now calculate x limits
    //y = r/2ME and x = Q2/r 
    double r = 2*_mTar*_ebeam*y;
    double Q2min = M2_el()*y*y/(1-y);
    Q2min = Q2min<_requestQ2min ?_requestQ2min : Q2min ; //override with given min if exists
    
    double xmin =escat::M2_el()*y/(1-y)/(2*_mTar*_ebeam); //Q2min/r;
    
    double Q2max = r + _mTar*_mTar - _Wthresh2;  
    Q2max = Q2max>_requestQ2max ?_requestQ2max : Q2max ; //override with given max if exists
    double xmax = Q2max/r;  //= 1 when _Wthres2=_mTar*_mTar

    double lnx = gRandom->Uniform(TMath::Log(xmin),TMath::Log(xmax)); 
   

   */
