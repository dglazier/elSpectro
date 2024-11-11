#include "DistTH2Slice.h"
#include <TRandom.h>
#include <TMath.h>
#include <iostream>

namespace elSpectro{

  DistTH2Slice::DistTH2Slice(const TH2D& ff):
    DistTH2{ff}{
    auto& th2 = GetTH2();
 
    //cache integrals
    _binIntegrals.resize(th2.GetNbinsX(), std::vector<double >(th2.GetNbinsY() + 1) );
    _binMax.resize(th2.GetNbinsX());
    _binTotal.resize(th2.GetNbinsX());
    
    _binIntegrals[0][0]=0;
    _maxBinVal=0;
    for(auto ix=1;ix<=th2.GetNbinsX();++ix){//loop over x
      auto vx = ix-1;
      //_binIntegrals[vx][0]=0;
      _binMax[vx]=0;
      for(auto iy=1;iy<=th2.GetNbinsY();++iy){//loop over y

	double cont = th2.GetBinContent(ix,iy)*th2.GetYaxis()->GetBinWidth(iy);
	//running sum for binary search
	//	_binIntegrals[vx][iy] = _binIntegrals[vx][iy-1] + cont;
	_binIntegrals[vx][iy] = _binIntegrals[vx][iy-1] + cont;

	if(th2.GetBinContent(ix,iy)>_binMax[vx])_binMax[vx]=th2.GetBinContent(ix,iy);//save highest bin value for weight
	//	if(cont>_maxBinVal)_maxBinVal=cont;
	 
      }//x
    }//y
    _maxBinVal=th2.GetMaximum();
    _maxIntegral=0;
    //normalise total integral in each x bin = 1
    for(auto ix=1;ix<=th2.GetNbinsX();++ix){//loop over x
      auto vx = ix-1;
      auto total = _binIntegrals[vx][th2.GetNbinsY()];
      _binTotal[vx]= total;
      if(_maxIntegral<total)_maxIntegral = total;
      // std::cout<<"DistTH2Slice::DistTH2Slice total "<<total<<std::endl;
      for(auto iy=1;iy<=th2.GetNbinsY();++iy){//loop over y
	
	if(total==0||TMath::Abs(total)==TMath::Infinity()){ //if zero total, just make equidistant integrals so samples uniformly
	  _binIntegrals[vx][iy] = _binIntegrals[vx][iy-1] +  th2.GetYaxis()->GetBinWidth(iy)/(th2.GetYaxis()->GetXmax()-th2.GetYaxis()->GetXmin());
	  
	}
	else{
	  //running sum for binary search
	  _binIntegrals[vx][iy] = _binIntegrals[vx][iy]/total;
	}
	//	std::cout<<"DistTH2Slice::DistTH2Slice "<< _binIntegrals[vx][iy]<<std::endl;
	if(TMath::IsNaN(_binIntegrals[vx][iy])) _binIntegrals[vx][iy] =0;
 	//std::cout<<"DistTH2Slice::DistTH2Slice "<< vx<<" "<<iy<<" "<<_binIntegrals[vx][iy]<<std::endl;
     }//x
    }//y

    //distribution with highest x values
    //used to get max for given x,
    //provides interpolation rather than binned value
    TH1D hxhigh("xhigh","xhigh",th2.GetNbinsX(),th2.GetXaxis()->GetXbins()->GetArray());
    for(UInt_t i=1;i<=th2.GetNbinsX();++i){
      hxhigh.SetBinContent(i,_binMax[i-1]);
    }
    _distHighX  =  DistTH1{hxhigh};

   }
  
  void DistTH2Slice::RandomXY() noexcept
  {
    //fix slice in x, x value must be fixed aready
    auto x = GetX();
    auto& th2 = GetTH2();
    auto binx = th2.GetXaxis()->FindBin(x) -1 ; //we start from 0 not 1 like TH1

     if(binx>=th2.GetNbinsX())binx=th2.GetNbinsX()-1;
    //first check on integrated x
    // if(gRandom->Uniform()>_binTotal[binx]/_maxIntegral){
    //   SetVal(0);
    //   return;
    // }
       
    //get random number to find how far along y
    Double_t r1 = gRandom->Rndm(); //between 0 and 1
    //need to use the slice in x y-bins
    //_binIntegrals[binx] = cumulative y vector at x
    Int_t ibin = TMath::BinarySearch(_binIntegrals[binx].size(),_binIntegrals[binx].data(),r1);

    // if(ibin == (_binIntegrals[binx].size()) ){
    //   //all bins empty, just return random flat
    //   SetY(gRandom->Uniform( GetMinY(), GetMaxY()));
    // }
    // std::cout<<"DistTH2Slice "<<ibin<<" "<<binx<<" "<<x<<std::endl;
    Double_t y = th2.GetYaxis()->GetBinLowEdge(ibin+1);
    //Double_t y = th2.GetYaxis()->GetBinLowEdge(ibin);
    // std::cout<<"y from low edge "<<y<<" "<<th2.GetYaxis()->GetBinLowEdge(ibin+2)<<" "<<_binIntegrals[binx][ibin+1]<<" "<<_binIntegrals[binx][ibin]<<" "<<_binIntegrals[binx][ibin-1]<<std::endl;
    if (r1 > _binIntegrals[binx][ibin]){
      if(_binIntegrals[binx][ibin+1] - _binIntegrals[binx][ibin] != 0 ){
	y += th2.GetYaxis()->GetBinWidth(ibin+1)*(r1-_binIntegrals[binx][ibin])/
	  (_binIntegrals[binx][ibin+1] - _binIntegrals[binx][ibin]);
      }
    }
     SetY(y);
    
    //_current_max = _binMax[binx];
    _current_max =_distHighX.GetValueFor(x);;
    if(_current_max==0){ //zero bin entries, just returning a uniform random 
      //make GetCurrentWeight return 1
      _current_max = 1;
      SetVal(1);
    }
    else{
      //auto val1 = LinearInterpolateAlongX(x,ibin+1,binx+1,binx+2);
      //auto val2 = LinearInterpolateAlongX(x,ibin+1,binx,binx+1);
      //auto val = val1>val2 ? val1:val2;
     SetVal(th2.GetBinContent(th2.FindFixBin(GetX(),GetY())));//ibin starts from 0
     // SetVal(val);//ibin starts from 0
      _current_max = 1;
    }

    // std::cout<<" DistTH2Slice::SampleSingle() x "<<x<<" y "<<y<<" "<<CurrentValue()<<" "<<MaxValue()<<" "<<GetCurrentWeight()<<" next w "<<th2.GetBinContent(th2.FindFixBin(GetX(),GetY())+1)<<" prev w "<<th2.GetBinContent(th2.FindFixBin(GetX(),GetY())-1)<<" interp "<<val1<<" "<<val2<<std::endl;
    // std::cout<<"DistTH2Slice done "<<std::endl;
 
  }

 
  // DistTH2Slice::DistTH2Slice(const TH2D& ff):
  //   DistTH2{ff}{
  //   auto& th2 = GetTH2();
 
  //   //cache integrals
  //   _binIntegrals.resize(th2.GetNbinsX(), std::vector<double >(th2.GetNbinsY() + 1) );
  //   _binMax.resize(th2.GetNbinsX());
  //   _binTotal.resize(th2.GetNbinsX());
    
  //   _binIntegrals[0][0]=0;
  //   _maxBinVal=0;
  //   for(auto ix=1;ix<=th2.GetNbinsX();++ix){//loop over x
  //     auto vx = ix-1;
  //     _binIntegrals[vx][0]=0;
  //     _binMax[vx]=0;
  //     for(auto iy=1;iy<=th2.GetNbinsY();++iy){//loop over y

  // 	double cont = th2.GetBinContent(ix,iy)*th2.GetYaxis()->GetBinWidth(iy);
  // 	//running sum for binary search
  // 	//	_binIntegrals[vx][iy] = _binIntegrals[vx][iy-1] + cont;
  // 	_binIntegrals[vx][iy] = _binIntegrals[vx][iy-1] + cont;

  // 	if(cont>_binMax[vx])_binMax[vx]=cont;//save highest bin value for weight
  // 	//	if(cont>_maxBinVal)_maxBinVal=cont;
	 
  //     }//x
  //   }//y
  //   _maxBinVal=th2.GetMaximum();
  //   _maxIntegral=0;
  //   //normalise total integral in each x bin = 1
  //   for(auto ix=1;ix<=th2.GetNbinsX();++ix){//loop over x
  //     auto vx = ix-1;
  //     auto total = _binIntegrals[vx][th2.GetNbinsY()];
  //     _binTotal[vx]= total;
  //     if(_maxIntegral<total)_maxIntegral = total;
  //     // std::cout<<"DistTH2Slice::DistTH2Slice total "<<total<<std::endl;
  //     for(auto iy=1;iy<=th2.GetNbinsY();++iy){//loop over y
	
  // 	if(total==0||TMath::Abs(total)==TMath::Infinity()){ //if zero total, just make equidistant integrals so samples uniformly
  // 	  _binIntegrals[vx][iy] = _binIntegrals[vx][iy-1] +  th2.GetYaxis()->GetBinWidth(iy)/(th2.GetYaxis()->GetXmax()-th2.GetYaxis()->GetXmin());
	  
  // 	}
  // 	else{
  // 	  //running sum for binary search
  // 	  _binIntegrals[vx][iy] = _binIntegrals[vx][iy]/total;
  // 	}
  // 	//	std::cout<<"DistTH2Slice::DistTH2Slice "<< _binIntegrals[vx][iy]<<std::endl;
  // 	if(TMath::IsNaN(_binIntegrals[vx][iy])) _binIntegrals[vx][iy] =0;
  // 	//std::cout<<"DistTH2Slice::DistTH2Slice "<< vx<<" "<<iy<<" "<<_binIntegrals[vx][iy]<<std::endl;
  //    }//x
  //   }//y
  //  }
  
  // void DistTH2Slice::RandomXY() noexcept
  // {
  //   //fix slice in x, x value must be fixed aready
  //   auto x = GetX();
  //   auto& th2 = GetTH2();
  //   auto binx = th2.GetXaxis()->FindBin(x) -1 ; //we start from 0 not 1 like TH1

  //   //first check on integrated x
  //   if(gRandom->Uniform()>_binTotal[binx]/_maxIntegral){
  //     SetVal(0);
  //     return;
  //   }
       
  //   //get random number to find how far along y
  //   Double_t r1 = gRandom->Rndm(); //between 0 and 1
  //   //need to use the slice in x y-bins
  //   //_binIntegrals[binx] = cumulative y vector at x
  //   Int_t ibin = TMath::BinarySearch(_binIntegrals[binx].size(),_binIntegrals[binx].data(),r1);

  //   // if(ibin == (_binIntegrals[binx].size()) ){
  //   //   //all bins empty, just return random flat
  //   //   SetY(gRandom->Uniform( GetMinY(), GetMaxY()));
  //   // }
  //   // std::cout<<"DistTH2Slice "<<ibin<<" "<<binx<<" "<<x<<std::endl;
  //   Double_t y = th2.GetYaxis()->GetBinLowEdge(ibin+1);
  //   //Double_t y = th2.GetYaxis()->GetBinLowEdge(ibin);
  //   // std::cout<<"y from low edge "<<y<<" "<<th2.GetYaxis()->GetBinLowEdge(ibin+2)<<" "<<_binIntegrals[binx][ibin+1]<<" "<<_binIntegrals[binx][ibin]<<" "<<_binIntegrals[binx][ibin-1]<<std::endl;
  //   if (r1 > _binIntegrals[binx][ibin]){
  //     if(_binIntegrals[binx][ibin+1] - _binIntegrals[binx][ibin] != 0 ){
  // 	y += th2.GetYaxis()->GetBinWidth(ibin+1)*(r1-_binIntegrals[binx][ibin])/
  // 	  (_binIntegrals[binx][ibin+1] - _binIntegrals[binx][ibin]);
  //     }
  //   }
  //   //std::cout<<"y after correct "<<y<<_binIntegrals[binx][ibin+1] - _binIntegrals[binx][ibin]<<" "<<_binIntegrals[binx][ibin+1]<<" "<<_binIntegrals[binx][ibin]<<" "<<th2.GetYaxis()->GetBinWidth(ibin+1)<<std::endl;
  //   // std::cout<<"DistTH2Slice y "<<y<<std::endl;
 
  //   SetY(y);
    
  //   //_current_max = _binMax[binx];
  //   //_current_max = _maxIntegral;//*(th2.GetYaxis()->GetXmax()-th2.GetYaxis()->GetXmin()) ;//*(_binMax[binx])/_maxBinVal;//need to normalise to integral of slice at x
  //   _current_max = _maxBinVal;//need to normalise to integral of slice at x
  //   if(_current_max==0){ //zero bin entries, just returning a uniform random 
  //     //make GetCurrentWeight return 1
  //     _current_max = 1;
  //     SetVal(1);
  //   }
  //   else{
  //    SetVal(th2.GetBinContent(th2.FindFixBin(GetX(),GetY())));//ibin starts from 0
  //    // SetVal(th2.GetBinContent(th2.FindFixBin(GetX(),GetY())) / th2.GetYaxis()->GetBinWidth(th2.GetYaxis()->FindFixBin(GetY())) );//ibin starts from 0
  //    //std::cout<<"DistTH2Slice::SampleSingle() ybin "<<th2.GetYaxis()->FindFixBin(GetY())<<" "<<ibin+1<<" value "<<CurrentValue()<<" next "<<th2.GetBinContent(binx,ibin+2)<<" xbin "<<binx<<" "<<GetX()<<std::endl;
  //  }

  //   /*
  //   SetY(gRandom->Uniform(th2.GetYaxis()->GetXmin(),th2.GetYaxis()->GetXmax()));
  //   SetVal(th2.GetBinContent(th2.FindFixBin(GetX(),GetY())));
  //   if(gRandom->Uniform()>GetCurrentWeight())
  //     SetVal(0);//reject
  //   */
    
  //   //  _weight = CurrentValue()/_binTotal[binx];
  //   //std::cout<<" DistTH2Slice::SampleSingle() x "<<x<<" y "<<y<<" r1 "<<r1<<"x bin "<<binx<<" ybin "<<ibin<<"bin integrals "<<_binIntegrals[binx][ibin]<<" "<<_binIntegrals[binx][ibin+1]<<std::endl;
  //   // std::cout<<"DistTH2Slice done "<<std::endl;
 
  // }



}
