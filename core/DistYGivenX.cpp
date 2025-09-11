#include "DistYGivenX.h"
#include "VectorUtils.h"
#include <TRandom.h>
#include <TMath.h>
#include <iostream>

namespace elSpectro{

  DistYGivenX::DistYGivenX(const bins& x,const std::vector<bins>& y,const std::vector<bins>& vals){
    ///x = x-values, first and last elements are bounds, inner values are taken as midpoints
    ///y = y[x]-values, i.e. values of y @ x. similar first and last are bounds
    ///vals = value at x,y = vals[x][y]

    //calculate low edges. We assume x contains upper and lower limit
    //as first and last entries, other entries are bin centers
    //so we now calculate the midpoints to get low edges.
    //////Helper function to calculate vector low edges
    auto makeLowEdges= [](const bins& centres){
      
      bins lowedges(centres.size()+1);
      lowedges[0]=centres[0]; //first number is lower limit
      for(uint i = 1 ; i <centres.size(); ++i){
	lowedges[i] =(centres[i]+centres[i-1])/2;//midpoint
      }
      lowedges[centres.size()]=centres[centres.size()-1]; //last number is upper limit
      return lowedges;
    };
    //////Helper function to calculate vector widths from low edges
    auto makeWidths= [](const bins& edges){
      
      bins widths(edges.size());
      for(uint i = 0 ; i < edges.size() -1; ++i){
	widths[i] =(edges[i+1]-edges[i]);//midpoint
      }
      return widths;
    };
 
    _xbins=x;
    _ybins=y;
    _vals=vals;

    
    //we are going to allow unsorted vectors and so we must sort here
    //first make sure xbins is sorted
    //    std::cout<<" print x "<<_xbins.size()<<" "<<_ybins.size()<<std::endl;
    vector_utils::synched_sort(_xbins,_ybins,_vals);
     // std::cout<<" done x "<<_xbins.size()<<" "<<_ybins.size()<<" "<<_ybins[1].size()<<std::endl;
    //vector_utils::print(_xbins);
    //now for each x bin sort and synch y and vals
    for(uint i=0;i<_ybins.size();++i){
      //    std::cout<<" print y "<< i <<" "<<_ybins[i].size()<<std::endl;
      vector_utils::synched_sort(_ybins[i],_vals[i]);
      vector_utils::zero_nans(_vals[i]);
      //vector_utils::print(_ybins[i]);
      //vector_utils::print(_vals[i]);
      // std::cout<<" done y "<< i <<_ybins[i].size()<<std::endl;
    }
    _nBinsX = _xbins.size(); //e.g. if x={-1,0,1},lowedges = -1,0.25,0.25 for 3 bins 
    
    //for each x bin there must be a vector of ybins in y
    if(_nBinsX!=_ybins.size()){
      std::cerr<<" DistYGivenX::DistYGivenX x and y size not the same "<<_nBinsX<<" "<<_ybins.size()<<std::endl;
      exit(0);
    }
     
    _xbins_lowedges = makeLowEdges(_xbins);
      
    //cache integrals
    for(const auto& ybin:_ybins){
      //create an integral for each y bin in x
      _xbins_integrals.push_back(bins(ybin.size()));
      _ybins_lowedges.push_back(makeLowEdges(ybin));
      //make widths from the just calculated low edges, and the final value.
      _ybins_widths.push_back(makeWidths( _ybins_lowedges.back() ));
    }
    
    _xbins_max=bins(_nBinsX);

    //e.g y={0,1,2,3} ylowedges = {0,0.5,1.5,2.5} vals ={a,b,c,d}
    //xbins_integrals[x]={a*w,bw+aw,cw+bw+aw,dw+cw+bw+aw}
    for(auto ix=0;ix<_nBinsX;++ix){//loop over x
      //Find the max difference between 2 ybins and set this as minimum in each ybin
      //this accounts for changing value across bin when it comes to sampling weights
      // auto max_diff=vector_utils::max_consecutive_diff(_vals[ix]);//TESTING
      vector_utils::take_max_of_neighbours(_vals[ix]);
      
      _xbins_integrals[ix][0]=0;
      _xbins_max[ix]=0;
      for(auto iy=0;iy<y[ix].size();++iy){//loop over y
	//	_vals[ix][iy]+=max_diff;//add max difference//TESTING
	_vals[ix][iy]+=_xbins_max[ix]*0.5;//add max difference//TESTING
	
	double contribution = _vals[ix][iy]*_ybins_widths[ix][iy];
	//runing sums for binary search
	_xbins_integrals[ix][iy] = iy==0 ?
	  contribution : _xbins_integrals[ix][iy-1]+contribution;

	if(_vals[ix][iy]>_xbins_max[ix])_xbins_max[ix]=_vals[ix][iy];//save highest bin value for weight
      }//y
    }//x

    //normalise total integral in each x bin = 1
    for(auto ix=0;ix<_nBinsX;++ix){//loop over x
      //the last bin_integrals is the max by definition
      auto total = _xbins_integrals[ix].back();
      for(auto iy=0;iy<y[ix].size();++iy){//loop over y
	//running sum for binary search
	_xbins_integrals[ix][iy]/=total;
	//normalise vals to max in x
	//	_vals[ix][iy]/=_xbins_max[ix];
	if(TMath::IsNaN(_xbins_integrals[ix][iy])) _xbins_integrals[ix][iy] = 0.;
      }//x
    }//y
    _maxXHist=GetMaxValVersusX();
   }
  
  void DistYGivenX::RandomXY() noexcept
  {
    //fix slice in x, x value must be fixed aready
    auto x = GetX();
    auto xbin =_current_xbin;

    // std::cout<<"DistYGivenX::RandomXY() "<<x <<" "<<xbin<<" "<<_xbins_integrals.size()<<" should be OK  "<<_xbins_lowedges.size()<<" "<<_xbins_lowedges[0]<<" "<<_xbins_lowedges[1]<< " " <<_xbins_integrals[xbin].size()<<" "<<_xbins_integrals[xbin].back()<<std::endl;

    if(_xbins_integrals[xbin].back()==0){
      SetY(gRandom->Uniform(GetMinY(),GetMaxY()));
      SetVal(1);
      _current_max = 1.;
      std::cout<<"Warning :: DistYGivenX::RandomXY() warning zero distribution at x = "<<x<<" bin "<<xbin<<std::endl;
      return;
    }
    
    //get random number to find how far along y
    Double_t r1 = gRandom->Rndm(); //between 0 and 1

    //need to use the slice in x y-bins
    //_binIntegrals[binx] = cumulative y vector at x
    // std::cout<<"DistYGivenX::RandomXY() "<<_xbins_integrals[xbin].size()<<" "<<_xbins_integrals[xbin][0]<<" "<<" "<<_xbins_integrals[xbin].back()<<" r1 "<<r1<<std::endl;
   Int_t ybin = TMath::BinarySearch(_xbins_integrals[xbin].size(),_xbins_integrals[xbin].data(),r1)+1;
   // std::cout<<"DistYGivenX::RandomXY() "<<_xbins_integrals[xbin].size()<<" "<<_xbins_integrals[xbin][0]<<" "<<_xbins_integrals[xbin][ybin]<<" "<<_xbins_integrals[xbin].back()<<" r1 "<<r1<<" ybin "<<ybin<<std::endl;
    
    //start y from low edge of bin
    Double_t y = _ybins_lowedges[xbin][ybin];
    auto low_integral = ybin>0 ? _xbins_integrals[xbin][ybin-1] :0;
    auto high_integral = _xbins_integrals[xbin][ybin];
    //and add on residual fraction of integral for this bin
    y += _ybins_widths[xbin][ybin]*(r1-low_integral)/(high_integral - low_integral);
    /*
    if(xbin>0&&xbin<_nBinsX-1){
      std::cout<<"DistYGivenX::RandomXY()x   "<< x<<" "<<xbin<<" = "<<_xbins[xbin]<<" y "<<y<<" "<<ybin<<" vals "<<_vals[xbin][ybin]<<" + "<<_vals[xbin+1][ybin]<<" - "<<_vals[xbin-1][ybin]<<" interp "<<_maxXHist.Interpolate(x)<<" "<<_vals[xbin].back()<<" "<<MaxValue()<<"   "<<_vals[xbin][ybin]/_vals[xbin].back()<<" +  "<<_vals[xbin][ybin+1]/_vals[xbin].back()<<std::endl;
      auto ybin2 = TMath::BinarySearch(_xbins_integrals[xbin+1].size(),_xbins_integrals[xbin+1].data(),r1);
      std::cout<<"DistYGivenX::RandomXY()x   "<< x<<" "<<_xbins[xbin+1]<<" y "<<y<<" "<<ybin2<<" vals "<<_vals[xbin+1][ybin2]<<" "<<_vals[xbin+1].back()<<"   "<<_vals[xbin+1][ybin2]/_vals[xbin+1].back()<<" "<< _vals[xbin+1][ybin2+1]/_vals[xbin+1].back()<<std::endl;
      ybin2 = TMath::BinarySearch(_xbins_integrals[xbin-1].size(),_xbins_integrals[xbin-1].data(),r1);
      std::cout<<"DistYGivenX::RandomXY()x   "<< x<<" "<<_xbins[xbin-1]<<" y "<<y<<" "<<ybin2<<" vals "<<_vals[xbin-1][ybin2]<<" "<<_vals[xbin-1].back()<<"  weight  "<<_vals[xbin-1][ybin2]/_vals[xbin-1].back()<<" "<< _vals[xbin-1][ybin2+1]/_vals[xbin-1].back()<<std::endl;

    }
     
     if(ybin>0&&ybin<_ybins[xbin].size()-1){
       std::cout<<"DistYGivenX::RandomXY()y   "<< x<<" "<<xbin<<" y "<<y<<" "<<ybin<<" size "<<_ybins[xbin].size()<<" vals "<<_vals[xbin][ybin]<<" + "<<_vals[xbin][ybin+1]<<" - "<<_vals[xbin][ybin-1]<<" interp "<<_maxXHist.Interpolate(x)<<" "<<_vals[xbin].back()<<" "<<MaxValue()<<" + "<<_vals[xbin][ybin+1]/_vals[xbin].back()<<std::endl;
  
     }
    */
    SetY(y);
    SetVal(_vals[xbin][ybin]);
    _current_max=_xbins_max[xbin];
    // std::cout<<"DistYGivenX::RandomXY() "<< x<<" "<<xbin<<" y "<<y<<" "<<ybin<<" vals "<<_vals[xbin][ybin]<<" vals+1 "<<_vals[xbin][ybin+1]<<" current max "<<_current_max<<" bins "<<_xbins[xbin]<<" "<<_ybins[xbin][ybin-1]<<" "<<_ybins[xbin][ybin]<<" "<<_ybins[xbin][ybin+1]<<" width "<<_ybins_widths[xbin][ybin]<<" frac "<<(r1-low_integral)/(high_integral - low_integral)<<" "<<r1<<" "<<low_integral<<" "<<high_integral<<_xbins_integrals[xbin][ybin+1]<<std::endl;
  }

 
}
