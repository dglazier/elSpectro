#include "DecayChannel.h"

namespace elSpectro{

  uint DecayChannel::ChooseDecay(double W) const{
    if(N()==1) return _idecay;
    //choose random number between 0 and sum of all branch ratios
    //find index corresponding to that value
    //std::cout<<"DecayChannel::ChooseDecay W " <<W<< " thresh "<<Threshold()<<" max thresh "<<_maxThreshold<<std::endl;
    if(W>_maxThreshold){  //no need to worry, all channels allowed
      return _brRatioSum.size() == 1 ? _idecay :
	(_idecay = (std::lower_bound(_brRatioSum.begin(),_brRatioSum.end(),gRandom->Uniform(0,_brRatioSum.back()))) - _brRatioSum.begin() -1) ; //return 0 if 1 decay, if not choose.
    }
    /////////////////////////////////////////////////////
    //Not all decays above threshold, shorten choices
    std::vector<double> physicalBRSums={0.};
    //create branching ratio sums with threshold above W
    std::vector<short> shiftIdx;
    for(auto i=0;i<N();++i){
      if(W>_models[i]->MinimumMassPossible()){
	//ideally we should adjust the br depending on how far above threshold we are
  	physicalBRSums.push_back(physicalBRSums.back()+_brRatios[i]);
	shiftIdx.push_back(i);
	//	std::cout<<"DecayChannel::ChooseDecay add channel "<<i<<" "<<_brRatios[i]<<" "<<W<<" > threshold "<<_models[i]->MinimumMassPossible()<<std::endl;
      }
    }
    if(physicalBRSums.size()==1){
      std::cerr<<"DecayChannel::ChooseDecay , no decays have sufficent mass. Should not be here. W = "<<W<<std::endl;exit(0);
    }
     _idecay = shiftIdx[(std::lower_bound(physicalBRSums.begin(),physicalBRSums.end(),gRandom->Uniform(0,physicalBRSums.back()))) - physicalBRSums.begin() -1];
     //std::cout<<"DecayChannel::ChooseDecay number of channels " <<physicalBRSums.size()-1<<" "<<physicalBRSums.front()<<" "<<physicalBRSums.back()<<" at W "<<W<<" choose "<<_idecay<<std::endl;
    return  _idecay;
  }
}
