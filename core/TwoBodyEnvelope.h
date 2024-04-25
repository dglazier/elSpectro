//////////////////////////////////////////////////////////////
///
///Class:		TwoBodyEnvelope
///Description:
///           Provide 2D histogram distribution for W and cosTheta
///            
#pragma once

#include "TwoBodyFlat.h"
#include "DistTH2Slice.h"

namespace elSpectro{


  class TwoBodyEnvelope : public TwoBodyFlat {

  public:

  //   virtual ~TwoBodyEnvelope()=default;
  //   TwoBodyEnvelope(const TwoBodyEnvelope& other)=default; //need the virtual destructor...so rule of 5
  //   TwoBodyEnvelope(TwoBodyEnvelope&&)=default;
  //   TwoBodyEnvelope& operator=(const TwoBodyEnvelope& other)=default;
  //   TwoBodyEnvelope& operator=(TwoBodyEnvelope&& other) = default;
 
  TwoBodyEnvelope()=default;
  TwoBodyEnvelope(const DistTH2Slice& dist);

    double RandomCosTh() noexcept final{
      _dist.SetX(W());
      _dist.SampleSingle();
      _weight =  _dist.GetCurrentWeight();
      //  std::cout<<"TwoBodyEnvelopeRandomCosTh() "<<_dist.GetY()<<" weight = "<<_weight<<" xs "<<_dist.CurrentValue()<<" "<<_dist.MaxValue()<<" W "<<W()<<" cosTh "<<_dist.GetY()<<std::endl;
      return _dist.GetY();
    }
   
   
    void PostInit(ReactionInfo* info) final {std::cout<<"TwoBodyEnvelope:: PostInit"<<std::endl; TwoBodyFlat::PostInit(info);};

    const DistTH2Slice& GetDist() const {return _dist;}
    void SetDist(const DistTH2Slice& dist){_dist = dist;}
    
  private:

   DistTH2Slice _dist={TH2D()};
     
   ClassDef(elSpectro::TwoBodyEnvelope,1); //class DecayVectors
 

  };

}
  
