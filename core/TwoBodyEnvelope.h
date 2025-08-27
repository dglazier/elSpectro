//////////////////////////////////////////////////////////////
///
///Class:		TwoBodyEnvelope
///Description:
///           Provide 2D histogram distribution for W and cosTheta
///            
#pragma once

#include "TwoBodyFlat.h"
#include "DistYGivenX.h"

namespace elSpectro{


  class TwoBodyEnvelope : public TwoBodyFlat {

  public:

  //   virtual ~TwoBodyEnvelope()=default;
  //   TwoBodyEnvelope(const TwoBodyEnvelope& other)=default; //need the virtual destructor...so rule of 5
  //   TwoBodyEnvelope(TwoBodyEnvelope&&)=default;
  //   TwoBodyEnvelope& operator=(const TwoBodyEnvelope& other)=default;
  //   TwoBodyEnvelope& operator=(TwoBodyEnvelope&& other) = default;
 
  TwoBodyEnvelope()=default;
  TwoBodyEnvelope(const DistYGivenX& dist);

    double RandomCosTh() noexcept final{
      std::cout<<"TwoBodyEnvelopeRandomCosTh() "<<std::endl;
      _dist.SetX(W());
      _dist.SampleSingle();
      _weight = _dist.GetCurrentWeight();
      auto y = _dist.GetY();
      // std::cout<<"TwoBodyEnvelopeRandomCosTh() "<<_dist.GetY()<<" weight = "<<_weight<<" xs "<<_dist.CurrentValue()<<" "<<_dist.MaxValue()<<" W "<<W()<<" cosTh "<<_dist.GetY()<<" "<<(_dist.GetY()>1.0)<<" "<<((y/TMath::Abs(y))>1.0)<<std::endl;
      //protect against rounding errors giving |y|>1
      return TMath::Abs(y)>1.0 ? y/TMath::Abs(y) : y;
    }
   
   
    void PostInit(ReactionInfo* info) final {std::cout<<"TwoBodyEnvelope:: PostInit"<<std::endl; TwoBodyFlat::PostInit(info);};

    const DistYGivenX& GetDist() const {return _dist;}
    void SetDist(const DistYGivenX& dist){_dist = dist;}
    
  private:

    //DistYGivenX _dist={TH2D()};
   DistYGivenX _dist;
     
   ClassDef(elSpectro::TwoBodyEnvelope,1); //class DecayVectors
 

  };

}
  
