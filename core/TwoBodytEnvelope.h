//////////////////////////////////////////////////////////////
///
///Class:		TwoBodytEnvelope
///Description:
///           Provide 2D histogram distribution for W and cosTheta
///            
#pragma once

#include "TwoBodyFlat.h"
#include "DistYGivenX.h"

namespace elSpectro{


  class TwoBodytEnvelope : public TwoBodyFlat {

  public:

    TwoBodytEnvelope()=default;
    TwoBodytEnvelope(const DistYGivenX& dist);

    double RandomCosTh() noexcept final{
      _dist.SetX(W());
      //sample t
      //std::cout<<"TwoBodytEnvelopeRandomCosTh() "<<std::endl;
      _dist.SampleSingle();
      _weight = _dist.GetCurrentWeight();
      auto t = _dist.GetY();
      //convert t to cos(theta)
      auto cosTh=kine::costhFromt(t, W(),0.0,0.938272,_a.M(),_b.M());
      // std::cout<<"TwoBodytEnvelopeRandomCosTh() "<<_dist.GetY()<<" weight = "<<_weight<<" xs "<<_dist.CurrentValue()<<" "<<_dist.MaxValue()<<" W "<<W()<<" cosTh "<<_dist.GetY()<<" "<<(_dist.GetY()>1.0)<<" "<<((y/TMath::Abs(y))>1.0)<<std::endl;
      //protect against rounding errors giving |y|>1
      return TMath::Abs(cosTh)>1.0 ? cosTh/TMath::Abs(cosTh) : cosTh;
    }
   
   
    void PostInit(ReactionInfo* info) final {std::cout<<"TwoBodytEnvelope:: PostInit"<<std::endl; TwoBodyFlat::PostInit(info);};

    const DistYGivenX& GetDist() const {return _dist;}
    void SetDist(const DistYGivenX& dist){_dist = dist;}
    
  private:

    //DistYGivenX _dist={TH2D()};
   DistYGivenX _dist;
     
   ClassDef(elSpectro::TwoBodytEnvelope,1); //class DecayVectors
 

  };

}
  
