//////////////////////////////////////////////////////////////
///
///Class:		DecayVectors
///Description:
///            Abstract class, derive methods for 
///            Generating LorentzVectors for S->1,2,3...
#pragma once

#include "Particle.h"
#include "DecayModel.h"
#include <Math/RotationZYX.h>
#include <Math/RotationZ.h>
#include <vector>

namespace elSpectro{
  
  using ROOT::Math::VectorUtil::boost;

 
  class DecayVectors {

  public:

    DecayVectors()=default;
    virtual ~DecayVectors()=default;
    DecayVectors(const DecayVectors& other); //need the virtual destructor...so rule of 5
    DecayVectors(DecayVectors&&)=default;
    DecayVectors& operator=(const DecayVectors& other);
    DecayVectors& operator=(DecayVectors&& other) = default;

    ///Function will decay parent to products and return a weight that
    ///can be used to get phase space distribution
    virtual double  Generate(const LorentzVector& parent,
			     const particle_ptrs& products)  = 0;

    
    virtual void PostInit(ReactionInfo* info){};

    virtual double dsigma() const {return 1;}

    void ForIntegrate(bool integ){_forIntegral=integ;}
    virtual double Probability() const {return 1;}

  protected:
    
    mutable double _weight={1};

    virtual void BoostToParent(const LorentzVector& parent, const particle_ptrs& products){
  
      if(_cachedParent!=parent){ //SetAngle is expensive (sin,cos calls) only call if necessary
	_cachedParent = parent;
	//	_rotateToZaxis.SetAngle(_cachedParent.Theta());
	_rotateToZaxis.SetComponents(-_cachedParent.Phi(),-_cachedParent.Theta(),0);
        }
      
      //Apply random phi angle when z-axis is in correct direction
      _rotateAroundZaxis.SetAngle(RandomPhi());
      //boost from parent rest frame 
      auto boostFromParent=-_cachedParent.BoostToCM();

      for(auto& child : products){
	auto p4=child->P4();
	p4=_rotateToZaxis*p4;
	p4=_rotateAroundZaxis*p4; 
	p4=boost(p4,boostFromParent);//ROOT::Math::VectorUtil::boost;
	child->SetP4(p4);
      }
      
    }
   virtual void BoostToParent(const LorentzVector& parent, LorentzVector& child){
  
      if(_cachedParent!=parent){ //SetAngle is expensive (sin,cos calls) only call if necessary
	_cachedParent = parent;
	_rotateToZaxis.SetComponents(-_cachedParent.Phi(),-_cachedParent.Theta(),0);
      }
      
      //Apply random phi angle when z-axis is in correct direction
      _rotateAroundZaxis.SetAngle(RandomPhi());
 
      //boost from parent rest frame 
      auto boostFromParent=-_cachedParent.BoostToCM();
      
      child=_rotateToZaxis * child;
      child=_rotateAroundZaxis * child; 
      child=boost(child,boostFromParent);//ROOT::Math::VectorUtil::boost;
      
      auto check = _cachedParent;
      check=_rotateToZaxis*check;
   }
   virtual void RotateToParent(const LorentzVector& parent, LorentzVector& child){
  
      if(_cachedParent!=parent){ //SetAngle is expensive (sin,cos calls) only call if necessary
	_cachedParent = parent;
	//	_rotateToZaxis.SetAngle(_cachedParent.Theta());
	_rotateToZaxis.SetComponents(-_cachedParent.Phi(),-_cachedParent.Theta(),0);
        }
      
      //Apply random phi angle now z-axis is in correct direction
      _rotateAroundZaxis.SetAngle(RandomPhi());
       
      child=_rotateToZaxis * child;
      child=_rotateAroundZaxis * child; 
      
   }

    virtual double RandomPhi() const noexcept { return gRandom->Uniform(-TMath::Pi(),TMath::Pi()); }
 
  protected:
    
    bool _forIntegral=false;

  private:
    
    LorentzVector _cachedParent;
    ROOT::Math::RotationZYX _rotateToZaxis; //save memory allocation
    ROOT::Math::RotationZ _rotateAroundZaxis;//save memory allocation
 
 
    ClassDef(elSpectro::DecayVectors,1); //class DecayVectors
 

  };

}
  
