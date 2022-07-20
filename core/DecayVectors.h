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
#include <Math/AxisAngle.h>
#include <Math/Vector3Dfwd.h>
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

    virtual double Probability() const {return 1;}



   virtual void BoostToParentWithRandPhi(const LorentzVector& parent, LorentzVector& child){
     if(parent.P()==0){ return;} //no boost to be done, or direction
     //std::cout<<"BoostToParentW "<<child<<" "<<child.M()<<" "<<parent.M()<<std::endl;
   
    //Need axis to rotate around to align
     //z-axis with with parent vector
     //This is just  parent X z-axis
     auto axis1 = parent.Vect().Cross(ROOT::Math::XYZVector(0,0,1)).Unit();
     //align z axis by rotating by parent theta
     ROOT::Math::AxisAngle rot1(axis1,-parent.Theta());
     child=rot1*child;
     //random rotation around parent axis = phi angle of decay
     ROOT::Math::AxisAngle rotAroundParent(parent.Vect().Unit(),RandomPhi());
     child=  rotAroundParent* child; 
 
     //auto check2 = LorentzVector(0,0,1,1);
     //check2=rot1*check2;
     //std::cout<<" BoostToParentWithRandPhi "<<check2.Vect().Unit()<<" "<<parent.Vect().Unit()<<"      axis "<<axis1<<" rot "<<rot1<<std::endl;

     //boost from parent rest frame 
     auto boostFromParent=-parent.BoostToCM();
     child=boost(child,boostFromParent);//ROOT::Math::VectorUtil::boost;
     //std::cout<<"Done BoostToParentW "<<child<<" "<<child.M()<<std::endl;
   }

  protected:
    
    mutable double _weight={1};
    virtual double RandomPhi() const noexcept { return gRandom->Uniform(-TMath::Pi(),TMath::Pi()); }
 
     
  private:
    
   
    ClassDef(elSpectro::DecayVectors,1); //class DecayVectors
 

  };

}
  
