//////////////////////////////////////////////////////////////
///
///Class:		Particle
///Description:
///             Control behaviour of particles
///             Particle is defined by
///             1) its instaneous LorentzVector
///             2) any subsequent Decays


#pragma once

#include "LorentzVector.h"
#include "Distribution.h"
#include <TObject.h> //for ClassDef
#include <TMath.h> //for Sqrt
#include <TRandom.h> //for Sqrt
#include <vector>
#include <memory>

namespace elSpectro{
  
  class DecayModel; //so can make friend
 
  enum class DistType {kMass, kMassSquared};

  class Particle {

  public:

 
    Particle()=default;
    //only declaring default constructor
    //so other 5 constructors also defaulted(rule of 5)

    Particle(int pdg);

    const LorentzVector& P4() const {return _vec;}//can be changed by others
    int Pdg()const{return _pdg;}

 
    void SetXYZT(double xx,double yy,double zz, double tt){
      _vec.SetXYZT(xx,yy,zz,tt);
      _dynamicMass=_vec.M();
    }
    
    void SetXYZ(double xx,double yy,double zz){
      auto m2=_vec.M2();auto P2=xx*xx+yy*yy+zz*zz;
      _vec.SetXYZT(xx,yy,zz,TMath::Sqrt(P2+m2));
    }

    void Boost(const  elSpectro::BetaVector& vboost ){
      _vec=ROOT::Math::VectorUtil::boost(_vec,vboost);
    }

    double M2() const {
      return _dynamicMass*_dynamicMass;
    }
    double Mass() const {
      return _dynamicMass;
    }
    
    void SetMassDist(Distribution* dist){
      _massDist=dist;
      //_distType=DistType::kMass;
    }
    const Distribution* MassDistribution() const{return _massDist;}
    
    /* void SetMassSquaredDist(Distribution* dist){ */
    /*   _massDist=dist; */
    /*   _distType=DistType::kMassSquared; */
    /* } */

  
    double PdgMass()const {return _pdgMass;}

    virtual double MinimumMassPossible()const {return PdgMass(); }
    

    virtual void Print() const;
    
  protected:
 

  private:
    friend DecayModel; //for  DetermineDynamicMass()
    
    //if mass comes from a distribution
    void  DetermineDynamicMass(){
      _dynamicMass=-1;
      auto min = MinimumMassPossible();
      while(_dynamicMass<min){
	//	std::cout<<_pdg<<"  DetermineDynamicMass( "<<MinimumMassPossible()<<" "<<_dynamicMass<<std::endl;
	_dynamicMass=
	  (_massDist==nullptr ? _pdgMass : _massDist->SampleSingle());
      }
      
     }

    LorentzVector _vec;
    double _pdgMass={0};
    double _dynamicMass={0};
    
    int _pdg={0};

    Distribution* _massDist={nullptr};

    // DistType _distType={DistType::kMass};
    
    ClassDef(elSpectro::Particle,1); //class Particle
    
  };//class Particle

  using particle_uptr = std::unique_ptr<Particle>;

  
}//namespace elSpectro
