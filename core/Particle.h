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
#include "SDME.h"
#include <TObject.h> //for ClassDef
#include <TMath.h> //for Sqrt
#include <TRandom.h> //for Sqrt
#include <vector>
#include <memory>

namespace elSpectro{
  
  class DecayModel; //so can make friend
 
  enum class DistType {kMass, kMassSquared};
  enum class DecayType{ Stable, Attached, Detached, Production };

  class Particle {

  public:

 
    Particle()=default;
    virtual ~Particle()=default;
    Particle(const Particle& other); //need the virtual destructor...so rule of 5
    Particle(Particle&&)=default;
    Particle& operator=(const Particle& other);
    Particle& operator=(Particle&& other) = default;

    Particle(int pdg);

    LorentzVector const& P4() const {return _vec;}//const be changed by others
    LorentzVector* P4ptr() {return &_vec;}
    
    int Pdg()const{return _pdg;}

 
    void SetXYZT(double xx,double yy,double zz, double tt){
      _vec.SetXYZT(xx,yy,zz,tt);
      _dynamicMass=_vec.M();
    }
    
    void SetXYZ(double xx,double yy,double zz){
      auto m2=_vec.M2();auto P2=xx*xx+yy*yy+zz*zz;
      _vec.SetXYZT(xx,yy,zz,TMath::Sqrt(P2+m2));
    }
    void SetP4(const LorentzVector& p4){
      _vec=p4;
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

  
    double PdgMass()const  noexcept{return _pdgMass;}

    virtual double MinimumMassPossible()const  noexcept{return PdgMass(); }
    
    double MassWeight() const noexcept {
      return _massWeight;
    }

    virtual void Print()  const;

    void SetVertex(int vertexID,const LorentzVector* v){
      _vertexID=vertexID;
      //_vertex.SetXYZT(v.X(),v.Y(),v.Z(),v.T());
      _vertex=v;
    }
    const LorentzVector* VertexPosition()const noexcept{return _vertex;}
    int VertexID()const noexcept{return _vertexID;}

    virtual DecayType IsDecay() const noexcept {return DecayType::Stable;}
  

    SDME* InitSDME(uint J,uint alphaMax){
      _sdme=SDME(J,alphaMax);
      return &_sdme;
    }
    const SDME* GetSDME() const noexcept{ return &_sdme; }
    
  protected:
 
    void SetP4M(double mm){
      auto P2=_vec.P2();
      _vec.SetXYZT(_vec.X(),_vec.Y(),_vec.Z(),TMath::Sqrt(P2+mm*mm));
    }

  private:

    friend DecayModel; //for  DetermineDynamicMass()
    
    //if mass comes from a distribution
    void  DetermineDynamicMass(){
      
      if(_massDist==nullptr ) return; //stick at pdgMass
      
      _dynamicMass=-1;
      auto min = MinimumMassPossible();
      while(_dynamicMass<min){
	
	_dynamicMass= _massDist->SampleSingle();
	//need a weight for "envelope"
	_massWeight =_massDist->GetCurrentWeight();
	//	std::cout<<_pdg<<"  DetermineDynamicMass( "<<MinimumMassPossible()<<" "<<_dynamicMass<<" "<<_massWeight<<std::endl;
	
      }
      SetP4M(_dynamicMass);

    }

    
    LorentzVector _vec;
    SDME _sdme;
    double _pdgMass={0};
    double _dynamicMass={0};
    double _massWeight={1};
    
    int _pdg={0};
    int _vertexID={0};
    const LorentzVector* _vertex={nullptr};
    
    Distribution* _massDist={nullptr};

    // DistType _distType={DistType::kMass};
    
    ClassDef(elSpectro::Particle,1); //class Particle
    
  };//class Particle

  using particle_uptr = std::unique_ptr<Particle>;

  
}//namespace elSpectro
