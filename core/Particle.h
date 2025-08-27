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
  class DistFlatMassMaster; //so can make friend
 
  enum class DistType {kMass, kMassSquared};
  enum class DecayType{ Stable, Attached, Detached, Production };

  class Particle {

  public:

 
    Particle()=default;
    virtual ~Particle()=default;
    Particle(const Particle& other)=default; //need the virtual destructor...so rule of 5
    Particle(Particle&&)=default;
    Particle& operator=(const Particle& other)=default;
    Particle& operator=(Particle&& other) = default;

    Particle(int pdg);
    
    virtual bool IsDecaying() const {return false;}

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
      _dynamicMass=_vec.M();
    }
    void TakeMaximumMass(){
      SetP4M( MaximumMassPossible() );
    }
   void TakeMinimumMass(){
      SetP4M( MinimumMassPossible() );
    }
    void TakePdgMass(){
      SetP4M( PdgMass() );
    }

    void Boost(const  elSpectro::BetaVector& vboost ){
      _vec=ROOT::Math::VectorUtil::boost(_vec,vboost);
    }

    double M2() const {
      if(Pdg()==22) return 0.;
      return _dynamicMass*_dynamicMass;
    }
    double Mass() const {
      if(Pdg()==22) return 0.;
      return _dynamicMass;
    }
    
    void SetMassDist(std::shared_ptr<Distribution> dist){
      _massDist = dist;
    }
   
    Distribution* MassDistribution() const{return _massDist.get();}
    
    void SetParent(Particle* parent);
 
    void SetPdgMass(double val){ _pdgMass=val; SetP4M(val); }
    
    double PdgMass()const  noexcept{
      return _pdgMass;
    }

    virtual double MinimumMassPossible()const  noexcept{
      // std::cout<<"Particle::MinimumMassPossible() "<< _pdgMass<<std::endl;
      return  PdgMass();
    }
    virtual double MaximumMassPossible()const  noexcept{
      return  PdgMass();
    }
    virtual double MinimumMassForChannel() const  noexcept{
      return MinimumMassPossible();
    }
 
    double MassWeight() const noexcept {
      return _massWeight;
    }

    virtual void Print()  const;

    void SetVertexID(int vertexID){
      _vertexID=vertexID;
    }
    void SetVertexPosition(const LorentzVector& v){
      _vertex=v;
    }
    // void SetVertex(int vertexID,const LorentzVector* v){
    //   _vertexID=vertexID;
    //   _vertex=v;
    //  }
    const LorentzVector& VertexPosition()const noexcept{return _vertex;}
    int VertexID()const noexcept{return _vertexID;}

    virtual DecayType IsDecay() const noexcept {return DecayType::Stable;}
  

    SDME* InitSDME(uint J,uint alphaMax){
      _sdme=SDME(J,alphaMax);
      return &_sdme;
    }
    const SDME* GetSDME() const noexcept{ return &_sdme; }

    void LockMass(){_massLocked=true;}
    void UnlockMass(){_massLocked=false;}
    
    void SetP4M(double mm){
      auto P2=_vec.P2();
      _vec.SetXYZT(_vec.X(),_vec.Y(),_vec.Z(),TMath::Sqrt(P2+mm*mm));
      _dynamicMass=mm;
    }

  private:

    friend DecayModel; //for  DetermineDynamicMass()
    friend DistFlatMassMaster; //for  DetermineDynamicMass()
    
    //if mass comes from a distribution sample it
    void  DetermineDynamicMass(double xmin=-1,double xmax=-1){
      // std::cout<<"Particle::DetermineDynamicMass "<<Pdg()<<" "<<_dynamicMass<<" "<<_massDist<<" "<<xmin<<" "<<xmax<<std::endl;
      if(_massDist==nullptr ){
	TakePdgMass();
	return; //stick at pdgMass
      }
      if(_massLocked==true) return; //someone else in charge...
      _dynamicMass=-1;
      _massWeight=0;
      auto minposs = MinimumMassPossible();
   
      auto minRange = (xmin==-1)?minposs:xmin;
      if(minRange>minposs)minRange=minposs;
      auto maxRange = (xmax==-1)?_massDist->GetMaxX():xmax;
      if((maxRange-minRange)<0) {
	if((maxRange-minRange)>-1E-6) {
	  _dynamicMass=minRange;
	  return;
	}
      }

      if(minRange>maxRange){//unphysical
	std::cout<<"Warning  Particle::DetermineDynamicMass min "<<minRange<<" greater than max "<<maxRange<<" for "<<_pdg<<" minposs "<<minposs<<" "<<xmax<<" "<<_massDist->GetMaxX()<<" equal "<<(minposs==_massDist->GetMaxX())<<std::endl;
	_dynamicMass=minRange; 
	return ;
	//	exit(0);
	 
      }
      while(_dynamicMass<minposs){
	// if((maxRange-minRange)<1E-6) {
	//   _dynamicMass=minRange;
	//   return;
	// }
	//        std::cout<<_pdg<<"  DetermineDynamicMass( "<<MinimumMassPossible()<<" "<<_dynamicMass<<" "<<_massWeight<<" "<<minRange<<" "<<maxRange<<" check "<<minposs-minRange<<"check "<<_massDist->GetMinX()-minRange<<std::endl;

	_dynamicMass= _massDist->SampleSingle(minRange,maxRange);
	

	//std::cout<<"DONE "<<_pdg<<"  DetermineDynamicMass( "<<_dynamicMass<<" "<<MinimumMassPossible()<<" diff "<<_dynamicMass-minRange<<" "<<maxRange-minRange<<std::endl;
	//need a weight for "envelope"
	_massWeight =_massDist->GetCurrentWeight();

	if(_dynamicMass==0) {
	  std::cout<<"Error  Particle::DetermineDynamicMass zero mass"<<std::endl;
	  std::cout<<_pdg<<"  DetermineDynamicMass( "<<MinimumMassPossible()<<" "<<_dynamicMass<<" "<<_massWeight<<" "<<minRange<<" "<<maxRange<<std::endl;
	  exit(0);
	}
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
    LorentzVector _vertex;
    
    std::shared_ptr<Distribution> _massDist={nullptr};
    bool _massLocked={false};
    
    
    ClassDef(elSpectro::Particle,1); //class Particle
    
  };//class Particle

  using particle_uptr = std::unique_ptr<Particle>;

  
}//namespace elSpectro
