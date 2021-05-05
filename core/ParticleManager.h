//////////////////////////////////////////////////////////////
///
///Class:		ParticleManager
///Description:
///            Class to manage particles
///            Should only be accessed via Manager::Particles()
///             1) Keep ownership of all particles
///             2) Construct particles via Create
///             3) Control output file particles
#pragma once

#include "DecayingParticle.h"
#include "CollidingParticle.h"
#include "Distribution.h"
#include <Math/RotationZYX.h>
#include <Math/RotationZ.h>
#include <map>

namespace elSpectro{

  class Manager;
  
  using particle_ptrs= std::vector<Particle*> ;
  using particle_constptrs= std::vector<const Particle*> ;
  using decaying_ptrs= std::vector<DecayingParticle*> ;
  using decaying_constptrs= std::vector<const DecayingParticle*> ;
  using colliding_ptrs= std::vector<CollidingParticle*> ;
  using colliding_constptrs= std::vector<const CollidingParticle*> ;

  using ROOT::Math::VectorUtil::boost;

  class ParticleManager{

  public:

    //take ownership of a particle
    //do not use p after calling this
    Particle* Take(Particle* p);

    void RegisterMassDistribution(int pdg, Distribution* dist){
      _massDist[pdg]=dist_uptr{dist};
      dist=nullptr;
      //Check if any predefined particles need this distribution
      for(auto& particle:_particles){
	if(pdg==particle->Pdg()){
	  particle->SetMassDist( _massDist.at(pdg).get());
	}
      }

    }
    
    int RegisterNewPdgParticle(double nominalMass,Distribution* dist=nullptr){
      std::cout<<"RegisterNewPdgParticle "<<" "<<_nextPdg<<std::endl;
      AddToPdgTable(_nextPdg,nominalMass);
      if(dist != nullptr){
	RegisterMassDistribution(_nextPdg,dist);
      }
      _nextPdg++;
      return _nextPdg-1;
    }

    void AddToPdgTable(int pdg,double mass);
    Double_t GetMassFor(int pdg);
    /* void RegisterMassSquaredDistribution(int pdg, Distribution* dist){ */
    /*   _mass2Dist[pdg]=dist_uptr{dist}; */
    /* } */
    Distribution* GetMassDist(int pdg)const noexcept
    {return _massDist.at(pdg).get();}

    void BoostStable(const BetaVector& vboost ){
      std::for_each(_stables.begin(),_stables.end(),[&vboost](Particle* p){p->Boost(vboost);});
    }

    void MoveStableToLab(Particle* particle){
      //So the particle does not get boosted
      //but is written out to final state
      _stables.erase(std::remove(_stables.begin(),_stables.end(),particle),_stables.end());
      _stableslab.push_back(particle);
    }
   void RemoveStable(Particle* particle){
      //So the particle does not get boosted
      //or written out to final state
      _stables.erase(std::remove(_stables.begin(),_stables.end(),particle),_stables.end());
    }
    
    const decaying_ptrs UnstableParticles()const {return _unstables;}
    const particle_ptrs StableParticles()const {
      particle_ptrs allstables;
      allstables.insert(std::end(allstables), std::begin(_stables), std::end(_stables));
      allstables.insert(std::end(allstables), std::begin(_stableslab), std::end(_stableslab));
      return allstables;
    }
    
    void BoostToFrame(const BetaVector& vboost,const LorentzVector& parent){

      //vboost defines  z-axis
      /*if(_cachedBoost!=vboost){ //SetAngle is expensive (sin,cos calls) only call if necessary
	_cachedBoost = vboost;
	_rotateToZaxis.SetComponents(-_cachedBoost.Phi(),-_cachedBoost.Theta(),0);
	}*/
      // _rotateToZaxis.SetComponents(0,-(TMath::Pi()-parent.Theta()),0);  
      //Apply random phi angle when z-axis is in correct direction
      //_rotateAroundZaxis.SetAngle(0);

  
      for(auto& child : _stables){
	auto p4=child->P4();
	p4=_rotateToZaxis*p4;
	//p4=_rotateAroundZaxis*p4; 
	p4=boost(p4,vboost);//ROOT::Math::VectorUtil::boost;
	child->SetP4(p4);
      }
      
    }
  private:
    
    friend Manager; //only Manager can construct a ParticleManager
    ParticleManager();
    
    //all the particles in the generator
    std::vector<particle_uptr> _particles;
    
    particle_ptrs _stables; //products which are stable (to be detected)
    particle_ptrs _stableslab; //products which are stable and in lab frame
    decaying_ptrs _unstables; //products which decay
    colliding_ptrs _initials; //intial interactin particles

    //distribitions for generating dynamic particle masses
    std::map<int,dist_uptr> _massDist; 
    // std::map<int,dist_uptr> _mass2Dist; 

    int _nextPdg={10000};

    //For boost to lab
    BetaVector _cachedBoost;
    ROOT::Math::RotationZYX _rotateToZaxis; //save memory allocation
    ROOT::Math::RotationZ _rotateAroundZaxis;//save memory allocation
 
    
    ClassDef(elSpectro::ParticleManager,1); //class ParticleManager
  };

}//namespace elSpectro
