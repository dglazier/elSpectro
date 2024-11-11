#include "Writer.h"
#include "Manager.h"

namespace elSpectro{

  void Writer::InitEvent(const particle_ptrs& iptrs,const particle_ptrs& sptrs,const std::vector<const LorentzVector*>& vers){
    
    _finalParticles.clear();
    _initialParticles.clear();
    _vertices.clear();
    
    //get copies of the particle pointers
 
    for(auto* p:sptrs){
      _finalParticles.push_back(p);
    }

    
    for(auto* p:iptrs){
      _initialParticles.push_back(p);
    }

    for(auto* v:vers){
      _vertices.push_back(v);
    }

    
  }
  void Writer::Init(const particle_ptrs& iptrs){
    _initialParticles.clear();
      for(auto* p:iptrs){
	_initialParticles.push_back(p);
      }
  }
}
