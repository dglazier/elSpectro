#include "Writer.h"
#include "Manager.h"

namespace elSpectro{

  void Writer::Init(){
    
     //get copies of the particle pointers
    auto& sptrs=Manager::Instance().Particles().StableParticles();
    for(const auto* p:sptrs){
      _finalParticles.push_back(p);
     }
    auto& iptrs=Manager::Instance().Reaction()->InitialParticles();
    for(const auto* p:iptrs)
      _initialParticles.push_back(p);

    _vertices = (&(Manager::Instance().GetVertices()));
   
  }
}
