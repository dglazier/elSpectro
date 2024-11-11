//////////////////////////////////////////////////////////////
///
///Class:		Writer
///Description:
///             Interface to different output formats


#pragma once

#include "DecayModel.h"

#include <TObject.h> //for ClassDef

namespace elSpectro{
  
  class Writer {
    
    
  public:

    Writer()=default;
    virtual ~Writer()=default;
    Writer(const Writer& other); //need the virtual destructor...so rule of 5
    Writer(Writer&&)=default;
    Writer& operator=(const Writer& other);
    Writer& operator=(Writer&& other) = default;

    
    
    virtual void Init(const particle_ptrs& iptrs);
    virtual void InitEvent(const particle_ptrs& iptrs,const particle_ptrs& sptrs,const std::vector<const LorentzVector*>& vers);
    virtual void WriteHeader()=0;
    virtual void FillAnEvent()=0;
    virtual void Write()=0;
    virtual void End()=0;


  protected :
    
    particle_ptrs _initialParticles;
    particle_ptrs _finalParticles;
    std::vector<const LorentzVector*> _vertices;
     
    ClassDef(elSpectro::Writer,1); //class Writer
    
  };
}

