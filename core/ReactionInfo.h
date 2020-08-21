//////////////////////////////////////////////////////////////
///
///Class:		ReactionInfo
///Description:
///             Interface to any reaction dependent info
///             For example pointers to Lorentz Vectors
#pragma once

#include "LorentzVector.h"
#include "PhotonPolarisationVector.h"

namespace elSpectro{
  
  class ReactionInfo{

 
  public:

    virtual ~ReactionInfo()=default;

    
  };

  class ReactionPhotoProd : public ReactionInfo {


  public:
    virtual ~ReactionPhotoProd()=default;

    LorentzVector* _photon={nullptr};   
    LorentzVector* _target={nullptr};   
    LorentzVector* _photoN{nullptr};   
    LorentzVector* _meson={nullptr};   
    LorentzVector* _baryon={nullptr};

    PhotonPolarisationVector* _photonPol={nullptr};
    
    mutable double _sWeight = {1}; //s=W^2 excitation function weight
    
  };

  class ReactionElectroProd : public ReactionPhotoProd {

 
  public:
    virtual ~ReactionElectroProd()=default;

    LorentzVector* _scattered={nullptr}; //scattered electron   
    LorentzVector* _ebeam={nullptr}; //beam electron   
    
  };


}
