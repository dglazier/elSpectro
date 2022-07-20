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

    double _Wmax=0;
    double _Wmin=0;

  };

  class ReactionPhotoProd : public ReactionInfo {


  public:
    virtual ~ReactionPhotoProd()=default;

    LorentzVector* _photon={nullptr};   
    LorentzVector* _target={nullptr};   
    LorentzVector* _photoN{nullptr};   
    LorentzVector* _meson={nullptr};   
    LorentzVector* _baryon={nullptr};
    const LorentzVector* _ebeam={nullptr}; //beam electron   

    PhotonPolarisationVector* _photonPol={nullptr};
    
    mutable double _sWeight = {1}; //s=W^2 excitation function weight
  };

  class ReactionElectroProd : public ReactionPhotoProd {

 
  public:
    virtual ~ReactionElectroProd()=default;

    LorentzVector* _scattered={nullptr}; //scattered electron   
    //Distribution* _Wdist={nullptr}; //W photoproduction
  };


}
