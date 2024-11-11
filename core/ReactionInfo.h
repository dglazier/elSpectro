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

  class ProductionProcess;
  
  class ReactionInfo{

 
  public:

    virtual ~ReactionInfo()=default;

    // double _Wmax=0;
    double _Wmin=0;
    //  ProductionProcess* _process={nullptr};
    virtual double Wmax() const {return 0;}
  };

  class ReactionPhotoProd : public ReactionInfo {


  public:
    virtual ~ReactionPhotoProd()=default;

    mutable LorentzVector _photon;   
    mutable LorentzVector _target;   
    // mutable LorentzVector _photoN;   
    mutable LorentzVector _meson;   
    mutable LorentzVector _baryon;
    mutable LorentzVector _ebeam; //beam electron   

    mutable PhotonPolarisationVector _photonPol;
    
    mutable double _sWeight = {1}; //s=W^2 excitation function weight

    double Wmax() const override {return (_target + _ebeam).M();}
  };

  class ReactionElectroProd : public ReactionPhotoProd {

 
  public:
    virtual ~ReactionElectroProd()=default;

    mutable LorentzVector _scattered; //scattered electron   
    //Distribution* _Wdist={nullptr}; //W photoproduction
  };


}
