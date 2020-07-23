//////////////////////////////////////////////////////////////
///
///Class:		ReactionInfo
///Description:
///             Interface to any reaction dependent info
///             For example pointers to Lorentz Vectors
#pragma once

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

    mutable double _sWeight = {1};
    
  };

  class ReactionElectroProd : public ReactionPhotoProd {

 
  public:
   virtual ~ReactionElectroProd()=default;

    LorentzVector* _scattered={nullptr}; //scattered electron   
    LorentzVector* _ebeam={nullptr}; //beam electron   
    
  };


}
