//////////////////////////////////////////////////////////////
///
///Class:		Formation
///Description:
///            Control production and decay of state formed by
///            initial state reaction, e.g. e- scattering or photoprodction
///            this is equivalent to the gamma+N centre of mass system
///            which decays to meson + baryon final state

#pragma once

#include "DecayChannel.h"
#include "DecayingParticle.h"
#include "ExcitationSpectra.h"
#include "ReactionInfo.h"
#include "PhotonPolarisationVector.h"
#include "DistTH1.h"
#include <TH1D.h>

namespace elSpectro{

  
  class Formation : public DecayModel {

  public:
    
    Formation()=default;
    virtual ~Formation()=default;

    Formation( double thresh, const decaying_objs& decs, const particle_objs& stables)
      :DecayModel(decs,stables),   _threshold{thresh} {}
    
 
    bool RegenerateOnFail() const  noexcept override {return false;}

    const DecayingParticle& GetGammaN() const noexcept{return *static_cast<const DecayingParticle*>(Product(_idxGStarNuc)); }
   
    double getW() const noexcept{ return GetGammaN().Mass();}
    double getThreshold() const noexcept{return _threshold;}
    
    void setThreshold(double val) const noexcept{
      if(val<_threshold) return;
      _threshold=val;
    }
    

    double dsigma() const override { return GetGammaN().Model()->dsigma();}// * Q2 factor }

    //Any preliminaries required
    bool ReadyForDecay() override{
      ChooseDecay(); //need to choose from all TwoBodyProductions based on CrossSection
      return true;
      // or return CheckThreshold(); ?
    }

    const DecayChannel& Channels() const {
      //auto dp =  dynamic_cast<DecayingParticle*>( _gstarNuc );
      //if(dp) return &(dp->Channels());
      //return nullptr;
      return  GetGammaN().Channels();
    }
    
  protected:
    
    uint _idxGStarNuc=0;
    mutable LorentzVector _gamma;
    mutable PhotonPolarisationVector _photonPol;
    mutable double _threshold = {0};
    ReactionInfo* _reactInfo=nullptr;
    const ReactionInfo* GetReactionInfo() const noexcept {return _reactInfo;}

  public:
    
    void ZeroPhoton() const noexcept{
      _gamma.SetXYZT(0,0,0,0);
    }

    const DistTH1& TotalCrossSection()const  noexcept{return _spectra.TotalCrossSection();}

    const ExcitationSpectra& ExciteSpectra() const noexcept{return _spectra;}

    void CreateExcitationSpectra(){
      _spectra.CreateSpectra( Channels(), _reactInfo->Wmax());
    }
    
  private:


    ExcitationSpectra _spectra;
    DistTH1 _distTotalXS;
 
    ClassDefOverride(elSpectro::Formation,1); //class Formation
    
  };//class Formation

}//namespace elSpectro
