//////////////////////////////////////////////////////////////
///
///Class:		DecayGammaN_Test
///Description:
///             Calculate intensity for phase space decay (=1 )!

#pragma once

#include "DecayModel.h"
#include "DistTF1.h"
#include "FunctionsForElectronScattering.h"
#include <TMath.h>

namespace elSpectro{

 
  class DecayGammaN_Test : public DecayModel {

  public:
    
    DecayGammaN_Test()=default;
    //only declaring default constructor
    //so other 5 constructors also defaulted(rule of 5)
    //constructor to decay into particles
    DecayGammaN_Test( particle_ptrs , const std::vector<int> pdgs );

    // Each model must define its intensity
    double Intensity() const final{

      double weight=_massModel.GetWeightFor(_meson->Mass())/_meson->MassWeight();
      //myInfo._weight=_meson->MassReWeight(_massModel);
      
      //std::cout<< _meson->Mass()<<" "<<_massModel.GetWeightFor(_meson->Mass())<<" "<<_meson->MassWeight()<<std::endl;
      //double t = (_target - _proton->P4()).M2();
      
      //_myInfo._weight=TMath::Exp( 10*t );
      // _myInfo._weight=1;
      return weight;
    }

    bool RegenerateOnFail() const noexcept final {return false;}

  private:
 
    //mutable CurrentEventInfo _myInfo;//!

    DistTF1 _massModel;

    Particle* _meson={nullptr};
    Particle* _proton={nullptr};
    LorentzVector _target={0,0,0,escat::M_pr()};
    
    ClassDef(elSpectro::DecayGammaN_Test,1); //class DecayGammaN_Test
    
  };//class DecayGammaN_Test

}//namespace elSpectro
