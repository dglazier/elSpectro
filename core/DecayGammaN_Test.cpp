#include "DecayGammaN_Test.h"


namespace elSpectro{

  DecayGammaN_Test::DecayGammaN_Test( particle_ptrs ps, const std::vector<int> pdgs):
    DecayModel{ps,pdgs},_massModel{TF1("hh","TMath::BreitWigner(x,0.78,0.1)",0.3,3.5)}
  {

    _name={"DecayGammaN_Test"};

    if(Products()[0]->Pdg()==2212){
      _meson= Products()[1]; //
      _proton=Products()[0]; //
    }
    else   if(Products()[1]->Pdg()==2212){
      _meson= Products()[0]; //
      _proton=Products()[1]; //
 
    }
    else {
      std::cerr<<" DecayGammaN_Test::DecayGammaN_Test no proton"<<std::endl;
      exit(0);
    }
      

    
  }

}
