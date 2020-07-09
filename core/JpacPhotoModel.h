//////////////////////////////////////////////////////////////
///
///Class:		JpacPhotoModel
///Description:
///             Control behaviour of Particle decay to Particle products
///             Defined by
///             1) list of Particle products
///             2) Intensity function given the produt LorentzVectors
///
///            Note drevied classes should include a constructor to initialise
///            JpacPhotoModel( particle_ptrs , const std::vector<int> pdgs );
#pragma once

#include "DecayModel.h"

namespace elSpectro{

 
  class JpacPhotoModel : public DecayModel {

  public:
    
    JpacPhotoModel()=default;
    //only declaring default constructor
    //so other 5 constructors also defaulted(rule of 5)

    //constructor to decay into particles
    JpacPhotoModel( particle_ptrs , const std::vector<int> pdgs );
    
    // Each model must define its intensity
    virtual double Intensity() const=0;


    
  private:

    
     
    ClassDef(elSpectro::JpacPhotoModel,1); //class JpacPhotoModel
    
  };//class JpacPhotoModel

}//namespace elSpectro
