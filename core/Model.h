//////////////////////////////////////////////////////////////
///
///Class:		Model
///Description:
///             Control behaviour of particles
///             Model is defined by its instaneous LorentzVector
///             and any subsequent Decays


#pragma once

#include <Math/Vector4D.h> //For XyZTVector
#include <TObject.h> //for ClassDef
#include <vector>

namespace elSpectro{

  
  using LorentzVector= ROOT::Math::XYZTVector;

  class Model {

  public:
    
    Model()=default;
    //only declaring default constructor
    //so other 5 constructors also defaulted(rule of 5)


  private:

    
    
    std::vector<Particle>_children;


    
    ClassDef(elSpectro::Model,1); //class Model
    
  };//class Model

}//namespace elSpectro
