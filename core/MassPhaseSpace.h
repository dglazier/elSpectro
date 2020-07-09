//////////////////////////////////////////////////////////////
///
///Class:		MassPhaseSpace
///Description:
///            Class to manage full reaction phase space
///            Should only be accessed via Manager::PhaseSpace;
///            Needs to be given the primary decay to calculate
///            the Mass phase space element for all children
///            and allocate masses for all particles in the chain
#pragma once

#include "DecayModel.h"
#include <TRandom.h>

namespace elSpectro{

  class Manager;
  
  class MassPhaseSpace{

  public:
    
  private:
    
    friend Manager; //only Manager can construct and use a MassPhaseSpace
    MassPhaseSpace()=default;

    double PhaseSpaceWeight(double parentM) const noexcept{
      double result=TMath::Sqrt(_model->PhaseSpaceWeightSq(parentM));
      // std::cout<<"mass phase space "<<result<<std::endl;
      return result;
    }

    void Find(double parentM,const  DecayModel* amodel){
      //sample all masses according to overall decay phase space
      if(_model==nullptr) return;
      if(_model!=amodel) return; //only 1 model controls phasespace

      //in case decay chain may change each event coulsd get the masses each time
      // std::vector<double> masses;
      //_model->GetStableMasses(masses); //fill vector with masses of final particles

      double max= kine::PhaseSpaceWeightMax(parentM,_masses);
      //  std::cout<<"MassPhaseSpace() max "<< max<<std::endl;

      //accept or reject mass combinations until got one
      //PhaseSpaceWeight will try alternative masses
      while( PhaseSpaceWeight(parentM) < gRandom->Uniform()*max ){
	//reject this combintation
      }
    
    }
    
    void SetModel(DecayModel* amodel){
      _model=amodel;
      _model->GetStableMasses(_masses); //fill masses vector
    }
    
 
    DecayModel* _model=nullptr;
    std::vector<double> _masses;
    
    ClassDef(elSpectro::MassPhaseSpace,1); //class MassPhaseSpace
  };

}//namespace elSpectro
