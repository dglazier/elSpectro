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
    void Print(){
      std::cout<<"MassPhaseSpace number calcs= "<<_weightCalcN<<" number of successes = "<<_successN<<" ratio  ="<< double(_successN)/_weightCalcN <<std::endl;
    }
    void SuppressPhaseSpace(double val){_suppressPhaseSpace=val;}
  private:
    
    friend Manager; //only Manager can construct and use a MassPhaseSpace
    MassPhaseSpace()=default;

    double PhaseSpaceWeight(double parentM) const noexcept{
      double result=TMath::Sqrt(_model->PhaseSpaceWeightSq(parentM));
      _weightCalcN++;
      return result;
    }

    void Find(double parentM,const  DecayModel* amodel){
      //sample all masses according to overall decay phase space
      if(_model==nullptr) return;
      if(_model!=amodel) return; //only 1 model controls phasespace

      //in case decay chain may change each event coulsd get the masses each time

      // double max= kine::PhaseSpaceWeightMax(parentM,_masses);//TGenPhaseSpace max . Note this is too high an estimate

      //Note PhaseSpaceWeightMaxFromEquDist does not quite get to max, so increase by 10% to be safe....will give warning if event weight is above this....
      double max= kine::PhaseSpaceWeightMaxFromEquDist(parentM,_masses)*1.05;
      //double max= kine::PhaseSpaceWeightMaxFromEquDist(parentM,_masses);
      //accept or reject mass combinations until got one
      //as W dependence accounted for elsewere
      //PhaseSpaceWeight will try alternative masses
      double wee=0;
      while( (wee=PhaseSpaceWeight(parentM)) < gRandom->Uniform()*max*_suppressPhaseSpace )
	{
	  // if(wee>max){
	  // _sampledMax=wee;
	  //  std::cout<<" weight >  max "<<wee<<" "<<max<<" normal max "<<kine::PhaseSpaceWeightMax(parentM,_masses)<<std::endl;
	  // }
	//reject this combintation
	  //	if(wee==0)	std::cout<<"ps "<<wee <<" "<<max<<" W "<<parentM<<" sample max"<<_sampledMax<<std::endl;
      }
      if(wee>max){
	// _sampledMax=wee;
	std::cerr<<"MassPhaseSpace check weight >  max,  W "<<parentM<<" this "<<wee<<" "<<max<<" normal max "<<kine::PhaseSpaceWeightMax(parentM,_masses)<<std::endl;
      }

      _successN++;
    }
    bool AcceptPhaseSpace(double parentM){
      if(_model==nullptr) {
	std::cerr<<"AcceptPhaseSpace  wrong model "<<std::endl;
	exit(0);
      }
      _masses.clear();
      _model->GetStableMasses(_masses); //fill masses vector
      double max= kine::PhaseSpaceWeightMax(parentM,_masses);

      auto weight = PhaseSpaceWeight(parentM);
      std::cout<<"AcceptPhaseSpace "<<parentM<<" "<< weight<<" "<<max<<std::endl;
      
      return ( weight > gRandom->Uniform()*max ) ?
	true : false;  
    }
    
    void SetModel(DecayModel* amodel){
      _model=amodel;
      _masses.clear();
      _model->GetStableMasses(_masses); //fill masses vector
    }
    
 
    DecayModel* _model=nullptr;
    std::vector<double> _masses;
    double _sampledMax = 0;
    mutable long _weightCalcN=0;
    mutable long _successN=0;
    double _suppressPhaseSpace=1;
    
    ClassDef(elSpectro::MassPhaseSpace,1); //class MassPhaseSpace
  };


}//namespace elSpectro
