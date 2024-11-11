#include "NBodyPhaseSpace.h"

namespace elSpectro{


    void SampleNBodyPhaseSpace(double parentM, DecayModel* model){
 
      //Note PhaseSpaceWeightMaxFromEquDist does not quite get to max, so increase by 10% to be safe....will give warning if event weight is above this....
      std::vector<double> final_masses;
      model->GetStableMasses(final_masses);
      double max= kine::PhaseSpaceWeightMaxFromEquDist(parentM,final_masses)*1.2;

      double wee=0;

      //calculation of Phase space for model
      //Model PhaseSpaceWeightSq includes all
      //subsequent decays to final particles
      auto calcPhaseSpaceWeight = [model](double M){
         double result=TMath::Sqrt(model->PhaseSpaceWeightSq(M));
	 //_weightCalcN++;
	 return result;
      };
      
      //accept or reject mass combinations until got one
      //as W dependence accounted for elsewere
      while( (wee=calcPhaseSpaceWeight(parentM)) < gRandom->Uniform()*max ){}
      
      if(wee>max){
	// _sampledMax=wee;
	std::cerr<<"MassPhaseSpace check weight >  max,  W "<<parentM<<" this "<<wee<<" "<<max<<" normal max "<<kine::PhaseSpaceWeightMax(parentM,final_masses)<<std::endl;
	std::cout<<" parent "<<parentM<<" ";
	for(auto& m:final_masses)
	  std::cout<<m<<" ";
	std::cout<<std::endl;
      }
 
      // _successN++;
     }
    
}//namespace elSpectro
