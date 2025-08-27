#include "NBodyPhaseSpace.h"

namespace elSpectro{


    void SampleNBodyPhaseSpace(double parentM, DecayModel* model){
       //Note PhaseSpaceWeightMaxFromEquDist does not quite get to max, so increase by 10% to be safe....will give warning if event weight is above this....
      std::vector<double> final_masses;
      model->GetStableMasses(final_masses);
      double max= kine::PhaseSpaceWeightMaxFromEquDist(parentM,final_masses)*1.2;
      //double max= kine::PhaseSpaceWeightMax(parentM,final_masses)*1.02;
      // std::cout<<"SampleNBodyPhaseSpace max "<<max<<std::endl;
      // if(parentM>50)max*=0.000001;

      double wee=0;

      //calculation of Phase space for model
      //Model PhaseSpaceWeightSq includes all
      //subsequent decays to final particles
      // Long64_t weightCalcN=0;
      // Double_t maxSample=0.;
      bool resample=false;
      auto calcPhaseSpaceWeight = [model,&resample,max](double M){
	//std::cout<<"SampleNBodyPhaseSpace calcPhaseSpaceWeight   "<<M<<std::endl;
	
	double result=TMath::Sqrt(model->PhaseSpaceWeightSq(M,resample));
	 //weightCalcN++;
	//	std::cout<<"SampleNBodyPhaseSpace  "<<M<<" "<<result<< "< "<<max<<std::endl;
	 //if(maxSample<result)maxSample=result;
	 return result;
      };
      
      //accept or reject mass combinations until got one
      //as W dependence accounted for elsewere
      auto tryN=0;
      while( (wee=calcPhaseSpaceWeight(parentM)) < gRandom->Uniform()*max ){
	//we failed, must resample masses
	resample=true;
	//std::cout<<"SampleNBodyPhaseSpace "<< wee<<" < "<<max<<" "<<kine::PhaseSpaceWeightMaxFromEquDist(parentM,final_masses)<<" "<<tryN++<<std::endl;
      }
      
      if(wee>max){
	// _sampledMax=wee;
	std::cerr<<"MassPhaseSpace check weight >  max,  W "<<parentM<<" this "<<wee<<" "<<max<<" normal max "<<kine::PhaseSpaceWeightMax(parentM,final_masses)<<" equidist "<< kine::PhaseSpaceWeightMaxFromEquDist(parentM,final_masses)<<std::endl;
	std::cout<<" parent "<<parentM<<" "<<final_masses.size()<<" ";
	for(auto& m:final_masses)
	  std::cout<<m<<" ";
	std::cout<<std::endl;
	exit(0);
      }
      // std::cout<<"SampleNBodyPhaseSpace nattepmts = "<<weightCalcN<<" "<<max<<" W "<<parentM<<" final "<<wee<<" max sample "<<maxSample<<std::endl;
      //    for(auto m: final_masses)std::cout<<m<<" ";
      //std::cout<<std::endl;

      //if(parentM>50)exit(0);
      // _successN++;
     }
    
}//namespace elSpectro
