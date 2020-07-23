#include "DecayingParticle.h"
#include "TwoBodyFlat.h"
#include "Manager.h"
#include <TRandom.h>

namespace elSpectro{

  DecayingParticle::DecayingParticle(DecayModel* model):
    Particle{0},_decay{model},_decayer{new TwoBodyFlat}{

  }
  DecayingParticle::DecayingParticle(int pdg,DecayModel* model,DecayVectors* decayer):
    Particle{pdg},_decay{model},_decayer{decayer}{

  }
  DecayingParticle::DecayingParticle(int pdg,DecayVectors* decayer,DecayModel* model):
    Particle{pdg},_decayer{decayer},_decay{model}{

  }
  //////////////////////////////////////////////////////////////////////
  DecayStatus   DecayingParticle::GenerateProducts(){

    //std::cout<<"DecayingParticle::GenerateProducts "<<Model()->MinimumMassPossible()<<" "<<Mass()<<" p4 "<<P4().M()<<std::endl;
    if(Model()->CheckThreshold()==false) return DecayStatus::ReGenerate;
    
    bool decayed=false;

    double _maxWeight=1;

  
    // std::cout<<"DecayingParticle::GenerateProducts "<<Pdg()<<std::endl;
    //if in charge of phase space calculate masses for full decay chain
    Manager::Instance().FindMassPhaseSpace(Mass(),Model());
    
    // while(decayed==false){

      //generate decay product vectors
      //samplingWeight = 1 for phase space decay
      //for others it allows to weigth phase space back in 
    auto samplingWeight=Decay();

    
    //std::cout<<MassWeight()<<" "<<samplingWeight<<std::endl;
    //    /MassWeight();
    //auto samplingWeight=Decay();
      
  
    //samplingWeight ==0 => not physical (below threshold)
    if(samplingWeight==0) return DecayStatus::ReGenerate;
    
    // std::cout<<" DecayPArticle W = "<<Mass()<<" "<<P4()<<std::endl; 
    //evaluate the model intensity for the product vectors
    double weight = 1;
    if(Model()!=nullptr)  weight = Model()->Intensity();
    if(weight==0)  return DecayStatus::ReGenerate;
    if(samplingWeight<weight){
      std::cout<<"DecayingParticle::GenerateProducts model weight is greater than envelope " <<Mass()<<" "<<Model()->GetName()<<" "<<Class_Name()<<" "<<samplingWeight <<" "<<weight<<" "<<Model()->Products()[0]->Mass()<<" "<<Model()->Products()[1]->Mass()<<std::endl;
    //exit(0);
    }
    //if event info use its weight, if not assume phse space model = 1.
    weight/=samplingWeight;
    
 
    //accept/reject this decay
    //if decay depends on variable chosen by parent need to regenerate on fail
    //if decay indendent of parent variables can just try for another
    decayed = weight > gRandom->Uniform()*_maxWeight ;
    if (decayed == false && (Model()->RegenerateOnFail()==false) )
      return DecayStatus::TryAnother;
    else if (decayed == false && (Model()->RegenerateOnFail()==true) )
      return DecayStatus::ReGenerate;

    //else true
    
   
    // std::cout<<"+++++++++++++++++DecayingParticle "<<decayed<<std::endl;
    auto& unproducts=_decay->UnstableProducts();
    for(auto* prod: unproducts){
      DecayStatus prodStatus=DecayStatus::ReGenerate;
      while((prodStatus=prod->GenerateProducts()) != DecayStatus::Decayed){
	if(prodStatus==DecayStatus::ReGenerate) return DecayStatus::ReGenerate;
	//else ==DecayStatus::TryAnother
      }
    }
    
    return DecayStatus::Decayed;
    
  }
 //////////////////////////////////////////////////////////////////////
  void DecayingParticle::Print() const {
    Particle::Print();
    if(Model()) Model()->Print();
  }


}
