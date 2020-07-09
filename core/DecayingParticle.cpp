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
  double DecayingParticle::GenerateProducts(const CurrentEventInfo* parentInfo){
    bool decayed=false;

    double _maxWeight=1;

    //in case I have any useful info to pass to decay products
    //e.g. polarisations, SDMEs, moments,...
    const CurrentEventInfo* myInfo ={nullptr}; 

    Manager::Instance().FindMassPhaseSpace(Mass(),Model());

    
    while(decayed==false){
      //generate decay product vectors
      //samplingWeight = 1 for phase space decay
      //for others it allows to weith phase space back in 
      auto samplingWeight=Decay();
      if(samplingWeight==0) continue;

      
      //evaluate the model intensity for the product vectors
      if(Model()!=nullptr) myInfo = Model()->Intensity(parentInfo);

      //if event info use its weight, if not assume phse space model = 1.
      samplingWeight*=(myInfo==nullptr ? 1. :  myInfo->_weight);

      // if(!samplingWeight )std::cout<<"DecayingParticle::GenerateProducts " <<Model()->GetName()<<" "<<Class_Name()<<" "<<samplingWeight <<" "<<myInfo->_weight<<std::endl;
      
      //accept/reject this decay
      decayed = samplingWeight > gRandom->Uniform()*_maxWeight ;
      
    }
    auto& unproducts=_decay->UnstableProducts();
    for(auto* prod: unproducts){
      prod->GenerateProducts(myInfo);
    }

    return 1.;
    
  }
 //////////////////////////////////////////////////////////////////////
  void DecayingParticle::Print() const {
    Particle::Print();
    if(Model()) Model()->Print();
  }


}
