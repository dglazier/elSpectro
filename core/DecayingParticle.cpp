#include "DecayingParticle.h"
#include "TwoBodyFlat.h"
#include "Manager.h"
#include <TRandom.h>
#include <TDatabasePDG.h>

namespace elSpectro{

  DecayingParticle::DecayingParticle(DecayModel* model):
    Particle{0},_decay{model},_decayer{new TwoBodyFlat}{
      _decay->SetParent(this);
    }
  DecayingParticle::DecayingParticle(int pdg,DecayModel* model,DecayVectors* decayer):
    Particle{pdg},_decay{model},_decayer{decayer}{
      _decay->SetParent(this);
      
  }
  DecayingParticle::DecayingParticle(int pdg,DecayVectors* decayer,DecayModel* model):
    Particle{pdg},_decayer{decayer},_decay{model}{
      _decay->SetParent(this);

  }
  
//////////////////////////////////////////////////////////////////////
  void DecayingParticle::PostInit(ReactionInfo* info) {
    //Decay type depends on Lifetime
    if(TDatabasePDG::Instance()->GetParticle(Pdg())){
      double meanFreePath=TDatabasePDG::Instance()->GetParticle(Pdg())->Lifetime();
      meanFreePath*=TMath::C()*100; //in mm
      if( meanFreePath>0.01 ){ //10umm
    	_decayType=DecayType::Detached;
      }
      else _decayType=DecayType::Attached;
    }
    else _decayType=DecayType::Attached;
   
   
    //decay vertex position
    auto& products=_decay->Products();
    if(IsDecay()==DecayType::Detached||IsDecay()==DecayType::Production){
      //create a new detached vertex
      _decayVertexID=Manager::Instance().AddVertex(&_decayVertex);
      for(auto* prod: products){
	prod->SetVertex(_decayVertexID,&_decayVertex);
      }
    }
    else{//same vertex as parent
      for(auto* prod: products){
	prod->SetVertex(VertexID(),VertexPosition());
      }
    }

    if(_decay)_decay->PostInit(info);
    if(_decayer)_decayer->PostInit(info);
 
    std::cout<<"DecayingParticle::PostInit pdg "<<Pdg()<<" vertexID "<<_decayVertexID<<std::endl;
    std::cout<<"DecayingParticle::PostInit  min mass "<<MinimumMassPossible()<<std::endl;
  };
  //////////////////////////////////////////////////////////////////////
  DecayStatus   DecayingParticle::GenerateProducts(){

    
    if(Model()->CheckThreshold()==false) return DecayStatus::ReGenerate;
    
    bool decayed=false;

    double _maxWeight=1;

  
    // std::cout<<"DecayingParticle::GenerateProducts "<<Pdg()<<" "<<Mass()<<" "<<P4().M()<<" "<<_decay->Products().size()<<" "<<" "<<_decay->Products()[0]->Pdg()<<" "<<_decay->Products()[1]->Pdg()<<std::endl;
    //if in charge of phase space calculate masses for full decay chain
    Manager::Instance().FindMassPhaseSpace(Mass(),Model());
  
    //generate decay product vectors
    //samplingWeight = 1 for phase space decay
    //for others it allows to weigth phase space back in 
    auto samplingWeight=Decay();
    
    //samplingWeight ==0 => not physical (below threshold)
    if(samplingWeight==0) return DecayStatus::ReGenerate;
  
    //evaluate the model intensity for the product vectors
    double weight = 1;
    if(Model()!=nullptr)  weight = Model()->Intensity();
    if(weight==0)  return DecayStatus::ReGenerate;
    if(samplingWeight<weight){
      std::cout<<"DecayingParticle::GenerateProducts model weight is greater than envelope " <<Mass()<<" "<<Model()->GetName()<<" "<<Class_Name()<<" weights "<<samplingWeight <<" "<<weight<<" masses "<<Model()->Products()[0]->Mass()<<" "<<Model()->Products()[1]->Mass()<<std::endl;
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

    //decay vertex position
    GenerateVertexPosition();

 
    // std::cout<<"+++++++++++++++++DecayingParticle "<<decayed<<std::endl;
    auto& unproducts=_decay->UnstableProducts();
    for(auto* prod: unproducts){
      DecayStatus prodStatus=DecayStatus::ReGenerate;
      while((prodStatus=prod->GenerateProducts()) != DecayStatus::Decayed){
	if(prodStatus==DecayStatus::ReGenerate) return DecayStatus::ReGenerate;
      }
    }

    /*
    //stable products need a vertex too
    auto& stproducts=_decay->StableProducts();
    // std::cout<<"DecayingPArticle size "<<stproducts.size()<<" "<<Pdg()<<" "<<_decayVertexID<<" "<<Model()->GetName()<<std::endl;
    for(auto* prod: stproducts){
      prod->SetVertex(_decayVertexID,_decayVertex);
    }
    */
    
    return DecayStatus::Decayed;
    
  }
 //////////////////////////////////////////////////////////////////////
  void DecayingParticle::Print() const {
    Particle::Print();
    if(Model()) Model()->Print();
  }


}
