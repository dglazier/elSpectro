#include "DecayingParticle.h"
#include "TwoBodyFlat.h"
#include "Manager.h"
#include "Interface.h"
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
      double lifetime=TDatabasePDG::Instance()->GetParticle(Pdg())->Lifetime();
      double meanFreePath=lifetime*TMath::C()*1000; //in mm
      if( meanFreePath>0.1 ){ //0.1mm
    	_decayType=DecayType::Detached;
	_decVertexDist = new DistTF1(TF1("MFP","TMath::Exp(-x/[0])",0,25*lifetime));//in mm
	_decVertexDist->GetTF1().SetParameter(0,lifetime);
	_decVertexDist->GetTF1().SetNpx(500);
      }
      else _decayType=DecayType::Attached;
    }
    else _decayType=DecayType::Attached;
   
    //decay vertex position
    auto& products=_decay->Products();
    //  std::cout<<" DecayingParticle::PostInit "<<Pdg()<<" "<<products.size()<<std::endl;
    //    if(IsDecay()==DecayType::Detached||IsDecay()==DecayType::Production){
    if(IsDecay()==DecayType::Detached){
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

    if(Pdg()!=0&&Pdg()!=-2211){//for real particles
      Double_t productMasses = 0.0;
      for(auto* prod: products){
	productMasses+=prod->MinimumMassPossible();
      }
      auto mmp=MaximumMassPossible();
      if(mmp<productMasses){
	std::cerr<<"DecayingParticle::PostInit  insufficient mass to decay to its products, max mass = "<<mmp<<" while product masses "<<productMasses<<" for "<<Pdg()<<std::endl;
	Print();
	std::cerr<<"DecayingParticle::PostInit  EXITING "<<std::endl;
	exit(0);
      }
    }

 
    

     //std::cout<<"DecayingParticle::PostInit pdg "<<Pdg()<<" vertexID "<<_decayVertexID<<std::endl;
    // std::cout<<"DecayingParticle::PostInit  min mass "<<MinimumMassPossible()<<std::endl;
  };
  //////////////////////////////////////////////////////////////////////
  DecayStatus   DecayingParticle::GenerateProducts(){

    _generateCalls++;
    
    if(Model()->CheckThreshold()==false) return DecayStatus::ReGenerate;
    
    bool decayed=false;

    double _maxWeight=1;
  
    // std::cout<<"DecayingParticle::GenerateProducts "<<Pdg()<<" "<<Mass()<<" "<<P4().M()<<" "<<_decay->Products().size()<<" "<<" "<<_decay->Products()[0]->Pdg()<<" "<<_decay->UnstableProducts().size()<<" "<<_decay->StableProducts().size()<<" "<<std::endl;
    //if in charge of phase space calculate masses for full decay chain
    Manager::Instance().FindMassPhaseSpace(Mass(),Model());
  
    //generate decay product vectors
    //samplingWeight = 1 for phase space decay
    //for others it allows to weigth phase space back in 

    auto samplingWeight= Decay();

    if(Model()->HasAngularDistribution()==false)samplingWeight=1; //Model has no angular distribution
     
    //samplingWeight ==0 => not physical (below threshold)
    if(samplingWeight==0) return DecayStatus::ReGenerate;
  
    //evaluate the model intensity for the product vectors
    double weight = 1;
    if(Model()!=nullptr)  weight = Model()->Intensity();
 
    if(weight==0)  return DecayStatus::ReGenerate;
    if(samplingWeight - weight < -1E-4 ){//tolerance 0.0001
      std::cout<<"DecayingParticle::GenerateProducts model weight is greater than envelope " <<Mass()<<" "<<Model()->GetName()<<" "<<Class_Name()<<" weights "<<samplingWeight <<" "<<weight<<" masses "<<Model()->Products()[0]->Mass()<<" "<<Model()->Products()[1]->Mass()<<" difference in weights "<<samplingWeight-weight <<std::endl;
      // exit(0);
    }
    //if event info use its weight, if not assume phse space model = 1.
    weight/=samplingWeight;
  
   
 
    //accept/reject this decay
    //if decay depends on variable chosen by parent need to regenerate on fail
    //if decay indendent of parent variables can just try for another
    //   std::cout<<"DecayingParticle::GenerateProducts "<<Pdg()<<" "<<weight <<" "<<_maxWeight<<" "<<samplingWeight<<std::endl;
    decayed = weight > gRandom->Uniform()*_maxWeight ;
    if (decayed == false && (Model()->RegenerateOnFail()==false) )
      return DecayStatus::TryAnother;
    else if (decayed == false && (Model()->RegenerateOnFail()==true) )
      return DecayStatus::ReGenerate;

    //else true

    //decay vertex position
    GenerateVertexPosition();

 
    //std::cout<<"+++++++++++++++++DecayingParticle "<<decayed<<std::endl;
    auto& unproducts=_decay->UnstableProducts();
    for(auto* prod: unproducts){
      DecayStatus prodStatus=DecayStatus::ReGenerate;
      // std::cout<<"+++++DecayingParticle start checking unstable "<<unproducts.size()<<" "<<_decay->StableProducts().size()<<" "<<_decay->Products().size()<<" "<<Pdg()<<std::endl;
  
      while((prodStatus=prod->GenerateProducts()) != DecayStatus::Decayed){

	if(prodStatus==DecayStatus::ReGenerate) return DecayStatus::ReGenerate;
      }
    }

  
    return DecayStatus::Decayed;
    
  }
 //////////////////////////////////////////////////////////////////////
  void DecayingParticle::Print() const {
    Particle::Print();
    std::cout<<"\t DecayParticle GenerateProducts calls "<<_generateCalls<<std::endl;
    if(Model()) Model()->Print();
    
  }

   void DecayingParticle::GenerateVertexPosition()  noexcept{
      //auto old=VertexPosition();
      if( IsDecay()==DecayType::Detached){
	Double_t t0=_decVertexDist->SampleSingle();//in s
	//Need lab 4-vector
	LorentzVector lab=P4();
	generator().BoostToLab(lab);
	
	Double_t r= t0 * lab.Gamma() * TMath::C() * lab.Beta() *1000; //Lorentz contraction , mm
	Double_t labP=lab.P();
	//Set in direction of particle momentum
	//with length of decay
	_decayVertex.SetXYZT(lab.X()/labP*r,lab.Y()/labP*r,lab.Z()/labP*r,r/1000/TMath::C());
	//add production vertex
	_decayVertex+=*VertexPosition(); 
      }
    }
 
}
