#include "ProductionProcess.h"
#include "Manager.h"

namespace elSpectro{

   //////////////////////////////////////////////////////////////////
  ProductionProcess::ProductionProcess(CollidingParticle* p1,CollidingParticle* p2,DecayModel* model):
    DecayingParticle{0,nullptr,model}, //will set decayer in Init of derived class
    _in1{p1},
    _in2{p2}
  {
    Init();
  }
  
   //////////////////////////////////////////////////////////////////
    ProductionProcess::ProductionProcess(DecayModel* model):
    DecayingParticle{model}
  {

    Init();
  }
   //////////////////////////////////////////////////////////////////
  ProductionProcess::ProductionProcess(int pdg,DecayVectors* decayer, DecayModel* model):
    DecayingParticle{pdg,decayer,model}{
    Init();

  }
 void ProductionProcess::Init() {
   auto decayVertexID=Manager::Instance().AddVertex(&DecayVertexPosition());
    //Set my vertex position to my "decay vertex"
    SetVertex(decayVertexID,&DecayVertexPosition());
    //decay vertex position for colliding
    if(_in1){
      _in1->SetVertex(VertexID(),VertexPosition());
    }
    if(_in2){
      _in2->SetVertex(VertexID(),VertexPosition());
    }

    
    //decay vertex position for products
    auto& products=Model()->Products();
    for(auto* prod: products){
      std::cout<<" ProductionProcess::Init() "<<prod->Pdg()<<" "<<VertexPosition()<<std::endl;
	prod->SetVertex(VertexID(),VertexPosition());
    }

 }
  void ProductionProcess::PostInit(ReactionInfo* info) {
    
 
    if(_in1){
      _in1->PostInit(info);
    }
    if(_in2){
      _in2->PostInit(info);
    }
    if(_in1&&_in2){
      //Note for now take maximum from sum of both incident particles
      //this ignores that one might be a quasifree nucleon
      //but kinematically at high intiial nucleon momentum
      //we can approach this maximum....
      info->_Wmax=(*_in1->GetInteracting4Vector()+*_in2->GetInteracting4Vector()).M();
      std::cout<<"ProductionProcess::PostInit maximum W = "<< info->_Wmax <<_in1->GetInteracting4Vector()->M()<<" "<<_in2->GetInteracting4Vector()->M()<<std::endl;
    }
     DecayingParticle::PostInit(info);
  
  }
}
