#include "ProductionProcess.h"
#include "Manager.h"

namespace elSpectro{

   //////////////////////////////////////////////////////////////////
  ProductionProcess::ProductionProcess(const CollidingParticle& p1,const CollidingParticle& p2,decaymodel_ptr  model,decayer_ptr  decayer):
    DecayingParticle{0,decayer,model}, //will set decayer in Init of derived class
    _in1{p1},
    _in2{p2}
  {
   }
  
 // void ProductionProcess::Init() {
 //   auto decayVertexID=Manager::Instance().AddVertex(&DecayVertexPosition());
 //    //Set my vertex position to my "decay vertex"
 //    SetVertex(decayVertexID,&DecayVertexPosition());
 //    //decay vertex position for colliding
 //     // // if(_in1){
 //    _in1.SetVertex(VertexID(),VertexPosition());
 //    // }
 //    //if(_in2){
 //    _in2.SetVertex(VertexID(),VertexPosition());
 //    //}
    
  void ProductionProcess::InitVertex() {
    //clear previous event
    auto* info = GetReactionInfo();
    info->ClearVertices();
    
    // auto decayVertexID=generator.AddVertex(&DecayVertexPosition());
    SetDecayVertexID(info->AddVertex());
 
    //Set my vertex position to my "decay vertex"
    SetVertexID(DecayVertexID());
    //decay vertex position for colliding
    // // if(_in1){
    _in1.SetVertexID(VertexID());
    // }
    //if(_in2){
    _in2.SetVertexID(VertexID());
    //}
    //    std::cout<< "ProductionProcess::InitVertex" <<decayVertexID<<" "<<&_in1<<" "<<&_in2<<" "<<InitialParticles()[0]<<" this "<<this<<std::endl;

    //decay vertex position for products
    auto& products=Model()->Products();
    for(auto* prod: products){
       prod->SetVertexID(VertexID());
    }
  }
  
  void ProductionProcess::PostInit(ReactionInfo* info) {
    
     // if(_in1){
    _in1.PostInit(info);
    //}
    // if(_in2){
    _in2.PostInit(info);
    //}
    //if(_in1&&_in2){
      //Note for now take maximum from sum of both incident particles
      //this ignores that one might be a quasifree nucleon
      //but kinematically at high intiial nucleon momentum
      //we can approach this maximum....
      //info->_Wmax=(*_in1.GetInteracting4Vector()+*_in2.GetInteracting4Vector()).M();
      std::cout<<"ProductionProcess::PostInit maximum W = "<< info->Wmax() <<_in1.GetInteracting4Vector().M()<<" "<<_in2.GetInteracting4Vector().M()<<std::endl;
      std::cout<<"ProductionProcess::PostInit "<<dynamic_cast<ReactionElectroProd*>(info) <<std::endl; 
      //}
      //info->_process=this;
    DecayingParticle::PostInit(info);
  
  }

  //////////////////////////////////////////////////////////////////////////
  ///Loop over all channels and integrate cross section
  ///Use cross sections set branching ratio to each 2-body process
  double ProductionProcess::IntegrateCrossSections(){
    auto Nchannels = Product().Channels().N();
    std::vector<double> branchRatios;
    double total=0.0;
    for(uint i=0;i<Nchannels;++i){
      std::cout<<"ProductionProcess::IntegrateCrossSections, Integrate cross section for channel "<<i<<" out of "<<Nchannels<<std::endl;
      Product().Channels().SetCurrChannel(i);
      auto production_model = dynamic_cast<TwoBodyProduction*>(Product().Model());
      if(production_model==nullptr){
	std::cerr<< "  ProductionProcess::IntegrateCrossSections() model not TwoBodyProduction, invalid "<< std::endl;
	exit(0);
      }
      auto xs = IntegrateCrossSection(production_model);
      branchRatios.push_back(xs);
      total+=xs;
    }
    Product().Channels().SetCurrChannel(0);
    Product().Channels().SetBranchRatios(branchRatios);
    return total;
  }
 
}
