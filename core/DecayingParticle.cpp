#include "DecayingParticle.h"
#include "TwoBodyFlat.h"
#include "TwoBodyProduction.h"
#include "Manager.h"
#include "ProductionProcess.h"
#include <TRandom.h>
#include <TDatabasePDG.h>

namespace elSpectro{

//////////////////////////////////////////////////////////////////////
  DecayingParticle::DecayingParticle(decaymodel_ptr model):
    Particle{0}{
    _channels.AddDecay(this,1.,std::move(model),CloneDecayer(TwoBodyFlat()));
     }
//////////////////////////////////////////////////////////////////////
  DecayingParticle::DecayingParticle(int pdg,decaymodel_ptr model,decayer_ptr decayer):
    Particle{pdg}{
    _channels.AddDecay(this,1.,std::move(model),std::move(decayer));
       
  }
//////////////////////////////////////////////////////////////////////
  DecayingParticle::DecayingParticle(int pdg,decayer_ptr decayer,decaymodel_ptr model):
    Particle{pdg}{
    _channels.AddDecay(this,1.,std::move(model),std::move(decayer));

  }
  ///////////////////////////////////////////////////////////////////
  void DecayingParticle::CopyOther(const DecayingParticle& other){
    //std::cout<<"DecayingParticle::CopyOther( "<<this <<" other "<<&other <<" "<<other.Pdg()<<" "<<other._channels.CurrModel()->Product(0)->Pdg()<<" "<<other._channels.CurrModel()->Product(0)<<std::endl;
    _channels=other._channels;
    // std::cout<<"DecayingParticle::CopyOther( "<<_channels.N()<<std::endl;
    _decVertexDist=other._decVertexDist;
    
    _minMass=other._minMass;
    _decayVertex=other._decayVertex;
    _decayVertexID=other._decayVertexID;
    _decayType=other._decayType;

    _generateCalls=other._generateCalls;
    _gtSample=other._gtSample;

    //important, my models must have me as their parent
    _channels.SetParent(this);
  
  }
  ////////////////////////////////////////////////////////////////////
  DecayingParticle::DecayingParticle(const DecayingParticle& other):Particle(other){
    CopyOther(other);  
  }
  ////////////////////////////////////////////////////////////////////
  DecayingParticle::DecayingParticle(DecayingParticle&& other):Particle(other){
    CopyOther(other);
  }
  ////////////////////////////////////////////////////////////////////
  DecayingParticle& DecayingParticle::operator=(const DecayingParticle& other){
    Particle::operator=(other);
    CopyOther(other);  
    return *this;
  }
  ////////////////////////////////////////////////////////////////////
  DecayingParticle& DecayingParticle::operator=(DecayingParticle&& other){
    Particle::operator=(other);
    CopyOther(other);
    return *this;
  }
  //////////////////////////////////////////////////////////////////////
  void DecayingParticle::PostInit(ReactionInfo* info) {
    //std::cout<<" DecayingParticle::PostInit "<<Pdg()<<" with n decay channels = "<<_channels.N()<<std::endl;
    if(_channels.N() == 0){
      throw std::runtime_error("DecayingParticle::PostInit, no decay channels given");
    }
    //_process = info->_process;
   
    
    //Decay type depends on Lifetime
    //if(TDatabasePDG::Instance()->GetParticle(Pdg())){
      //double lifetime=TDatabasePDG::Instance()->GetParticle(Pdg())->Lifetime();
      //double meanFreePath=lifetime*TMath::C()*1000; //in mm
      //std::cout<<"DecayingParticle::PostInit "<<Pdg()<<" "<<lifetime<<" "<<meanFreePath<<std::endl;
      //if( meanFreePath>0.1 ){ //0.1mm
      //	_decayType=DecayType::Detached;
	//	_decVertexDist = new DistTF1(TF1("MFP","TMath::Exp(-x/[0])",0,25*lifetime));//in mm
	//_decVertexDist->GetTF1().SetParameter(0,lifetime);
	//_decVertexDist->GetTF1().SetNpx(500);
      //}
      //else _decayType=DecayType::Attached;
      //}
      //else _decayType=DecayType::Attached;
   
    //decay vertex position


    
      /*
      //  std::cout<<" DecayingParticle::PostInit "<<Pdg()<<" "<<products.size()<<std::endl;
    //    if(IsDecay()==DecayType::Detached||IsDecay()==DecayType::Production){
    if(IsDecay()==DecayType::Detached){
      //create a new detached vertex
      //_decayVertexID=Manager::Instance().AddVertex(&_decayVertex);
      _decayVertexID=info->AddVertex();
      // _decayVertexID=0;
      for(auto* prod: products){
	prod->SetVertexID(_decayVertexID);
      }
    }
    else{//same vertex as parent
      for(auto* prod: products){
	prod->SetVertexID(VertexID());
      }
    }

    */

   
    _channels.PostInit(info);

    auto& products=Model()->Products();

  
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
  
  }
  
//////////////////////////////////////////////////////////////////////
  DecayStatus   DecayingParticle::GenerateProducts(const ProductionProcess* production){
    // if(Pdg()==-2211) std::cout<<"********************************DecayingParticle::GenerateProducts "<<Pdg()<<" "<<Mass()<<" "<<P4().M()<<" "<<Model()->Products().size()<<" "<<" "<<Model()->Products()[0]->Pdg()<<" prod mass "<<Model()->Products()[0]->Mass()<<" "<<Model()->UnstableProducts().size()<<" "<<Model()->StableProducts().size()<<" "<<std::endl;
    //std::cout<<"********************************DecayingParticle::GenerateProducts "<<Pdg()<<" "<<P4()<<" mass "<<Mass()<<" "<<P4().M()<<" minmass "<<MinimumMassPossible()<<" ? "<<(Mass()-MinimumMassPossible())<<" "<<Model()->Products().size()<<" "<<" pdg1 "<<Model()->Products()[0]->Pdg()<<" pdg2 "<<Model()->Products()[1]->Pdg()<<" prod mass "<<Model()->Products()[0]->Mass()<<" "<<Model()->UnstableProducts().size()<<" "<<Model()->StableProducts().size()<<" "<<std::endl;

    _generateCalls++;

    if(Mass()<MinimumMassPossible()) return DecayStatus::ReGenerate;
    //std::cout<<"DecayingParticle::GenerateProducts got here 1"<<std::endl;
    if(Model()->ReadyForDecay()==false) return DecayStatus::ReGenerate;
    //std::cout<<"DecayingParticle::GenerateProducts got here 2"<<std::endl;
   
    bool decayed=false;

    double _maxWeight=1;
  
   
    //generate decay product vectors
    //samplingWeight = 1 for phase space decay
    //for others it allows to weigth phase space back in 

    auto samplingWeight= Decay();
    // std::cout<<"DecayingParticle::GenerateProducts "<<Pdg()<<" sample weight "<<samplingWeight<<" "<<Model()->HasAngularDistribution()<<std::endl;
    if(Model()->HasAngularDistribution()==false)samplingWeight=1; //Model has no angular distribution
 

    //samplingWeight ==0 => not physical (below threshold)
    if(samplingWeight==0) return DecayStatus::ReGenerate;
    //if(samplingWeight==0) return DecayStatus::TryAnother;
      //evaluate the model intensity for the product vectors
    double weight = 1.;
    
    // std::cout<<"DecayingParticle::GenerateProducts GetIntensity"<<std::endl;
    if(Model()!=nullptr)  weight = Model()->Intensity();
    //if in charge of phase space calculate masses for full decay chain
    // std::cout<<"DecayingParticle::GenerateProducts GetIntensity "<<weight<<std::endl;

    // if(TMath::IsNaN(weight)||TMath::Abs(weight)==TMath::Infinity()){
    //   std::cout<<"$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$DecayingParticle::GenerateProducts "<<Pdg()<<" "<<Mass()<<" "<<P4().M()<<" "<<Model()->Products().size()<<" W "<<Model()->Products()[0]->Mass()<<" "<<Model()->Products()[0]->Pdg()<<" "<<Model()->UnstableProducts().size()<<" "<<Model()->StableProducts().size()<<" "<<std::endl;     //exit(0);
    // }
    if(weight==0)  return DecayStatus::ReGenerate;
    if(TMath::IsNaN(weight)||TMath::Abs(weight)==TMath::Infinity())  return DecayStatus::ReGenerate;
    //if(TMath::IsNaN(weight))  return DecayStatus::TryAnother;
    
    if(samplingWeight - weight < -1E-4 ){//tolerance 0.0001
      std::cout<<"DecayingParticle::GenerateProducts model weight is greater than envelope W =" <<Mass()<<" "<<Model()->GetName()<<" "<<Class_Name()<<" sampling weight "<<samplingWeight <<" current weight "<<weight<<" masses meson "<<Model()->Products()[0]->Mass()<<" baryon "<<Model()->Products()[1]->Mass()<<" ratio in weights "<<samplingWeight/weight <<std::endl;
      auto twoBody = dynamic_cast<TwoBodyProduction*>(Model());
      std::cout<<"W "<< twoBody->get_W()<<" t "<<twoBody->get_t()<<" th "<<twoBody->get_cosThCM()<<" Q2 "<<twoBody->GetPhoton().M2()<<" phase space correct Q2 "<<twoBody->Q2PhaseSpaceCorrect()<<" mass "<< twoBody->MassPhaseSpaceCorrect()<< std::endl; 
      // exit(0);
      _gtSample++;
    }
    // else{
    //   // std::cout<<"DecayingParticle::GenerateProducts model weight OK"<<std::endl;
    //   if(Pdg()==-2211){
    // 	std::cout<<"DecayingParticle::GenerateProducts  W =" <<Mass()<<" "<<Model()->GetName()<<" "<<Class_Name()<<" sampling weight "<<samplingWeight <<" current weight "<<weight<<" masses meson "<<Model()->Products()[0]->Mass()<<" baryon "<<Model()->Products()[1]->Mass()<<" ratio in weights "<<samplingWeight/weight <<std::endl;
    // 	auto twoBody = dynamic_cast<TwoBodyProduction*>(Model());
    // 	std::cout<<"W "<< twoBody->get_W()<<" t "<<twoBody->get_t()<<" th "<<twoBody->get_cosThCM()<<" "<<kine::tFromcosthW(twoBody->get_cosThCM(), twoBody->get_W(),twoBody->GetPhoton().M(),twoBody->GetTarget().M(),Model()->Product(0)->Mass(),Model()->Product(1)->Mass())<<" p1 "<<Model()->Product(1)->Mass()<<" p0 "<<Model()->Product(0)->Mass()<<" ph "<<twoBody->GetPhoton().M()<<" tar "<<twoBody->GetTarget().M()<<std::endl; 
    //   }
    // }
    // if(Pdg()==-2211)std::cout<<"DecayingParticle::GenerateProducts calc weihgt "<<Pdg()<<" "<<weight <<" "<<_maxWeight<<" "<<samplingWeight<<" "<<Model()->RegenerateOnFail()<<" "<<weight/samplingWeight<<std::endl;
    //if event info use its weight, if not assume phse space model = 1.
    weight/=samplingWeight;
        
 
    //accept/reject this decay
    //if decay depends on variable chosen by parent need to regenerate on fail
    //if decay indendent of parent variables can just try for another
    decayed = weight > gRandom->Uniform() ;
    //std::cout<<"+++++++++++++++++DecayingParticle "<<decayed<<std::endl;
    //if(Pdg()==-2211) exit(0);
    if (decayed == false && (Model()->RegenerateOnFail()==false) )
      return DecayStatus::TryAnother;
    else if (decayed == false && (Model()->RegenerateOnFail()==true) )
      return DecayStatus::ReGenerate;

    //else true

      //decay vertex position
    GenerateVertexPosition(production);

 
    auto& unproducts=Model()->MutableUnstableProducts();
    for(auto& prod: unproducts){
      DecayStatus prodStatus=DecayStatus::ReGenerate;
      // std::cout<<"+++++DecayingParticle start checking unstable "<<unproducts.size()<<" "<<Model()->StableProducts().size()<<" "<<Model()->Products().size()<<" "<<Pdg()<<std::endl;
  
      while((prodStatus=prod.GenerateProducts(production)) != DecayStatus::Decayed){
      	//std::cout<<"DEcayingPArticle prodStatus "<<Pdg()<<" "<<((Int_t)(prodStatus))<<" tryanother "<<((Int_t)(DecayStatus::TryAnother))<<" regen "<<((Int_t)(DecayStatus::ReGenerate))<<" decayed "<<((Int_t)(DecayStatus::Decayed))<<std::endl;
	if(prodStatus==DecayStatus::ReGenerate) return DecayStatus::ReGenerate;
      }
    }

  
    // std::cout<<"+++++++++++++++++DecayingParticle we did decay "<<Pdg()<<" "<<std::endl;
    return DecayStatus::Decayed;
    
  }
 //////////////////////////////////////////////////////////////////////
  void DecayingParticle::Print() const {
    Particle::Print();
    std::cout<<"\t DecayParticle GenerateProducts calls "<<_generateCalls<<" of which "<<_gtSample<<" were greater than sample weight "<< static_cast<float>(_gtSample)/(_generateCalls>0?_generateCalls:1) <<std::endl;
    if(Model()) Model()->Print();
    
  }

   void DecayingParticle::GenerateVertexPosition(const ProductionProcess* production)  noexcept{
      _decayVertex = VertexPosition();
      if( IsDecay()==DecayType::Detached){
	//	std::cout<<"\t DecayParticle GenerateVertexPosition "<<Pdg()<<" "<<" "<<production<<std::endl;
	Double_t t0=_decVertexDist.SampleSingle();//in s
	//Need lab 4-vector
	LorentzVector lab=P4();
	production->BoostToLab(lab);
	//std::cout<<"\t DecayParticle GenerateVertexPosition "<<Pdg()<<std::endl;
	Double_t r= t0 * lab.Gamma() * TMath::C() * lab.Beta() *1000; //Lorentz contraction , mm
	Double_t labP=lab.P();
	//Set in direction of particle momentum
	//with length of decay
	_decayVertex.SetXYZT(
			     _decayVertex.X()+lab.X()/labP*r,_decayVertex.Y()+lab.Y()/labP*r,
			     _decayVertex.Z()+lab.Z()/labP*r,_decayVertex.T()+r/1000/TMath::C());
	
	//add new decay vertex for writing
	_decayVertexID=production->GetReactionInfo()->AddVertex();
	auto products = Model()->MutableProducts();
	for(auto* prod: products){
	  prod->SetVertexID(_decayVertexID);
	}
      }
      //set vertex position of decays
      auto products = Model()->MutableProducts();
      for(auto prod:products){
	prod->SetVertexPosition(_decayVertex);
      }
      // std::cout<<"\t DecayParticle GenerateVertexPosition "<<VertexPosition()<<" "<<_decayVertex<<std::endl;

    }
 
}
