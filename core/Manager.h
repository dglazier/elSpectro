//////////////////////////////////////////////////////////////
///
///Class:		Manager
///Description:
///            Class to manage manangers!
///           1) Access ParticleManager via Manager::Instance()->Particles()
///          // 2) Access DecayManager via Manager::Instance()->Decays()
///           3) Access ProductionProcess via Manager::Instance()->Process()
#pragma once

//#include "DecayManager.h"
#include "ProductionProcess.h"
#include "Writer.h"
#include <TMath.h>
#include <TRandom3.h>
#include <TBenchmark.h>

namespace elSpectro{

  class Manager{

  public:
    
    
     void SetWriter(Writer* wr){
       _writer.reset(wr);
     }
     Writer* GetWriter()const {return _writer.get();}
     
     void Write(){
       if(_writer.get()==nullptr)return;
       //Check for bad 4-vectors and miss event
       if(CheckNans(_process->FinalParticles())==false)return;
       _writer->InitEvent(_process->InitialParticles(),_process->FinalParticles(),Reaction()->GetReactionInfo()->GetVertices());
       _writer->FillAnEvent();
       _writer->Write();
     }
    bool CheckNans(const particle_ptrs& sptrs){
      for(auto& p:sptrs){
	if(TMath::IsNaN(p->P4().E() )){
	  std::cout<<"Manager::CheckNans bad event " << _nEventsDone<< std::endl;
	  return false;
	}
      }
      return true;
    }
    
     void CountEvent(){_nEventsDone++;}
  
     bool Finished(){
       if(_nEventsDone==_nEventsToGen)
	 return true;
       return false;
     }
     
     double IntegratedXSection()const {return _integralXSection;}
     void SetNEvents(double n){_nEventsToGen=n;}
     long long GetNEvents()const noexcept{return _nEventsToGen;}
     long long GetNDone()const noexcept{return _nEventsDone;}
    
     void SetNEvents_via_LuminosityTime(double n_or_lum, double beamtime){
       if(beamtime==0){
	 SetNEvents(n_or_lum);
	 return;
       }
       _integralXSection=Reaction()->IntegrateCrossSections();
       _nEventsToGen=n_or_lum*1E-33*beamtime*_integralXSection*Reaction()->BranchingFraction();//1E-33(cm2tonb)
       std::cout<<"Manager::SetNEvents_via_LuminosityTime , going to generate "<<_nEventsToGen<<" events"<<std::endl;
       std::cout<<"\t based on an integrated cross section of "<<_integralXSection<<"; luminosity = "<<n_or_lum<<"; and beamtime of "<<beamtime <<" s "<<std::endl;
     }
     void SetNEvents_via_LuminosityTimeFast(double n_or_lum, double beamtime){
       //Note currently same as non Fast version!!!
       if(beamtime==0){
	 SetNEvents(n_or_lum);
	 return;
       }
       _integralXSection=Reaction()->IntegrateCrossSections();
       _nEventsToGen=n_or_lum*1E-33*beamtime*_integralXSection*Reaction()->BranchingFraction();//1E-33(cm2tonb)
       std::cout<<"Manager::SetNEvents_via_LuminosityTimeFast , going to generate "<<_nEventsToGen<<" events"<<std::endl;
       std::cout<<"\t based on an integrated cross section of "<<_integralXSection<<"; luminosity = "<<n_or_lum<<"; and beamtime of "<<beamtime <<" s "<<std::endl;
     }
     
     void Reaction(ProductionProcess* prod){
       _process.reset(prod);
       //_process->InitVertex(*this);
       prod=nullptr;
     }
    void BoostToLab(LorentzVector& boostme){
      _process->BoostToLab(boostme);
    }
    ProductionProcess* Reaction(){return _process.get();}

     void SetSeed(ULong_t seed = 0){gRandom->SetSeed(seed);}


    //void SetModelForMassPhaseSpace(DecayModel* amodel){_massPhaseSpace.SetModel(amodel);}
    //void SuppressPhaseSpace(double val){_massPhaseSpace.SuppressPhaseSpace(val);}
    // void  FindMassPhaseSpace(double parentM,const  DecayModel* amodel) {
    //   _massPhaseSpace.Find(parentM,amodel);
    // }
    // bool  AcceptPhaseSpace(double parentM) {
    //   return _massPhaseSpace.AcceptPhaseSpace(parentM);
    // }

    void InitGeneration(){
      _process->InitGen();
      if( _writer.get() )_writer->Init(_process->InitialParticles());
    }
    
    // int AddVertex(const LorentzVector* v){
     //   _vertices.push_back(v);
     //   return (_vertices.size()-1);
     // }
     // std::vector<const LorentzVector*>& GetVertices(){return _vertices;}
     
     void Clear(){
     
     }

     void Summary(){
       //_process->Print();
       //      _massPhaseSpace.Print();
      std::cout<<"Integrated Total Cross Section (nb) = "<<IntegratedXSection()<<std::endl;
      }

    void nextEvent(){
      Clear();
      Reaction()->GenerateProducts(nullptr);
      Write();
    }
    void processAll(){
      gBenchmark->Start("generator");//timer  
      while(Finished()==false){
       	//std::cout<<"Manager TEST event number "<<GetNDone()<<std::endl;
	nextEvent();
	CountEvent();
	if(GetNDone()%1000==0) std::cout<<"event number "<<GetNDone()<<std::endl;
      }
      gBenchmark->Stop("generator");//timer  
      gBenchmark->Print("generator");//timer  
    }
      
  private:

    // ParticleManager _particles;
    // DecayManager _decays;

    std::unique_ptr<ProductionProcess> _process;
    std::unique_ptr<Writer> _writer;

    //std::vector<const LorentzVector*> _vertices;
       
    //MassPhaseSpace _massPhaseSpace;


    double _integralXSection={0};
    long long _nEventsToGen={0};
    long long _nEventsDone={0};
    
    ClassDef(elSpectro::Manager,1); //class Manager
  };

}
