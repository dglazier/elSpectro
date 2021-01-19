//////////////////////////////////////////////////////////////
///
///Class:		Manager
///Description:
///            Class to manage manangers!
///           1) Access ParticleManager via Manager::Instance()->Particles()
///           2) Access DecayManager via Manager::Instance()->Decays()
///           3) Access ProductionProcess via Manager::Instance()->Process()
#pragma once

#include "ParticleManager.h"
#include "DecayManager.h"
#include "ProductionProcess.h"
#include "Writer.h"
#include "MassPhaseSpace.h"
#include <TRandom3.h>

namespace elSpectro{

  class Manager{

  public:
    
     static Manager& Instance() { static Manager instance; return instance; }

     ParticleManager& Particles() noexcept{return _particles;}
     DecayManager& Decays() noexcept{return _decays;}
     
     void SetWriter(Writer* wr){
       _writer.reset(wr);
     }
     Writer* GetWriter()const {return _writer.get();}
     
     void Write(){
       if(_writer.get()==nullptr)return;
       _writer->FillAnEvent();
       _writer->Write();
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
    
     void SetNEvents_via_LuminosityTime(double lum, double beamtime){
        _integralXSection=Reaction()->IntegrateCrossSection();
	_nEventsToGen=lum*1E-33*beamtime*_integralXSection*Reaction()->BranchingFraction();//1E-33(cm2tonb)
	std::cout<<"Manager::SetNEvents_via_LuminosityTime , going to generate "<<_nEventsToGen<<" events"<<std::endl;
	std::cout<<"\t based on an integrated cross section of "<<_integralXSection<<"; luminosity = "<<lum<<"; and beamtime of "<<beamtime <<" s "<<std::endl;
     }
     
     void Reaction(ProductionProcess* prod){
       _process.reset(prod);
     }
     
     ProductionProcess* Reaction(){return _process.get();}

     void SetSeed(ULong_t seed = 0){gRandom->SetSeed(seed);}


     void SetModelForMassPhaseSpace(DecayModel* amodel){_massPhaseSpace.SetModel(amodel);}
     void  FindMassPhaseSpace(double parentM,const  DecayModel* amodel) {
       _massPhaseSpace.Find(parentM,amodel);
     }
     bool  AcceptPhaseSpace(double parentM) {
       return _massPhaseSpace.AcceptPhaseSpace(parentM);
     }

     void InitGeneration(){
       _process->InitGen();
       if( _writer.get() )_writer->Init();
     }

     int AddVertex(const LorentzVector* v){
       _vertices.push_back(v);
       return (_vertices.size()-1);
     }
     std::vector<const LorentzVector*>& GetVertices(){return _vertices;}
     
     void Clear(){
     
     }

     void Summary(){
       _process->Print();
      _massPhaseSpace.Print();
      std::cout<<"Integrated Total Cross Section (nb) = "<<IntegratedXSection()<<std::endl;
      }
  private:

    ParticleManager _particles;
    DecayManager _decays;

    std::unique_ptr<ProductionProcess> _process;
    std::unique_ptr<Writer> _writer;

    std::vector<const LorentzVector*> _vertices;
    
    MassPhaseSpace _massPhaseSpace;


    double _integralXSection={0};
    long long _nEventsToGen={0};
    long long _nEventsDone={0};
    
    ClassDef(elSpectro::Manager,1); //class Manager
  };

}
