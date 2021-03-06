#include "EICSimpleWriter.h"
#include <TDatabasePDG.h>
#include <iostream>

namespace elSpectro{

  ///Constructor to create ouput file and intialise data structures
  EICSimpleWriter::EICSimpleWriter(const std::string &filename,long evPerFile):
    _filename(filename),
    _file(filename),
    _eventsPerFile(evPerFile)
  {
    if(!_file.is_open()){
      std::cerr<<"EICSimpleWriter::EICSimpleWriter file "<<filename<<" cannot be opened, exiting..."<<std::endl;
      exit(0);
    }
    _stream << "SIMPLE Event FILE"  << std::endl;
    _stream << "============================================" << std::endl;
    _stream << "    I, ievent, nParticles"  << std::endl;
    _stream << "============================================"  << std::endl;
    _stream << "I  K(I,1)  K(I,2)  K(I,3)  K(I,4)  K(I,5)  P(I,1)  P(I,2)  P(I,3)  P(I,4)  P(I,5)  V(I,1)  V(I,2)  V(I,3)"  << std::endl;
    _stream << "============================================"  << std::endl;
    
    Write();

  }
  
  EICSimpleWriter::~EICSimpleWriter(){
    End();
  }
  ///////////////////////////////////////////////////////////////
  ///Get beam and target, or e-/gamma and nucleus
  void EICSimpleWriter::Init(){

    Writer::Init();
    //need to find e- and baryon
    if(TDatabasePDG::Instance()->GetParticle(_initialParticles[0]->Pdg())->ParticleClass()==TString("Baryon") ){
      _inTarget=_initialParticles[0];
      _inBeam=_initialParticles[1]; 
    }
    else {
      _inTarget=_initialParticles[1];
      _inBeam=_initialParticles[0];
    }
 
    _beamPdg=_inBeam->Pdg();
    _targetPdg=_inTarget->Pdg();

    _photon.SetVertex(_inBeam->VertexID(),_inBeam->VertexPosition());
 
  }
  ///////////////////////////////////////////////////////////////
  ///Close the file stream
  void EICSimpleWriter::End(){
    if(!_file.is_open()) return;
    _file.close();

  }
  ///////////////////////////////////////////////////////////////
  ///Reached max events for this file start another
  void EICSimpleWriter::NewFile(){
    End();
    auto newfilename=TString(_filename);
    newfilename.ReplaceAll(".dat",Form("__%d.dat",_nFile++));
    _file.open(newfilename.Data());
  }

  /////////////////////////////////////////////////////////////
  //write all the info required for this event
  void EICSimpleWriter::FillAnEvent(){
     
    ////fill _stream
    StreamEventInfo();
 
    _id=1;//reset particle ID counter
 
    //initial particles
    int initial_status=21;
    StreamParticle(_inBeam,initial_status);
    StreamParticle(_inTarget,initial_status);
    StreamParticle(&_photon,initial_status);
   
    //final particles
    int final_status=1;
    for(const auto* p:_finalParticles){
     
      StreamParticle(p,final_status);
  
     }
    _stream << "=============== Event finished ===============\n";

    _nEvent++;
    
    if(_nEvent%_eventsPerFile==0)
      NewFile();
  }

  /////////////////////////////////////////////////////////
 
  /////////////////////////////////////////////////////////
  void EICSimpleWriter::Write(){
    //write stream
    _file<<_stream.rdbuf();
    //reset stream
    _stream.clear();
  }

 
}
