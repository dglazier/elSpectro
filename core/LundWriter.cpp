#include "LundWriter.h"
#include <iostream>

namespace elSpectro{

  ///Constructor to create ouput file and intialise data structures
  LundWriter::LundWriter(const std::string &filename,long evPerFile):
    _filename(filename),
    _file(filename),
    _eventsPerFile(evPerFile)
  {
    if(!_file.is_open()){
      std::cerr<<"LundWriter::LundWriter file "<<filename<<" cannot be opened, exiting..."<<std::endl;
      exit(0);
    }
 
  }
  
  LundWriter::~LundWriter(){
    End();
  }
  ///////////////////////////////////////////////////////////////
  ///Close the file stream
  void LundWriter::End(){
    if(!_file.is_open()) return;
    _file.close();

  }
  ///////////////////////////////////////////////////////////////
  ///Reached max events for this file start another
  void LundWriter::NewFile(){
    End();
    auto newfilename=TString(_filename);
    newfilename.ReplaceAll(".dat",Form("__%d.dat",_nFile++));
    _file.open(newfilename.Data());
  }

  /////////////////////////////////////////////////////////////
  //write all the info required for this event
  void LundWriter::FillAnEvent(){

    
    ////fill _stream
    StreamEventInfo();
 
    _id=1;//reset particle ID counter
 
    
    //final particles
    int final_status=1;
    for(const auto* p:_finalParticles){
      StreamParticle(p,final_status);
    }
      
    _nEvent++;
    
    if(_nEvent%_eventsPerFile==0)
      NewFile();
  }

  /////////////////////////////////////////////////////////
 
  /////////////////////////////////////////////////////////
  void LundWriter::Write(){
    //write stream
    _file<<_stream.rdbuf();
    //reset stream
    _stream.clear();
  }

 
}
