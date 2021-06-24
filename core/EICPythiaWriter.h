//////////////////////////////////////////////////////////////
///
///Class:		EICPythiaWriter
///Description:
///             Instance of Writer for EICPythia format

#pragma once

#include "Writer.h"
#include <string>
#include <fstream>
#include <sstream>

namespace elSpectro{

  class EICPythiaWriter : public Writer {

     
     
     EICPythiaWriter()=default;
     //don't want default contructor accessible
     //only declaring default constructor
     //so other 5 constructors also defaulted(rule of 5)

   public:
     EICPythiaWriter(const std::string& filename,long evPerFile=1E18);
     ~EICPythiaWriter() final;
     EICPythiaWriter(const EICPythiaWriter& other); //need the virtual destructor...so rule of 5
     EICPythiaWriter(EICPythiaWriter&&)=default;
     EICPythiaWriter& operator=(const EICPythiaWriter& other);
     EICPythiaWriter& operator=(EICPythiaWriter&& other) = default;
     
     void WriteHeader() final{};
     void FillAnEvent() final;
     void Write() final;
     void End() final;
     void Init() final;
     void NewFile();
     
   private:
     //streaming functions
     
     void StreamEventInfo(){
       // Header (Event Info):
       // # of Particles, # of Target Nucleons, # of Target Protons,
       // Pol. of Target, Pol. of Electron,
       // BeamType, BeamEnergy,Target ID, ProcessID, Weight
  
       _stream<<"0"<< "\t"<<_nEvent<< "\t"<<_finalParticles.size()<<"\n";
       _stream << "============================================\n";
    }
     void StreamParticle(const Particle* p,int status){
       auto p4=p->P4();
       auto ver=p->VertexPosition();
       //second entry 0. == lifetime Could add to Particle.h
       //note assume cm
       int parent=0;
       int daughter_first=0;
       int daughter_last=0;
       _stream<<_id++<<"\t"<<status
	      <<"\t"<<p->Pdg()<<"\t"<<parent<<"\t"
	      <<daugther_first<<"\t"<<daughter_last<<"\t"
	      <<p4.X()<<"\t"<<p4.Y()<<"\t"<<p4.Z()<<"\t"<<p4.T()<<"\t"
	      <<p4.M()<<"\t"<<ver->X()<<"\t"<<ver->Y()<<"\t"<<ver->Z()<<"\n";
       _stream << "============================================\n";
     }
   
     //data members
     std::ofstream _file; //! output file
     //std::ostream _stream; //! output stream
     std::stringstream _stream;
     std::string _filename;
     
     long _nEvent={0};
     int _nFile={1};
     long _eventsPerFile=static_cast<long>(1E18);
     
     int _id=1;
     int _beamPdg=0;
     int _targetPdg=0;
     
     const Particle* _inBeam={nullptr};
     const Particle* _inTarget={nullptr};
    Particle _photon={22}; //virtual photon pdg=22
    
     ClassDef(elSpectro::EICPythiaWriter,1); //class Writer
   };


}

