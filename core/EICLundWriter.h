//////////////////////////////////////////////////////////////
///
///Class:		EICLundWriter
///Description:
///             Instance of Writer for EICLund format

#pragma once

#include "Writer.h"
#include <string>
#include <fstream>
#include <sstream>

namespace elSpectro{

  class EICLundWriter : public Writer {

     
     
     EICLundWriter()=default;
     //don't want default contructor accessible
     //only declaring default constructor
     //so other 5 constructors also defaulted(rule of 5)

   public:
     EICLundWriter(const std::string& filename,long evPerFile=1E18);
     ~EICLundWriter() final;
     EICLundWriter(const EICLundWriter& other); //need the virtual destructor...so rule of 5
     EICLundWriter(EICLundWriter&&)=default;
     EICLundWriter& operator=(const EICLundWriter& other);
     EICLundWriter& operator=(EICLundWriter&& other) = default;
     
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
  
       _stream<< "\t "<<_finalParticles.size()<<" "<<1<<" "<<1
	      <<" "<<0.<<" "<<0.
	      <<" "<<_beamPdg<<" "<<_inBeam->P4().E()<<" "<<_targetPdg<<" "<< _inTarget->P4().E() <<" "<<0.<<"\n";
     }
     void StreamParticle(const Particle* p,int status){
       auto p4=p->P4();
       auto ver=p->VertexPosition();
       //second entry 0. == lifetime Could add to Particle.h
       //note convert vertex mm->cm
       _stream<<_id++<<" "<<0.<<" "<<status
	      <<" "<<p->Pdg()<<" "<<0<<" "<<0<<" "
	      <<p4.X()<<" "<<p4.Y()<<" "<<p4.Z()<<" "<<p4.T()<<" "
	      <<p4.M()<<" "<<ver->X()/10<<" "<<ver->Y()/10<<" "<<ver->Z()/10<<"\n";
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
     
     ClassDef(elSpectro::EICLundWriter,1); //class Writer
   };


}

