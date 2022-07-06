//////////////////////////////////////////////////////////////
///
///Class:		GlueXWriter
///Description:
///             Instance of Writer for Lund format

#pragma once

#include "Writer.h"
#include <TDatabasePDG.h>
#include <string>
#include <fstream>
#include <sstream>

namespace elSpectro{

  class GlueXWriter : public Writer {

     
     
     GlueXWriter()=default;
     //don't want default contructor accessible
     //only declaring default constructor
     //so other 5 constructors also defaulted(rule of 5)

   public:
     GlueXWriter(const std::string& filename,long evPerFile=1E18, int runnumber=72068);
     ~GlueXWriter() final;
     GlueXWriter(const GlueXWriter& other); //need the virtual destructor...so rule of 5
     GlueXWriter(GlueXWriter&&)=default;
     GlueXWriter& operator=(const GlueXWriter& other);
     GlueXWriter& operator=(GlueXWriter&& other) = default;
     
     void WriteHeader() final{};
     void FillAnEvent() final;
     void Write() final;
     void End() final;
     void Init() final;
     void NewFile();
     
   private:
     //streaming functions
     // output genr8 data format which can then be converted to hddm with genr8_2_hddm
     void StreamEventInfo(){
       // Header (Event Info):
       // run number, event number (start counting at 1), # of Particles
  
       _stream<< _runnumber << " " << _nEvent+1 << " " <<_finalParticles.size()<<"\n";
     }
     void StreamParticle(const Particle* p,int status){
       // index, PID, mass
       // charge, P.x, P.y, P.z, E
       auto p4=p->P4();
       int geantCode = TDatabasePDG::Instance()->ConvertPdgToGeant3(p->Pdg());
       int charge = TDatabasePDG::Instance()->GetParticle(p->Pdg())->Charge()/3;
       _stream<<_id++<<" "<<geantCode<<" "<<p4.M()<<"\n"
	      <<"\t"<< charge <<" "<<p4.X()<<" "<<p4.Y()<<" "<<p4.Z()<<" "<<p4.T()<<"\n";
     }
   
     //data members
     std::ofstream _file; //! output file
     //std::ostream _stream; //! output stream
     std::stringstream _stream;
     std::string _filename;
     
     long _nEvent={0};
     int _nFile={1};
     long _eventsPerFile=static_cast<long>(1E18);
     int _runnumber=72068;
     
     int _id=1;
     int _beamPdg=0;
     int _targetPdg=0;
     
     const Particle* _inBeam={nullptr};
     const Particle* _inTarget={nullptr};
     
     ClassDef(elSpectro::GlueXWriter,1); //class Writer
   };


}

