//////////////////////////////////////////////////////////////
///
///Class:		LundWriter
///Description:
///             Instance of Writer for Lund format

#pragma once

#include "Writer.h"
#include <string>
#include <fstream>
#include <sstream>

namespace elSpectro{

  class LundWriter : public Writer {

     
     
     LundWriter()=default;
     //don't want default contructor accessible
     //only declaring default constructor
     //so other 5 constructors also defaulted(rule of 5)

   public:
     LundWriter(const std::string& filename);
     ~LundWriter() final;
     LundWriter(const LundWriter& other); //need the virtual destructor...so rule of 5
     LundWriter(LundWriter&&)=default;
     LundWriter& operator=(const LundWriter& other);
     LundWriter& operator=(LundWriter&& other) = default;
     
     void WriteHeader() final{};
     void FillAnEvent() final;
     void Write() final;
     void End() final;

   private:
     //streaming functions
     
     void StreamEventInfo(){
       // Header (Event Info):
       // # of Particles, # of Target Nucleons, # of Target Protons,
       // Pol. of Target, Pol. of Electron, x, y, W2, Q2, nu

       _stream<< "\t "<<_finalParticles.size()<<" "<<1<<" "<<1<<" "<<0.<<" "<<0.
	      <<" "<<0.<<" "<<0.<<" "<<0.<<" "<<0.<<" "<<0.<<"\n";
     }
     void StreamParticle(const Particle* p,int status){
       auto p4=p->P4();
       auto ver=p->VertexPosition();

       //second entry 0. == charge Could add to Particle.h
       _stream<<_id++<<" "<<0.<<" "<<status
	      <<" "<<p->Pdg()<<" "<<0<<" "<<0<<" "
	      <<p4.X()<<" "<<p4.Y()<<" "<<p4.Z()<<" "<<p4.T()<<" "
	      <<p4.M()<<" "<<ver.X()<<" "<<ver.Y()<<" "<<ver.Z()<<"\n";
     }
   
     //data members
     std::ofstream _file; //! output file
     //std::ostream _stream; //! output stream
     std::stringstream _stream;

     long _nEvent={0};
   
     int _id=1;
     
     ClassDef(elSpectro::LundWriter,1); //class Writer
   };


}
