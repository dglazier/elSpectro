//////////////////////////////////////////////////////////////
///
///Class:		HepMC3Writer
///Description:
///             Instance of Writer for HepMC3 format
///             Requires HepMC3 installed and env variable HEPMC3 set

#pragma once

#include "Writer.h"
#include <string>
#include <fstream>
#include <sstream>

namespace elSpectro{

  class HepMC3Writer : public Writer {

     
     
     HepMC3Writer()=default;
     //don't want default contructor accessible
     //only declaring default constructor
     //so other 5 constructors also defaulted(rule of 5)

   public:
     HepMC3Writer(const std::string& filename);
     ~HepMC3Writer() final;
     HepMC3Writer(const HepMC3Writer& other); //need the virtual destructor...so rule of 5
     HepMC3Writer(HepMC3Writer&&)=default;
     HepMC3Writer& operator=(const HepMC3Writer& other);
     HepMC3Writer& operator=(HepMC3Writer&& other) = default;
     
     void WriteHeader() final{};
     void FillAnEvent() final;
     void Write() final;
     void End() final;
     
     void Init() final;
     
   private:
     decaying_constptrs _vertexParticles;
 
     //streaming functions
     
     void StreamEventInfo(){
       //E = event number, # vertices, # particles = initial+final
       _stream<< "E "<<" "<<_nEvent<<" "<<_vertices->size()<<
	 " "<<_initialParticles.size()+_finalParticles.size()+_vertexParticles.size()<<"\n";
     }
     /////////////////////////////////////////////////////////
     void StreamEventPosition();
     /////////////////////////////////////////////////////////
     void StreamUnits(){
       _stream<< "U GEV MM"<<"\n";
     }
     void StreamWeights(){}
     void StreamAttributes(){
       /*_stream<< "A 0 Q2 "<<"\n";
      _stream<< "A 0 W "<<"\n";
       */
     }
     void StreamParticle(const Particle* p,int vertex_id,int status){
       auto p4=p->P4();
       _stream<<"P "<<_id++<<" "<<-(vertex_id+1)<<" "<<p->Pdg()<<" "
	      <<p4.X()<<" "<<p4.Y()<<" "<<p4.Z()<<" "<<p4.T()
	      <<" "<<p->Mass()<<" "<<status<<"\n";
     }
     void StreamVertex(int vertex_id,int status,std::vector<int> in_pids){
       _stream<<"V "<<-(vertex_id+1)<<" "<<status<<" [";

       uint ip=0;
       
       if(!in_pids.empty()){
	 for(ip=0;ip<in_pids.size()-1;ip++){
	   _stream<<in_pids[ip]<<",";
	 }
	 _stream<<in_pids[ip]<<"]";
       }
 
       if(in_pids.empty())
	 _stream<<"]";
       
       const auto pos=_vertices->at(vertex_id);
       //if(pos->M()!=0){
	 _stream<<" @ "<<pos->X()<<" "<<pos->Y()<<" "
		<<pos->Z()<<" "<<pos->T();
	 // }
       _stream<<"\n";
     }
     
     //data members
     std::ofstream _file; //! output file
     //std::ostream _stream; //! output stream
     std::stringstream _stream;

     long _nEvent={0};
   
     int _id=1;
     
     ClassDef(elSpectro::HepMC3Writer,1); //class Writer
   };


}

