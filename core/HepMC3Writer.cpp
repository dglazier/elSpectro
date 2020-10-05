#include "HepMC3Writer.h"
#include "Manager.h"
#include <iostream>

namespace elSpectro{

  ///Constructor to create ouput file and intialise data structures
  HepMC3Writer::HepMC3Writer(const std::string &filename):
    _file(filename)
  {
    if(!_file.is_open()){
      std::cerr<<"HepMC3Writer::HepMC3Writer file "<<filename<<" cannot be opened, exiting..."<<std::endl;
      exit(0);
    }
    //Give version used when this code was written
    _stream << "HepMC::Version 3.02.02"  << std::endl;
    _stream << "HepMC::Asciiv3-START_EVENT_LISTING" << std::endl;
    Write();

  }
  HepMC3Writer::~HepMC3Writer(){
    End();
  }
  ///////////////////////////////////////////////////////////////
  ///Close the file stream
  void HepMC3Writer::End(){
    if(!_file.is_open()) return;
    _stream << "HepMC::Asciiv3-END_EVENT_LISTING" << std::endl << std::endl;
    Write();
    _file.close();

  }
  void HepMC3Writer::Init(){
    Writer::Init();
    //also
    //need to find detached vertices
    auto& unptrs=Manager::Instance().Particles().UnstableParticles();
    for(const auto* p:unptrs){
      //in case we have a detached vertex we need to write this particle
      if(p->IsDecay()==DecayType::Detached)_vertexParticles.push_back(p);
    }

  }
  /////////////////////////////////////////////////////////////
  //write all the info required for this event
  void HepMC3Writer::FillAnEvent(){

    
    
    ////fill _stream
    StreamEventInfo();
    StreamEventPosition();
    StreamUnits();
    StreamWeights();
  
    _id=1;//reset particle ID counter
    //initial particles
    auto nVer = _vertices->size();//last stored vertice = primary
    int initial_status=3;
    int initial_vertex_id=-1;

     for(const auto* p:_initialParticles)
      StreamParticle(p,initial_vertex_id,initial_status);
  
    //primary reaction vertex
    int primary_vertex_id=0;
    int primary_vertex_status=0;
    std::vector<int> primary_parent_ids={1,2};
    StreamVertex(primary_vertex_id,primary_vertex_status,primary_parent_ids);
  
    //final particles
    int final_status=1;
    std::vector<std::pair<int,int>> writtenVertexParticles;
    for(auto iver=0;iver<nVer;++iver){
  
      int final_vertex_id=iver;//numbers from -1,-2,...
      int final_vertex_status=0;

      //first stream vertex decaying particles
      for(const auto* p:_vertexParticles){
	if(p->VertexID()==iver){
	  int decay_particle_status=3;
	  //save decay vertex ID and this particle id
	  writtenVertexParticles.push_back(std::pair<int,int>(p->DecayVertexID(),_id));//before _id incremened
	  StreamParticle(p,final_vertex_id, decay_particle_status);
	}
      }
      bool first = (final_vertex_id != 0) ;//true if not production vertex
      for(const auto* p:_finalParticles){
	
	if(p->VertexID()==iver){ //is this particle from this vertex?
	  
	  //write vertex info first time there is a particle
	  if(first==true){//see if need to write other vertex
	   
	    auto findVertex=[&final_vertex_id, &writtenVertexParticles](){
	      for(auto const& findVertex: writtenVertexParticles){
		if(findVertex.first==final_vertex_id)
		  return findVertex.second;

	      }
	      return -1;//default, shouldn't ever happen!
	    };
	    auto decayParticleID=findVertex();
	    
	    StreamVertex(final_vertex_id,final_vertex_status,{decayParticleID});
	    first=false;
	  }

	  StreamParticle(p,final_vertex_id,final_status);
	}
      }
    }
    
    _nEvent++;
  }

  /////////////////////////////////////////////////////////
 
  /////////////////////////////////////////////////////////
  void HepMC3Writer::Write(){
    //write stream
    _file<<_stream.rdbuf();
    //reset stream
    _stream.clear();
  }

  /////////////////////////////////////////////////////////
  void HepMC3Writer::StreamEventPosition(){
    //from HepMC3::WriterAscii
    /* const FourVector &pos = evt.event_pos();
    if ( !pos.is_zero() ) {
        m_cursor += sprintf(m_cursor," @ %.*e",m_precision,pos.x());
        flush();
        m_cursor += sprintf(m_cursor," %.*e",  m_precision,pos.y());
        flush();
        m_cursor += sprintf(m_cursor," %.*e",  m_precision,pos.z());
        flush();
        m_cursor += sprintf(m_cursor," %.*e",  m_precision,pos.t());
        flush();
    }

    m_cursor += sprintf(m_cursor,"\n");
    */
  }
  
}
