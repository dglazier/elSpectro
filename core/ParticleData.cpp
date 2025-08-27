#include "ParticleData.h"
#include "ParticleFactory.h"
#include "PhaseSpaceDecay.h"

namespace elSpectro{
  namespace particles{
    
    ParticleData::ParticleData(const TString& name, const TString& type):TNamed(name,name){
      
	if(TDatabasePDG::Instance()->GetParticle(name)){
	
	  _pdg = TDatabasePDG::Instance()->GetParticle(name)->PdgCode();
	  ParticleData(_pdg,type);
	
	}
	else{
	  std::cerr<<"ParticleData::ParticleData Warning no pdgcode for "<<name<<std::endl;
	}    
      }
 
    ////////////////////////////////////////////////////////////
    void ParticleData::AddDecay(double branch_ratio,std::vector<ParticleData> hadrons){
      for(auto& h:hadrons)
	if(h.PDGCode()==0){
	  std::cout<<"ParticleData::AddDecay Unknown particle, will not add this decay "<<std::endl;
	  return;
	}
	  
      _decayChannels.push_back(hadrons);
      _branchRatio.push_back(branch_ratio);
      _decayPtrs.push_back(nullptr);
    }
    ///////////////////////////////////////////////////////////////
    decaymodel_ptr ParticleData::ChannelModel(uint ichn) const{
  	//determine what model decays into
	elSpectro::decaying_objs dec_parts;
	elSpectro::particle_objs stable_parts;
      
	//get products for requested channel
	auto chan_products = _decayChannels[ichn];

 	for(uint i=0;i<chan_products.size();++i){
	  auto h = chan_products[i];
	  if (h.ShouldItDecay() && h.IsDecaying() ){
	    auto dproduct = ParticleFactory::Instance().CreateDecayingParticle(h.PDGCode());
	    dec_parts.push_back(std::move(dproduct));
	  }
	  else{

	    stable_parts.push_back(h.PDGCode());
	  }
	}
	//CloneModel will copy all particles locally
	//dec_parts and ptrs will go out of scope here
	auto res = elSpectro::CloneModel(elSpectro::PhaseSpaceDecay(dec_parts,stable_parts) );
	
	return res;
  
      }
    ///////////////////////////////////////////////////////////////////
    //combined branching ratio of all decays
    double ParticleData::BranchRatio() const{
      
      if(_decayChannels.empty()) return 1.;
      
      uint ichn = 0;
      
      double br = _branchRatio[ichn];
      
      for(const auto& product :_decayChannels[ichn]){
	br *= product.BranchRatio();
      }
      
      return br;
    }
    //////////////////////////////////////////////////////////////////
     void ParticleData::Print(Option_t *option) const{
       std::cout<<"ParticleData "<< GetName()<<" "<<PDGCode() <<" mass = "<<Mass()<<" width = "<<Width()<<" lifetime = "<<Lifetime()<<"("<< MeanFreePath()<<" mm)"<<std::endl;
	std::cout<<"\t Decays :"<<std::endl;
	uint idec=0;
	double sumBranches = 0;
	for(auto& dec:_decayChannels){
	  std::cout<<"\t\t"<<idec<<" ";
	  for(auto& p:dec){
	    std::cout<<p.PDGCode()<<" ";
	  }
	  std::cout<<_branchRatio[idec]<<std::endl;
	  sumBranches += _branchRatio[idec];
	  idec++;
	}
	std::cout<<"\t total branching ratio : "<<sumBranches<<std::endl;
      }
  }
}
