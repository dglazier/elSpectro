//Need to automate this from downloadable PDG data
//https://pdg.lbl.gov/2022/html/computer_read.html
//which gives properties and decays and
//total photoproduction cross section
//**Note this turned out to not be viable at the moment
//** branching decays are not well defined

//Will not need specific classes for each particle
//Just Meson and Baryon with everything configured
//from the PDG.

//May need to add some exotic particles via similar
//machanism.

#include "Hadron.h"
#include "VectorMesons.h"
#include "XMesons.h"
#include "ZMesons.h"
#include "YMesons.h"
#include <TFile.h>


void ParticleFactory(){
}


namespace particles{

 
  class Factory {

  public:
    
    static Factory& Instance() { static Factory instance; return instance; }
 
    
    phoPro::Hadron  Get(int pdg){
      return Get(TDatabasePDG::Instance()->GetParticle(pdg)->GetName());
    }

    phoPro::Hadron  Get(const std::string& name){
      std::string type = "Meson";
      if(TDatabasePDG::Instance()->GetParticle(name.data()) == nullptr){
	std::cerr<< "phoPro::Hadron  Get " <<name<<" does not exist in database"<<std::endl;
	exit(0);
      }

      //if a new particle type register its properties
      if(std::find(_names.begin(),_names.end(),name) == _names.end()){
	
	auto p = phoPro::Hadron(name);
	auto rootPdg = TDatabasePDG::Instance()->GetParticle(name.data());
 	p.SetMass(rootPdg->Mass());
	p.SetWidth(rootPdg->Width());
	p.SetNameTitle(rootPdg->GetName(),rootPdg->GetName());
	if(p.ShouldItDecay()){ //check lifetime
	  
	  auto decays = rootPdg->DecayList();
	  for(auto* odecay:*decays){
	    auto decay  = static_cast<TDecayChannel*>(odecay);
	    std::vector<phoPro::Hadron> products;
	    for(auto ip = 0; ip<decay->NDaughters();++ip){
	      if(decay->DaughterPdgCode(ip)==1)exit(0);
	      products.push_back( Get(decay->DaughterPdgCode(ip)) );
	    }
	    p.AddDecay(decay->BranchingRatio(),products);
	  }
	  
	}
	p.Print();
	_existing.push_back(p);
	_names.push_back(p.GetName());
	//	elSpectro::particles().RegisterNewPdgParticle(p.Mass(),p.BreitWignerDistribution(),rootPdg->ParticleClass(),p.PDGCode());
	return p;
      }
      //particle already registered, return a copy
      auto it = find(_names.begin(), _names.end(), name); 
      auto entry = std::distance(_names.begin(),it);
      return _existing[entry];
    }
    
    DecayingParticle CreateDecayingParticle(int pdg){//,decaymodel)
      auto p = Get(pdg);
      //check width/lifetime and for decay channels
      //same condition should be in ChannelModel
      if( p.ShouldItDecay() && p.IsDecaying() ){

	DecayingParticle dp{pdg};
	for(uint i=0;i<p.NDecays();++i){
	  dp.AddDecay(BranchRatio(i),ChannelModel(i),ChannelDecayer(i));
	}
	return dp;
      }
      else{
	cerr<<"ParticleFactory::CreateDecayingParticle "<<pdg <<" does not decay!"<<std::endl;
	exit(0);
	return DecayingParticle(0);
      }
    }
    
  private:
    std::vector<phoPro::Hadron> _existing;
    std::vector<std::string> _names;
  };
  
}
