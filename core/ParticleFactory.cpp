#include "ParticleFactory.h"
#include <TDecayChannel.h>

namespace elSpectro{

  namespace particles{

    ParticleData  ParticleFactory::GetData(const std::string& name){
      std::string type = "Meson";
      if(TDatabasePDG::Instance()->GetParticle(name.data()) == nullptr){
	std::cerr<< "ParticleData  Get " <<name<<" does not exist in database"<<std::endl;
	exit(0);
      }

      //if a new particle type register its properties
      if(std::find(_names.begin(),_names.end(),name) == _names.end()){
	
	auto p = ParticleData(name);
	auto rootPdg = TDatabasePDG::Instance()->GetParticle(name.data());
	p.SetMass(rootPdg->Mass());
	p.SetWidth(rootPdg->Width());
	p.SetNameTitle(rootPdg->GetName(),rootPdg->GetName());
	if(p.ShouldItDecay()){ //check lifetime
	  
	  auto decays = rootPdg->DecayList();
	  for(auto* odecay:*decays){
	    auto decay  = static_cast<TDecayChannel*>(odecay);
	    std::vector<ParticleData> products;
	    for(auto ip = 0; ip<decay->NDaughters();++ip){
	      if(decay->DaughterPdgCode(ip)==1)exit(0);
	      products.push_back( GetData(decay->DaughterPdgCode(ip)) );
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
    //////////////////////////////////////////////////////////////////
    DecayingParticle ParticleFactory::CreateDecayingParticle(const std::string& name){
    if(TDatabasePDG::Instance()->GetParticle(name.data()) == nullptr){
	std::cerr<< "ParticleData  CreateDecayingParticle " <<name<<" does not exist in database"<<std::endl;
	exit(0);
      }
      return ParticleFactory::CreateDecayingParticle(TDatabasePDG::Instance()->GetParticle(name.data())->PdgCode());
    }
    //////////////////////////////////////////////////////////////////
    DecayingParticle ParticleFactory::CreateDecayingParticle(int pdg){//,decaymodel)
      auto pdata = GetData(pdg);
      std::cout<<"DecayingParticle ParticleFactory got pdata "<<pdg<<std::endl;
      //check width/lifetime and for decay channels
      //same condition should be in ChannelModel
      if( pdata.ShouldItDecay() && pdata.IsDecaying() ){
	
	DecayingParticle dp{pdg};
	dp.SetMassDist( cpp::MakeBaseShared<Distribution>(pdata.BreitWignerDistribution()) );
	std::cout<<"DecayingParticle ParticleFactory got pdata Add decays "<<std::endl;
 	
	//	for(uint i=0;i<pdata.NDecays();++i){
	for(uint i=0;i<1;++i){
	  std::cout<<"DecayingParticle ParticleFactory  "<<pdg<<" "<<i<<" "<<pdata.NDecays()<<" "<<pdata.BranchRatio(i)<<std::endl;
 
	  dp.AddDecay(pdata.BranchRatio(i),pdata.ChannelModel(i),pdata.ChannelDecayer(i));
	  std::cout<<"DecayingParticle ParticleFactory  decay added "<<std::endl;
	}
	std::cout<<"DecayingParticle ParticleFactory got pdata can return "<<std::endl;
	return dp;
      }
      else{
	std::cerr<<"ParticleFactory::CreateDecayingParticle "<<pdg <<" does not decay!"<<std::endl;
	exit(0);
	return DecayingParticle(0);
      }
    }
  }

}
