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
	p.SetLifetime(rootPdg->Lifetime());
	p.SetNameTitle(rootPdg->GetName(),rootPdg->GetName());

	if(p.ShouldItDecay()){ //check lifetime

	  //register all its decays
	  auto decays = rootPdg->DecayList();
	  for(auto* odecay:*decays){
	    auto decay  = static_cast<TDecayChannel*>(odecay);
	    std::vector<ParticleData> products;
	    //register all its children
	    for(auto ip = 0; ip<decay->NDaughters();++ip){
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
      // std::cout<<"DecayingParticle ParticleFactory  "<<pdata.ShouldItDecay()<<" "<<pdata.IsDecaying()<<" w "<<pdata.Width()<<" mfp "<<pdata.MeanFreePath()<<" "<<pdata.PDGCode()<<std::endl;
     //std::cout<<"DecayingParticle ParticleFactory got pdata "<<pdg<<std::endl;
      //check width/lifetime and for decay channels
      //same condition should be in ChannelModel
      if( pdata.ShouldItDecay() && pdata.IsDecaying() ){
	DecayingParticle dp{pdg};
	if(pdata.IsWide()){
	  //delta mass distribution if width < 0.1MeV

	  if(TMath::Abs(pdata.PDGCode())==82){
	    //special case of rndmflav particle
	    //this decays randomly to multi-mesons
	    //should have flat mass distribution
	    dp.SetMassDist( cpp::MakeBaseShared<Distribution>( elSpectro::DistTF1{TF1("Mass",Form("1"),pdata.Mass(),pdata.Mass()+pdata.Width())} ));
	  }
	  
	  else{
	    dp.SetMassDist( cpp::MakeBaseShared<Distribution>(pdata.BreitWignerDistribution()) );
	  }
	  
	}
	
	//	std::cout<<"DecayingParticle ParticleFactory got pdata Add decays "<<std::endl;
	if( pdata.MeanFreePath()>0.01 ){ //0.1mm
	  
	  auto decDist=DistTF1(TF1(Form("VertexFor%s",pdata.GetName()),"TMath::Exp(-x/[0])",0,25*pdata.Lifetime()));//in s
	  
	  decDist.GetTF1().SetParameter(0,pdata.Lifetime());
	  decDist.GetTF1().SetNpx(200);

	  dp.SetDecayVertexDist(std::move(decDist));
	}
	
	for(uint i=0;i<pdata.NDecays();++i){
	  //	for(uint i=0;i<1;++i){
	  //std::cout<<"DecayingParticle ParticleFactory  "<<pdg<<" "<<i<<" "<<pdata.NDecays()<<" "<<pdata.BranchRatio(i)<<std::endl;
 
	  dp.AddDecay(pdata.BranchRatio(i),pdata.ChannelModel(i),pdata.ChannelDecayer(i));
	  // std::cout<<"DecayingParticle ParticleFactory  decay added "<<std::endl;
	}
      	//std::cout<<"DecayingParticle ParticleFactory got pdata can return "<<pdata.GetName()<<std::endl;
	return (dp);
      }
      else{
	std::cerr<<"ParticleFactory::CreateDecayingParticle "<<pdg <<" does not decay!"<<std::endl;
	exit(0);
	return DecayingParticle(0);
      }
    }

  }
}
