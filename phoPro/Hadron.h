#pragma once
//#include "Interface.h"
#include "DecayChannel.h"
#include <TNamed.h>


namespace phoPro{

  using rawParticle =  elSpectro::Particle*;

  class Hadron : public TNamed {
    
  
  public:
    Hadron() = default;
    
    Hadron(const TString& name, const TString& type=""):TNamed(name,name){
      //Hadron(const TString& name){
      //Hadron(const char* name,const TString& type=""){
      // Hadron(_pdg,type);

      
      if(TDatabasePDG::Instance()->GetParticle(name)){
      
	_pdg = TDatabasePDG::Instance()->GetParticle(name)->PdgCode();
	Hadron(_pdg,type);
	
	//Hadron(_pdg,"");
     
      }
      else{
	std::cerr<<"Hadron::Hadron Warning no pdgcode for "<<name<<std::endl;
	//	exit(0);
      }
      
    }
    
    Hadron(int pdg,const TString& type=""):_pdg(pdg){
      //check if registered in PArticleManager
      // _registered = elSpectro::particles().CheckUsedPdg(pdg);
      //cout<<"Hadron "<<pdg<<endl;
      //Check if exits in PDG table
      // if(elSpectro::particles().CheckInPdgTable(pdg)){
      // 	auto rootPdg = TDatabasePDG::Instance()->GetParticle(pdg);
      // 	cout<<"Hadron "<<rootPdg<<endl;
      // 	_mass = rootPdg->Mass();
      // 	_width = rootPdg->Width();
      // 	SetNameTitle(rootPdg->GetName(),rootPdg->GetName());
      // 	auto decays = rootPdg->DecayList();
      // 	for(auto* odecay:*decays){
      // 	  auto decay  = static_cast<TDecayChannel*>(odecay);
      // 	  std::vector<Hadron> products;
      // 	  for(auto ip = 0; ip<decay->NDaughters();++ip){
      // 	    products.push_back( Hadron(decay->DaughterPdgCode(ip)) );
      // 	  }
      // 	  AddDecay(decay->BranchingRatio(),products);
      // 	}
      // }
      // else{ //need to set mass later
      // 	elSpectro::particles().RegisterNewPdgParticle(Mass(),nullptr,type,PDGCode());
      // }
 
      
    }

    // virtual const TString& Type() = 0;
    
    /* void AddDecay(double branch_ratio,elSpectro::particle_ptrs parts,std::vector<int> pdgs) { */
    /*   _decayChannels.push_back(elSpectro::particle(PDGCode(), */
    /* 		     model(new elSpectro::PhaseSpaceDecay(parts , pdgs) ) )); */
    /*   _decayChannels.back()->SetPdgMass(Mass()); */
    /*   _branchRatio.push_back(branch_ratio); */
    /* } */

    void AddDecay(double branch_ratio,std::vector<Hadron> hadrons){
      // std::cout<<"Hadron::AddDecay "<<branch_ratio<<" "<<hadrons.size()<<std::endl;
      for(auto& h:hadrons)
	if(h.PDGCode()==0){
	  std::cout<<"Hadron::AddDecay Unknown particle, will not add this decay "<<std::endl;
	  return;
	}
	  
      _decayChannels.push_back(hadrons);
      _branchRatio.push_back(branch_ratio);
      _decayPtrs.push_back(nullptr);
    }
    
    rawParticle  GetDecay(std::vector<uint> vchn={}) const{
      std::cout<<"GetDecay adding # "<<_decayChannels.size()<<std::endl;

      //Note if no list of channels is given will default
      //to {0,0,0,...} ie. the first specified channel in each decay
      if(_decayChannels.empty()) return nullptr;
      if(_decayPtrs[0]!=nullptr) return _decayPtrs[0];
      
      //_decayChannels vector of decays
      //decay is a vector of particles
      elSpectro::particle_ptrs dec_parts;
      std::vector<int> stable_parts; //if particles don't decay just use PDG
      
      uint ichn = 0; //first channel by default
      /* 
      //Below will not work. Need to be able to request
      //decay channels for each product in decay....
      //stick with only 1 decay channel for now.

      if(vchn.empty()==false){
      //take the first index for this channel
      auto ichn = vchn.front();
      if(ichn>=_decayChannels.size()){
      std::cerr<<"GetDecay invalid channel "<<ichn<<std::endl;
      exit(0);
      }
      //remove to pass on to next
      vchn.erase(vchn.begin());
      }
      */
      
      //get products for requested channel
      auto myproducts = _decayChannels[ichn];
      std::cout<<"GetDecay adding products "<<myproducts.size()<<std::endl;
      //loop over products, applying their decays if they have any
      for(const auto& product : myproducts ){
	std::cout<<"GetDecay product "<<product.PDGCode()<<std::endl; 
	auto decaying_particle  = product.GetDecay(vchn);
	if(decaying_particle!=nullptr){
	  dec_parts.push_back(decaying_particle);
	  std::cout<<"GetDecay adding decay particle "<<decaying_particle->Pdg()<<std::endl;
	}
	else{
	  stable_parts.push_back(product.PDGCode());
	}
      }
    
      
      rawParticle  p =
	elSpectro::particle(PDGCode(),
			    model(new elSpectro::PhaseSpaceDecay(dec_parts,stable_parts) ) );
      
      //set elSpectro particle mass to PDG
      p->SetPdgMass(Mass());
      _decayPtrs[0]=p;
      
      return p;
    }
    
    // rawParticle DecayChannel(uint i=0)const {return _decayChannels[i];}

    
    int PDGCode()const {return _pdg;}
    bool IsRegistered() const{return _registered;}

    bool ShouldItDecay() const {
      //particle may decay if short lived ( <1ns)
      auto rootPdg = TDatabasePDG::Instance()->GetParticle(GetName());
      return  rootPdg->Lifetime()<1E-9&&rootPdg->Lifetime()>0;
    }
    bool IsWide() const{
      return Width()>0.001; //minimum width 1MeV
    }
    void RegisterParticle(elSpectro::Distribution* dist=nullptr){
      if(IsRegistered()==false){
	std::cout<<" RegisterParticle "<< PDGCode()<<" "<<Mass()<<dist<<std::endl;
	if(Mass()==0){
	  std::cout<<" RegisterParticle can't have zero mass hadron"<<std::endl;
	  exit(0);
	}
	
	if(elSpectro::particles().CheckInPdgTable(PDGCode())){
	  if(dist!=nullptr)
	    elSpectro::particles().RegisterMassDistribution(PDGCode(),dist);
	  // else in pdg table and no mass distribution 
	}
	else{
	  elSpectro::particles().RegisterNewPdgParticle(Mass(),dist,PDGCode());
	}
	_registered = true;
      }
    }
    void SetMass(double m){_mass = m;}
    double Mass()const {return _mass;}

    void SetWidth(double w){_width=w;};
    double Width()const {return _width;}

    bool IsDecaying(){
      return _decayChannels.empty()!=true ;
    }
    uint NDecays() const{_decayChannels.size();}
    
    elSpectro::decaymodel_ptr ChannelModel(uint ichn) const{
      //determine what model decays into
      elSpectro::decaying_objs dec_parts;
      elSpectro::particle_objs stable_parts;
      
      //get products for requested channel
      auto chan_products = _decayChannels[ichn];

      for(uint i=0;i<chan_products.size();++i){
	auto h = chan_products[i];
	if (h.ShouldItDecay() && h.IsDecaying() ){
	  auto dproduct = CreateDecayingParticle(h.PDGCode());
	  dec_parts.push_back(dproduct);
	}
	else{
	  stable_parts.push_back(h.PDGCode());
	}
      }
      elSpectro::particle_ptrs ptrs;
      for(auto& p:dec_parts){
	ptrs.push_back(&p);
      }
      //CloneModel will copy all particles locally
      //dec_parts and ptrs will go out of scope here
      return elSpectro::CloneModel(elSpectro::PhaseSpaceDecay(ptrs,stable_parts) );
  
    }
 
    double BranchRatio(uint ichn) const{
      return _branchRatio[ichn];
    }
    
    std::unique_ptr<DecayVectors> ChannelDecayer(uint ichn) const{

    }
    
    //combined branching ratio of all decays
    double BranchRatio() const{
      if(_decayChannels.empty()) return 1.;
      /*
	uint ichn = 0; //first channel by default
	if(vchn.empty()==false){
	//take the first index for this channel
	auto ichn = vchn.front();
	if(ichn>=_decayChannels.size()){
	std::cerr<<"GetDecay invalid channel "<<ichn<<std::endl;
	exit(0);
	}
	//remove to pass on to next
	vchn.erase(vchn.begin());
	}
      */
      uint ichn = 0;
      
      double br = _branchRatio[ichn];
      
      cout<<"BranchRatio  "<<_decayChannels.size()<<endl;
      for(const auto& product :_decayChannels[ichn]){
	cout<<product.PDGCode()<<endl;
	br *= product.BranchRatio();
      }
      cout<<"BranchRatio br "<<br<<endl;
    
      return br;
    }
    
    elSpectro::DistTF1* BreitWignerDistribution(){
      if( IsWide()==false ) return nullptr;
      
      return new elSpectro::DistTF1{TF1("Mass",Form("TMath::BreitWigner(x,%lf,%lf)",Mass(),Width()),Mass()-5*Width(),Mass()+5*Width())};
    }

    void Print(){
      std::cout<<"Hadron "<< GetName()<<" "<<PDGCode() <<" "<<Mass()<<" "<<Width()<<std::endl;
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
      std::cout<<"\t total branching ration : "<<sumBranches<<std::endl;
    }
      
  private:
    
    std::vector< std::vector<Hadron> > _decayChannels;//
    mutable elSpectro::particle_ptrs _decayPtrs;//
    
    std::vector<double> _branchRatio;//
    double _mass=0;//
    double _width=0;//
    int _pdg=0;//
    bool _registered = false;//


    ClassDef(phoPro::Hadron,1);
  };
  
  /* class Meson : public Hadron { */

  /* public: */
    
  /* Meson(int pdg):Hadron(pdg,Type()){}; */

  /*   const TString& Type() override {return _type;} */
  
  /* private: */
  /*   TString _type="Meson"; */
  /* }; */

  /* class Baryon : public Hadron { */

  /* public: */
  /* Baryon(int pdg):Hadron(pdg,Type()){}; */

  /*   const TString& Type() override {return _type;} */

  /* private: */

  /*   TString _type="Baryon"; */
  /* }; */

 

}
#pragma link C++ class phoPro::Hadron+;
