#include "PhaseSpaceDecay.h"
#include "Interface.h"
#include "ParticleFactory.h"
#include "DistFlatMass.h"


namespace elSpectro{

  //  PhaseSpaceDecay::PhaseSpaceDecay( particle_ptrs ps, const std::vector<int> pdgs):
  PhaseSpaceDecay::PhaseSpaceDecay(const decaying_objs& decs, const particle_objs& stables):
    DecayModel{decs,stables}
  {
    //    std::cout<<"PhaseSpaceDecay::PhaseSpaceDecay, stable =  "<<StableProducts().size()<<" unstable = "<<UnstableProducts().size()<<" all = "<<Products().size()<<std::endl;
    _name={"PhaseSpaceDecay"};

  }

  void PhaseSpaceDecay::SetParentAndProducts(DecayingParticle* pa,  particle_objs& stable, decaying_objs& unstable){
    DecayModel::SetParent(pa);
    nBodyDecayer(pa,stable,unstable);

    //std::vector<double> masses;
    //GetStableMasses(masses);
    // for(auto& mass:masses){
    //   std::cout<<" stable mass "<<mass;
    // }
    // std::cout<<"\n\n\n*************************************"<<std::endl;
  }
  
  void PhaseSpaceDecay::SetParent(DecayingParticle* pa){
    DecayModel::SetParent(pa);
    //now can make cascade of decays
    bool size2=false;
    if( Products().size()>2 ){
      nBodyDecayer(pa,MutableStableProducts(),MutableUnstableProducts());
      size2=true;
    }
    // for(auto& p:Products()){
    // 	std::cout<<" PhaseSpaceDecay::SetParent "<<pa->Pdg()<<std::endl;
    // 	std::vector<double> masses;
    // 	GetStableMasses(masses);
    // 	for(auto& mass:masses){
    // 	  std::cout<<" stable mass "<<mass;
	
    // 	std::cout<<"\n\n\n*************************************"<<std::endl;
    //   }
    // }
  }
  void PhaseSpaceDecay::PostInit(ReactionInfo* info){
    auto maxmass=Parent()->MaximumMassPossible();
    if(maxmass==0) maxmass = info->Wmax(); //if parent has no max just use reaction W
    auto flatmass=dynamic_cast<DistFlatMass*>(Parent()->MassDistribution());
    if(flatmass) flatmass->SetMaxX(maxmass);
  
    for(auto& product:Products()){
      auto flatmass=dynamic_cast<DistFlatMass*>(product->MassDistribution());
      SumAllProducts();
      if(flatmass) flatmass->SetMaxX(maxmass);
    }
    
    if(_massMaster!=nullptr){
      _massMaster->SetMaxX(maxmass);
      _massMaster->SetParentPtr(Parent());
      SetMassMaster(_massMaster);
      // Print();
    }
 
    DecayModel::PostInit(info);

  }
  
  void PhaseSpaceDecay::nBodyDecayer(DecayingParticle* parent,particle_objs& stable,decaying_objs& unstable ) //take copies of particle vectors
  {
    // std::cout<<"PhaseSpaceDecay::nBodyDecayer  for "<<parent->Pdg()<<std::endl;
   //  for(auto& p:stable){
   //    std::cout<<p.Pdg()<<" ";
   //  }
   // for(auto& p:unstable){
   //    std::cout<<p.Pdg()<<" ";
   //  }
   // std::cout<<std::endl;
    //Make N-1 particles decaying into 2 bodies
    //Provide mass distribution

    //create and register all final particles
    //add unstable at start so more room for selecting their mass
    particle_ptrs ps;
    for(auto& unp: unstable)
      ps.push_back( dynamic_cast<Particle*>(&unp) );
    for(auto& sp: stable)
      ps.push_back( &sp );
    
    particle_ptrs ptrs; //new pointer list, when they are copied in sub-models
    decaying_ptrs Xm1s; //for sub X-n systems
 
    //create our flat phase space mass distribution master. 
    auto massMaster = DistFlatMassMaster{parent,ps};
  
    //    auto massMaster  = mass_distribution(_nextPdg,new DistFlatMassMaster(parent,ps));
    // _nextPDG++;
    //Start with X(2)
    //auto XN = particleMan.Take(new DecayingParticle{_nextPdg,new PhaseSpaceDecay(ps[0],ps[1]),new TwoBodyFlat()});
   
    //we are given X(N)
    // X(N)-> X(N-1) + p1
    //            ->X(N-2) + p2
    //                 -> X(N-3) + p3
    //                   ...
    //                         -> X(2) + pN-2
    //                                 ->pN-1  + pN
    
    //std::unique_ptr<Particle> XNminusM(ps[0]);

    //start with first product as X(N-1)
    std::unique_ptr<Particle> XNminusM;
    if(dynamic_cast<DecayingParticle*>(ps[0])) XNminusM.reset(new DecayingParticle(*dynamic_cast<DecayingParticle*>(ps[0])));
    else{
      XNminusM.reset(new Particle(*ps[0]));
    }
   
    auto N=ps.size();
    //if we have >3 products iterate through X(N-M)
    //The final one will be X(2).
    //If we only have 3 products then XNminusM is already X(2) and we do not need to iterate
    //there will be no intermediate X states
    for(uint m=1;m<=N-3;++m){//leave 2 particles
      auto p = ps[m];

      //PhaseSpaceDecay will copy the particles we give
      //so need to get ptrs to them for mass distribution after
      particle_objs pproducts;
      decaying_objs dproducts;
      //put copies into products vectors
      if(dynamic_cast<DecayingParticle*>(XNminusM.get())){
	dproducts.push_back(*dynamic_cast<DecayingParticle*>(XNminusM.get()));
	//	std::cout<< "***************************PhaseSpaceDecay:: add decaying product XNminusM "<<XNminusM->Pdg()<<std::endl;
      }
      else {
	pproducts.push_back(*XNminusM.get());
	//std::cout<< "***************************PhaseSpaceDecay:: add  product XNminusM "<<XNminusM->Pdg()<<std::endl;
       
      }
      if(dynamic_cast<DecayingParticle*>(p)){
	dproducts.push_back(*dynamic_cast<DecayingParticle*>(p));
	//std::cout<< "***************************PhaseSpaceDecay:: add decaying product p "<<p->Pdg()<<std::endl;
      }
      else{
	pproducts.push_back(*p);
	//std::cout<< "***************************PhaseSpaceDecay:: add  product p "<<p->Pdg()<<std::endl;
      }
      //create new iteration of X-m
      XNminusM.reset(new DecayingParticle{elSpectro::particles::inter_phasespace,
					  CloneModel(PhaseSpaceDecay{dproducts,pproducts} ),
					  CloneDecayer(TwoBodyFlat())});
      XNminusM->SetMassDist(cpp::MakeBaseShared<Distribution>(DistFlatMass(&massMaster)));
      
 
      //addd the standard particles to ptrs, that do not get decayed by this algorithm
      if(dynamic_cast<DecayingParticle*>(XNminusM.get())->Model()->MutableProduct(0)->Pdg()!=elSpectro::particles::inter_phasespace) {
	ptrs.push_back(dynamic_cast<DecayingParticle*>(XNminusM.get())->Model()->MutableProduct(0));
      }
      if(dynamic_cast<DecayingParticle*>(XNminusM.get())->Model()->MutableProduct(1)->Pdg()!=elSpectro::particles::inter_phasespace) {
	ptrs.push_back(dynamic_cast<DecayingParticle*>(XNminusM.get())->Model()->MutableProduct(1));
      }
      
      // std::cout<<"PhaseSpaceDecay Add ptr  "<<dynamic_cast<DecayingParticle*>(XNminusM.get())->Model()->MutableProduct(0)->Pdg()<<" "<<dynamic_cast<DecayingParticle*>(XNminusM.get())->Model()->MutableProduct(1)->Pdg()<<std::endl;
      // std::cout<<"PhaseSpaceDecay:: XNminusM products ";
      // for(auto& pp:dproducts) std::cout<<" d "<<pp.Pdg();
      // for(auto& pp:pproducts) std::cout<<" p "<<pp.Pdg();
      // std::cout<<std::endl;
    }

    
    //give master to the X(N-1) state
    auto p = ps[N-2]; //m=1
    //combine with the X(N-2) particle and give pdgXNm1 mass distribution
    //std::cout<<"PhaseSpaceDecay:: final  create new two-body "<<" from "<<p->Pdg()<<" and "<<XNminusM->Pdg()<<std::endl;
    particle_objs pproducts;
    decaying_objs dproducts;
    if(dynamic_cast<DecayingParticle*>(XNminusM.get())){
      dproducts.push_back(*dynamic_cast<DecayingParticle*>(XNminusM.get()));
      //std::cout<< "***************************PhaseSpaceDecay:: add decaying product XNminusM "<<XNminusM->Pdg()<<std::endl;
    }
    else{
      pproducts.push_back(*XNminusM.get());
      //std::cout<< "***************************PhaseSpaceDecay:: add stable product XNminusM "<<XNminusM->Pdg()<<std::endl;
    }
    if(dynamic_cast<DecayingParticle*>(p)) {
      dproducts.push_back(*dynamic_cast<DecayingParticle*>(p));
      //std::cout<< "***************************PhaseSpaceDecay:: add decaying product p "<<p->Pdg()<<" "<<dproducts.size()<<" ptr "<<std::endl;
    }
    else{
      pproducts.push_back(*p);
      //std::cout<< "***************************PhaseSpaceDecay:: add stable product p "<<p->Pdg()<<std::endl;
    }

    //OK sorted all the products now make final 2-body state
    auto XNminusMtest =  DecayingParticle{elSpectro::particles::inter_phasespace,
					  CloneModel(PhaseSpaceDecay{dproducts,pproducts} ),
					  CloneDecayer(TwoBodyFlat())};
    XNminusMtest.SetMassDist(cpp::MakeBaseShared<Distribution>( massMaster ) );
    //get the ptr to the new master mass distributon
    //the particle will be the owner of the distribution, not me.
    _massMaster = dynamic_cast< DistFlatMassMaster* >(XNminusMtest.MassDistribution());
    _massMaster->SetParentPtr(parent);

    
    p = ps[N-1]; //i.e. pN 

    //make parent two-body (final) 
    ResetProducts({&XNminusMtest,p});
    auto final2body=dynamic_cast<DecayingParticle*>(MutableProduct(0));

    //will only have 1 channel
     final2body->Channels().SetCurrChannel(0);
     //add final2body products to ptrs
     //std::cout<<"PhaseSpaceDecay final ptrs "<<final2body->Model()->MutableProduct(0)->Pdg()<<" "<<final2body->Model()->MutableProduct(1)->Pdg()<<" "<<MutableProduct(1)->Pdg()<<std::endl;
     if(final2body->Model()->MutableProduct(0)->Pdg()!=elSpectro::particles::inter_phasespace) ptrs.push_back(final2body->Model()->MutableProduct(0));
     if(final2body->Model()->MutableProduct(1)->Pdg()!=elSpectro::particles::inter_phasespace) ptrs.push_back(final2body->Model()->MutableProduct(1));
     ptrs.push_back(MutableProduct(1)); //the 2nd product, pN

  
     // std::cout<<"PhaseSpaceDecay Print pointers for "<<parent->Pdg()<<" = ";
     // for(const auto* pt:ptrs){
     //   std::cout<<" "<<pt->Pdg()<<" ";
     //   if(pt->Pdg()==0)exit(0);
     // }
     // if(parent->Pdg()==0)exit(0);
     // std::cout<<std::endl;

    //update mass distribution ptrs
    _massMaster->SetParticlePtrs(ptrs);
    //    std::cout<<"\n\n\n*************************************"<<std::endl;
        
    // XNminusM will go out of scope, but ResetProducts made a copy
  }
  
  void PhaseSpaceDecay::Print() const{
    DecayModel::Print();
    std::cout<<"PhaseSpaceDecay Print "<<std::endl;
    std::vector<double> masses;
    GetStableMasses(masses);
    for(auto& mass:masses){
      std::cout<<" stable mass "<<mass;
    }
    std::cout<<std::endl;
    if(_massMaster)_massMaster->Print();
    else{
       return;
    }
 
    auto& unstable = UnstableProducts() ;
    for(auto& un:unstable){
       if(un.MassDistribution())un.MassDistribution()->Print();
       un.Model()->Print();
    
    }
  }
  void PhaseSpaceDecay::SetMassMaster(DistFlatMassMaster* master){
    // _massMaster=master;
    auto& unstable = UnstableProducts() ;
    for(auto& un:unstable){
      if(dynamic_cast<DistFlatMass*>(un.MassDistribution())){
	dynamic_cast<DistFlatMass*>(un.MassDistribution())->SetMaster(master);
	if(dynamic_cast<PhaseSpaceDecay*>(un.Model()) ){
	  dynamic_cast<PhaseSpaceDecay*>(un.Model())->SetMassMaster(master);
	}
      }
    }
  }
  
}



