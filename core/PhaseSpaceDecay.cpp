#include "PhaseSpaceDecay.h"
#include "Interface.h"
#include "ParticleFactory.h"
#include "DistFlatMass.h"


namespace elSpectro{

  //  PhaseSpaceDecay::PhaseSpaceDecay( particle_ptrs ps, const std::vector<int> pdgs):
  PhaseSpaceDecay::PhaseSpaceDecay(const decaying_objs& decs, const particle_objs& stables):
    DecayModel{decs,stables}
  {
    std::cout<<"PhaseSpaceDecay::PhaseSpaceDecay, stable =  "<<StableProducts().size()<<" unstable = "<<UnstableProducts().size()<<" all = "<<Products().size()<<std::endl;
    _name={"PhaseSpaceDecay"};

  }

  void PhaseSpaceDecay::SetParentAndProducts(DecayingParticle* pa,  particle_objs& stable, decaying_objs& unstable){
    DecayModel::SetParent(pa);
    nBodyDecayer(pa,stable,unstable);
  }
  
  void PhaseSpaceDecay::SetParent(DecayingParticle* pa){
    DecayModel::SetParent(pa);
    //now can make cascade of decays
    if( Products().size()>2 )
      nBodyDecayer(pa,MutableStableProducts(),MutableUnstableProducts());
 
  }
  void PhaseSpaceDecay::PostInit(ReactionInfo* info){
    std::cout<<"PhaseSpaceDecay::PostInit " <<Parent()<<" "<<info<<" "<<Parent()->MaximumMassPossible()<<" pdg "<<Parent()->Pdg()<<std::endl;
    auto maxmass=Parent()->MaximumMassPossible();
    if(maxmass==0) maxmass = info->Wmax(); //if parent has no max just use reaction W
    auto flatmass=dynamic_cast<DistFlatMass*>(Parent()->MassDistribution());
    if(flatmass) flatmass->SetMaxX(maxmass);
  
    for(auto& product:Products()){
      auto flatmass=dynamic_cast<DistFlatMass*>(product->MassDistribution());
      SumAllProducts();
      std::cout<<"              PhaseSpaceDecay::PostInit product " <<product->Pdg()<<" flat mass dist ? "<<flatmass<<" "<<SumOfProductMasses()<<std::endl;
      std::cout<<"              PhaseSpaceDecay::PostInit set mass to  " <<maxmass<<std::endl;
      if(flatmass) flatmass->SetMaxX(maxmass);
    }
    
    if(_massMaster!=nullptr)_massMaster->SetMaxX(maxmass);
    std::cout<<"PhaseSpaceDecay::PostInit " <<Parent()<<" "<<info<<" "<<maxmass<<std::endl;
 
    DecayModel::PostInit(info);

    if(Parent()->MassDistribution()==nullptr){//production proces does not have mass distribution
      if(dynamic_cast<ProductionProcess*>(Parent())==nullptr&&Parent()->Pdg()!=-2211){
	std::cerr<<"PhaseSpaceDecay::PostInit parent needs a mass distribution for pdg = "<<Parent()->Pdg();
	std::cerr<<"\n  you need to use \n >>   mass_distribution(PDG,new DistTF1{TF1(\"massDist\",\"1\",MINMASS,MAXMASS)});";
	std::cerr<<" \n where PDG (9995-9999) is the pdg number you assigned the decaying particle, and MINMAMSS and MAXMASS is the mass limits it will be allowed to have, for pure phase space this must be at least the kinematically allowed range";
	std::cerr<<"\n NOTE eventually this will be automated! "<<std::endl;
	
	//	exit(0);
      }
    }

  }
  
  void PhaseSpaceDecay::nBodyDecayer(DecayingParticle* parent,particle_objs& stable,decaying_objs& unstable ) //take copies of particle vectors
  {
    std::cout<<"Start Particle* nBodyDecayer "<<parent<<std::endl;
    
    
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
 
 
   std::cout<<"Particle* nBodyDecayer "<<ps.size()<<std::endl;
 
    //Register our X-1 particle type
    //We will construct this particle at the end of the rest
   /* auto pdgXNm1=particleMan.RegisterNewPdgParticle(0,new DistFlatMassMaster(parent,ps));
    std::cout<<"Particle* nBodyDecayer pdg "<< pdgXNm1 <<std::endl;
 
    _massMaster  = dynamic_cast<DistFlatMassMaster*>(particleMan.GetMassDist(pdgXNm1));
   */
   auto massMaster = DistFlatMassMaster{parent,ps};
   std::cout<<"Particle* nBodyDecayer pdg "<< _massMaster <<std::endl;

    //if(ps.size() < 3){//in case ony two particle just use 2 body decay
    // return ;
    //}
 
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
    
   std::unique_ptr<Particle> XNminusM(ps[0]);
    auto N=ps.size();
    for(uint m=1;m<=N-3;++m){//leave 2 particles
      auto p = ps[m];
      //auto pdg=particleMan.RegisterNewPdgParticle(0,new DistFlatMass(_massMaster));
      std::cout<<"PhaseSpaceDecay:: create new two-body "<<" from "<<p->Pdg()<<" and "<<XNminusM->Pdg()<<std::endl;
      //PhaseSpaceDecay will copy the particles we give
      //so need to get ptrs to them for mass distribution after
      particle_objs pproducts;
      decaying_objs dproducts;
      //put copies into products vectors
      if(dynamic_cast<DecayingParticle*>(XNminusM.get())) dproducts.push_back(*dynamic_cast<DecayingParticle*>(XNminusM.get()));
      else pproducts.push_back(*XNminusM.get());
      if(dynamic_cast<DecayingParticle*>(p)) dproducts.push_back(*dynamic_cast<DecayingParticle*>(p));
      else pproducts.push_back(*p);
     
      /*XNminusM = particleMan.Take(new DecayingParticle{pdg,
						       CloneModel(PhaseSpaceDecay{dproducts,pproducts} ),
						       CloneDecayer(TwoBodyFlat())} );*/
      XNminusM.reset(new DecayingParticle{0,
					  CloneModel(PhaseSpaceDecay{dproducts,pproducts} ),
					  CloneDecayer(TwoBodyFlat())});
      XNminusM->SetMassDist(cpp::MakeBaseShared<Distribution>(DistFlatMass(&massMaster)));
      
      //  if(m==1) ptrs.push_back(XNminusM->Model()->MutableProduct(0)); //the 1st product, only first time
      ptrs.push_back(dynamic_cast<DecayingParticle*>(XNminusM.get())->Model()->MutableProduct(1)); //the 2nd product, as given 2nd above
    }

    
    //give master to the X(N-1) state
    auto p = ps[N-2]; //m=1
    //combine with the X(N-2) particle and give pdgXNm1 mass distribution
    std::cout<<"PhaseSpaceDecay:: final  create new two-body "<<" from "<<p->Pdg()<<" and "<<XNminusM->Pdg()<<std::endl;
    particle_objs pproducts;
    decaying_objs dproducts;
    if(dynamic_cast<DecayingParticle*>(XNminusM.get())) dproducts.push_back(*dynamic_cast<DecayingParticle*>(XNminusM.get()));
    else pproducts.push_back(*XNminusM.get());
    if(dynamic_cast<DecayingParticle*>(p)) dproducts.push_back(*dynamic_cast<DecayingParticle*>(p));
    else pproducts.push_back(*p);
    
    /* XNminusM = particleMan.Take(new DecayingParticle{pdgXNm1,
						     CloneModel(PhaseSpaceDecay{dproducts,pproducts} ),
						     CloneDecayer(TwoBodyFlat())} );*/
    // XNminusM.reset(new DecayingParticle{0,
    // 					CloneModel(PhaseSpaceDecay{dproducts,pproducts} ),
    // 					CloneDecayer(TwoBodyFlat())});

    auto XNminusMtest =  DecayingParticle{0,
      CloneModel(PhaseSpaceDecay{dproducts,pproducts} ),
      CloneDecayer(TwoBodyFlat())};
    XNminusMtest.SetMassDist(cpp::MakeBaseShared<Distribution>( massMaster ) );
    //get the ptr to the new master mass distributon
    //the particle will be the owner of the distribution, not me.
    _massMaster = dynamic_cast< DistFlatMassMaster* >(XNminusMtest.MassDistribution());
     std::cout<<"PhaseSpaceDecay:: test " <<XNminusMtest.Model()<<std::endl;
   
    //    ptrs.push_back(dynamic_cast<DecayingParticle*>(XNminusM.get())->Model()->MutableProduct(1)); //the 2nd product, as given 2nd above
   
    p = ps[N-1]; //i.e. pN 
    std::cout<<"PhaseSpaceDecay:: and final particle products "<<XNminusMtest.Pdg()<<" "<<p->Pdg()<<std::endl;
    //make parent two-body
    //    ResetProducts({XNminusMtest.get(),p});
    ResetProducts({&XNminusMtest,p});
    auto final2body=dynamic_cast<DecayingParticle*>(MutableProduct(0));
    
    std::cout<<"PhaseSpaceDecay:: " <<final2body<<std::endl;
    std::cout<<"PhaseSpaceDecay:: " <<final2body->Model()<<std::endl;
    std::cout<<"PhaseSpaceDecay:: chans " <<final2body->Channels().N()<<" "<<final2body->Channels().CurrChannel()<<std::endl;
    final2body->Channels().SetChannel(0);
    std::cout<<"PhaseSpaceDecay:: " <<final2body->Model()->Products().size()<<" "<<final2body->Model()->StableProducts().size()<<" "<<final2body->Model()->UnstableProducts().size()<<std::endl;
    ptrs.push_back(final2body->Model()->MutableProduct(0)); //from the 1st  product, XNminus1
    ptrs.push_back(final2body->Model()->MutableProduct(1)); //from the 1st  product, XNminus1
    ptrs.push_back(MutableProduct(1)); //the 2nd product, pN

    if(ptrs.size()!=N){
      std::cerr<<"Error : PhaseSpaceDecay nBody made more particles than we had"<<ptrs.size()<<" "<<N<<std::endl;
      exit(0);
    }
    std::cout<<"Print pointers ";
    for(const auto* pt:ptrs){
      std::cout<<pt<<" ";
    }
    std::cout<<std::endl;
    _massMaster->SetParticlePtrs(ptrs);
    // XNminusM will go out of scope, but ResetProducts made a copy
  }
  
  
}


