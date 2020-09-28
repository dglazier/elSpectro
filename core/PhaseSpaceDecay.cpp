#include "PhaseSpaceDecay.h"
#include "Interface.h"
#include "ParticleManager.h"
#include "DistFlatMass.h"


namespace elSpectro{

  PhaseSpaceDecay::PhaseSpaceDecay( particle_ptrs ps, const std::vector<int> pdgs):
    DecayModel{ps,pdgs}
  {

    _name={"PhaseSpaceDecay"};

  }

  void PhaseSpaceDecay::SetParent(DecayingParticle* pa){
    DecayModel::SetParent(pa);
    //now can make cascade of decays
    if( Products().size()>2 )
      nBodyDecayer(pa,StableProducts(),UnstableProducts());
 
  }
  void PhaseSpaceDecay::PostInit(ReactionInfo* info){
    DecayModel::PostInit(info);
    if(Parent()->MassDistribution()==nullptr){//production proces does not have mass distribution
      if(dynamic_cast<ProductionProcess*>(Parent())==nullptr){
	std::cerr<<"PhaseSpaceDecay::PostInit parent needs a mass distribution for pdg = "<<Parent()->Pdg();
	std::cerr<<"\n  you need to use \n >>   mass_distribution(PDG,new DistTF1{TF1(\"massDist\",\"1\",MINMASS,MAXMASS)});";
	std::cerr<<" \n where PDG (9995-9999) is the pdg number you assigned the decaying particle, and MINMAMSS and MAXMASS is the mass limits it will be allowed to have, for pure phase space this must be at least the kinematically allowed range";
	std::cerr<<"\n NOTE eventually this will be automated! "<<std::endl;
	
	exit(0);
      }
    }
  }
  
  void PhaseSpaceDecay::nBodyDecayer(DecayingParticle* parent, const particle_ptrs stable,  const decaying_ptrs unstable ) //take copies of particle vectors
  {
    std::cout<<"Start Particle* nBodyDecayer "<<parent<<std::endl;
    auto& particleMan = Manager::Instance().Particles();
    
    
    //Make N-1 particles decaying into 2 bodies
    //Provide mass distribution

    //create and register all final particles
    //add unstable at start so more room for selecting their mass
    particle_ptrs ps;
    for(const auto& unp: unstable)
      ps.push_back( unp );
    for(const auto& sp: stable)
      ps.push_back( sp );

 
   std::cout<<"Particle* nBodyDecayer "<<ps.size()<<std::endl;
 
    //Register our X-1 particle type
    //We will construct this particle at the end of the rest
    auto pdgXNm1=particleMan.RegisterNewPdgParticle(0,new DistFlatMassMaster(parent,ps));
      std::cout<<"Particle* nBodyDecayer pdg "<< pdgXNm1 <<std::endl;
 
    auto massMaster  = static_cast<DistFlatMassMaster*>(particleMan.GetMassDist(pdgXNm1));
    std::cout<<"Particle* nBodyDecayer pdg "<< massMaster <<std::endl;

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
    
    Particle* XNminusM=ps[0];
    auto N=ps.size();
    for(uint m=1;m<=N-3;++m){//leave 2 particles
      auto p = ps[m];
      auto pdg=particleMan.RegisterNewPdgParticle(0,new DistFlatMass(massMaster));
      XNminusM = particleMan.Take(new DecayingParticle{pdg,model(new PhaseSpaceDecay({XNminusM,p},{})),new TwoBodyFlat()});
    }
    //give master to the X(N-1) state
    auto p = ps[N-2]; //m=1
    //combine with the X(N-2) particle and give pdgXNm1 mass distribution
    XNminusM = particleMan.Take(new DecayingParticle{pdgXNm1,model( new PhaseSpaceDecay{{XNminusM,p},{}} ), new TwoBodyFlat() } );

    p = ps[N-1]; //i.e. pN 
    //make parent two-body
    ResetProducts({XNminusM,p});
  
  }
  
  
}


