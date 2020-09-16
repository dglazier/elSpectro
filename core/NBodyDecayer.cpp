
#include "NBodyDecay.h"
#include "Interface.h"


namespace elSpectro{

  NBodyDecayer::NBodyDecayer(DecayParticle* parent, particle_ptrs ps, const std::vector<int> pdgs)
  {
    auto particleMan = particles();

    
    //Make N-1 particles decaying into 2 bodies
    //Provide mass distribution

    //create and register all final particles
    for(const auto& pdg: pdgs)
      ps.push_front( particleMan.Take( new Particle{pdg}) );

    if(ps.size() < 3){
      auto pdg2=RegisterNewParticle(0);
      auto X2=particle(pdg2,model(new PhaseSpaceDecay{ps,{}}));
      return X2
    }

   
    auto pdgXNm1=RegisterNewParticle(0,new DistFlatMassMaster(parent,ps));
    auto massMaster  = particleMan().GetMassDist(pdgXNm1);
    
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
    for(uint m=1;m<N-1;++m){
      auto p = ps[m];
      pdg=RegisterNewParticle(0,new DistFlatMass(massMaster));
      XNminusM = particleMan.Take(new DecayingParticle{pdg,model(new PhaseSpaceDecay(XNminusM,p)),new TwoBodyFlat()});
    }
    //give master to the X(N-1) state
    auto p = ps[N-1]; //m=1
    XNminusM = particleMan.Take(new DecayingParticle{pdgXNm1,model(new PhaseSpaceDecay(XN,p)),new TwoBodyFlat()});

    return XNminusM;
  }
}


