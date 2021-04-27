#include "ParticleManager.h"
#include <TDatabasePDG.h>

namespace elSpectro{

  ParticleManager::ParticleManager(){
    
    TDatabasePDG *pdgDB = TDatabasePDG::Instance();
    //name,title,mass,stable,width,charge,type.code 
    pdgDB->AddParticle("gamma_star","gamma_star", 0.0, kFALSE,
		       0, 0, "virtual", -22);
    pdgDB->AddParticle("gamma_star_nucleon","gamma_star_nucleon",
		       pdgDB->GetParticle("proton")->Mass(), kFALSE,
		       0, 0, "virtual", -2211);

    //arbitrary resonaces for decaying
    pdgDB->AddParticle("resonance5","resonance1",
		       0, kFALSE,
		       0, 0, "virtual", 9995);
    pdgDB->AddParticle("resonance6","resonance1",
		       0, kFALSE,
		       0, 0, "virtual", 9996);
    pdgDB->AddParticle("resonance7","resonance1",
		       0, kFALSE,
		       0, 0, "virtual", 9997);
    pdgDB->AddParticle("resonance8","resonance1",
		       0, kFALSE,
		       0, 0, "virtual", 9998);
    pdgDB->AddParticle("resonance9","resonance1",
		       0, kFALSE,
		       0, 0, "virtual", 9999);

    pdgDB->AddParticle("deuteron","deuteron", 1.875612, kTRUE,0, 1, "Baryon", 45); //Jlab CLAS numbering
    pdgDB->AddParticle("deuteron","deuteron", 1.875612, kTRUE,0, 1, "Baryon", 1000010020); //PDG code numbering
    
  }
  Particle*  ParticleManager::Take(Particle* p){
    _particles.push_back(particle_uptr{p});
    //make p useable again
    p =_particles.back().get();
    auto pdg=p->Pdg();
    //assign mass distribution if exists
    //note maybe this should be just for unstable
    if(_massDist.count(pdg)!=0)
      p->SetMassDist(_massDist.at(pdg).get());
     
    auto dp=dynamic_cast<DecayingParticle* >( p );
      
    //distinguish between stable and unstable particles
    if(dp!=nullptr){
      _unstables.push_back(dp);
    }
    else{
      _stables.push_back(p);
    }
    
    
    return p;
  }
  
  void ParticleManager::AddToPdgTable(int pdg,double mass){
    TDatabasePDG *pdgDB = TDatabasePDG::Instance();
    pdgDB->AddParticle(Form("Particle%d",pdg),"resonance",
		       mass, kFALSE,
		       0, 0, "virtual", pdg);
  }
  Double_t ParticleManager::GetMassFor(int pdg){
    TDatabasePDG *pdgDB = TDatabasePDG::Instance();
    if(pdgDB->GetParticle(pdg)==nullptr){
      std::cerr<<"ParticleManager::GetMassFor no particle with pdg code "<<pdg<<std::endl;
    }
    return pdgDB->GetParticle(pdg)->Mass();
  }
  
}
