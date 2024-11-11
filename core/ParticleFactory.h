
//////////////////////////////////////////////////////////////
///
///Class:		ParticleFactory
///Description:
///            Class to manage particles
///            Should only be accessed via particles::Factory::Instance()
///   
///            Construct ParticleData objects for given pdg
///            elSpectro::Particles for generator are then created
///            from ParticleData
///            Relevent DecayModels are responsible for ownership/keeping alive
#pragma once



//Need to automate this from downloadable PDG data
//https://pdg.lbl.gov/2022/html/computer_read.html
//which gives properties and decays and
//total photoproduction cross section
//**Note this turned out to not be viable at the moment
//** branching decays are not well defined


#include "ParticleData.h"
#include "DecayingParticle.h"
#include "CppHelperFuncs.h"



namespace elSpectro{

  namespace particles{

 
    class ParticleFactory {

    public:
    
      static ParticleFactory& Instance() { static ParticleFactory instance; return instance; }
 
      static void  Init(){
	TDatabasePDG *pdgDB = new TDatabasePDG();
	//	pdgDB->ReadPDGTable(Form("%s/etc/el_pdg_table.txt",gSystem->Getenv("ELSPECTRO")));
	
	//name,title,mass,stable,width,charge,type.code 
	pdgDB->AddParticle("gamma_star","gamma_star", 0.0, kFALSE,
			   0, 0, "virtual", -22);
	pdgDB->AddParticle("gamma_star_nucleon","gamma_star_nucleon",
			   pdgDB->GetParticle("proton")->Mass(), kFALSE,
			   0, 0, "virtual", -2211);
      }
      
      ParticleData  GetData(int pdg){
	return GetData(TDatabasePDG::Instance()->GetParticle(pdg)->GetName());
      }

      ParticleData  GetData(const std::string& name);
    
      DecayingParticle CreateDecayingParticle(int pdg);
      DecayingParticle CreateDecayingParticle(const std::string& name);
    
    private:
      std::vector<ParticleData> _existing;
      std::vector<std::string> _names;
    };
  
  }
}
