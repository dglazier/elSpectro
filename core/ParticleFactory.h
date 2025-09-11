
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
///            Example, print data for a particle particelFactory.GetData("name").Print();

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
#include "TSystem.h"


namespace elSpectro{

  namespace particles{

    constexpr short inter_phasespace=-9;
    constexpr short gamma_star=-22;
    constexpr short gamma_star_nucleon=-2211;
 
    class ParticleFactory {

    public:
    
      static ParticleFactory& Instance() { static ParticleFactory instance; return instance; }
 
      static void  Init(){
	//	TDatabasePDG *pdgDB = new TDatabasePDG();
	auto pdgDB = TDatabasePDG::Instance();


	
        pdgDB->ReadPDGTable(Form("%s/etc/xyz_pdg_table.txt",gSystem->Getenv("ELSPECTRO")));
        pdgDB->ReadPDGTable(Form("%s/etc/rndmflav_pdg_table.txt",gSystem->Getenv("ELSPECTRO")));
	pdgDB->ReadPDGTable(Form("%s/etc/nuclei_table.txt",gSystem->Getenv("ELSPECTRO")));
	//	pdgDB->ReadPDGTable(Form("%s/etc/pdg_table.txt",gSystem->Getenv("ROOTSYS")));
	pdgDB->ReadPDGTable(Form("%s/etc/pdg_table.txt",gSystem->Getenv("ELSPECTRO")));

	//Must read table via files first,
	//as AddPArticle generates defualt table
  	//special elspectro particles first
	// //name,title,mass,stable,width,charge,type.code 
	pdgDB->AddParticle("intermediate_phasespace","intermediate_phasespace",
			   0.0, kFALSE, 0, 0, "virtual", inter_phasespace);
	pdgDB->AddParticle("gamma_star","gamma_star", 0.0, kFALSE, 0, 0,
			   "virtual", gamma_star);
	pdgDB->AddParticle("gamma_star_nucleon","gamma_star_nucleon",
	 		   pdgDB->GetParticle("proton")->Mass(), kFALSE,
	 		   0, 0, "virtual",gamma_star_nucleon );
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
