namespace elSpectro{};
using namespace elSpectro;

namespace elSpectro{namespace escat{}};
using namespace escat;

namespace jpacPhoto{};
using namespace jpacPhoto;



void Load(){
  gSystem->Load("libEG");
 
  TString JPAC = gSystem->Getenv("JPACPHOTO");
  TString ELSPECTRO = gSystem->Getenv("ELSPECTRO");
  TString HEPMC3 = gSystem->Getenv("HEPMC3");
  
  gInterpreter->AddIncludePath(JPAC+"/include/");
  gSystem->Load(JPAC+"/lib/libjpacPhoto."+gSystem->GetSoExt());

  gSystem->Load(ELSPECTRO+"/lib/libelSpectro."+gSystem->GetSoExt());
  gInterpreter->AddIncludePath(ELSPECTRO+"/core");

  if(HEPMC3.Length()){
    gSystem->AddDynamicPath(HEPMC3+"/hepmc3-install/lib64/:");
    gSystem->AddDynamicPath(HEPMC3+"/hepmc3-install/lib/:");//local install
    
    gInterpreter->AddIncludePath(HEPMC3+"/hepmc3-install/include/");
    gSystem->Load(TString("libHepMC3.")+gSystem->GetSoExt()); //local install
  
  }
  
  ROOT::Math::LorentzRotation ddd;

  // gROOT->ProcessLine("#include \"Interface.h\"");
  //gROOT->ProcessLine("#include \"Manager.h\"");
  gROOT->ProcessLine("elSpectro::Manager::Instance();");

}
