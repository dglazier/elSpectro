namespace elSpectro{};
using namespace elSpectro;

namespace elSpectro{namespace escat{}};
using namespace escat;

namespace jpacPhoto{};
using namespace jpacPhoto;



void Load(){

 
 
  TString JPAC = gSystem->Getenv("JPACPHOTO");
  
  TString ELSPECTRO = gSystem->Getenv("ELSPECTRO");

  //if not defined jpacPhoto look for elSpectro submodule
  if(JPAC.Length()==0) JPAC = ELSPECTRO+"jpacPhoto";
  
  gInterpreter->AddIncludePath(JPAC+"/include/");
  gSystem->Load(JPAC+"/lib/libjpacPhoto."+gSystem->GetSoExt());

  gSystem->Load(ELSPECTRO+"/lib/libelSpectro."+gSystem->GetSoExt());
  gInterpreter->AddIncludePath(ELSPECTRO+"/core");


  gROOT->ProcessLine("elSpectro::Manager::Instance();");

}
