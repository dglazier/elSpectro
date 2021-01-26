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
  
  gInterpreter->AddIncludePath(JPAC+"/include/");
  gSystem->Load(JPAC+"/lib/libjpacPhoto."+gSystem->GetSoExt());

  gSystem->Load(ELSPECTRO+"/lib/libelSpectro."+gSystem->GetSoExt());
  gInterpreter->AddIncludePath(ELSPECTRO+"/core");


  ROOT::Math::LorentzRotation ddd;
  RooFunctorPdfBinding forLinkingatJLAB; //not sure why, but will not load without this on ifarm...


  // gROOT->ProcessLine("#include \"Interface.h\"");
  //gROOT->ProcessLine("#include \"Manager.h\"");
  gROOT->ProcessLine("elSpectro::Manager::Instance();");

}
