{
  gSystem->Load("libEG");
 
  TString JPAC = gSystem->Getenv("JPACPHOTO");
  TString ELSPECTRO = gSystem->Getenv("ELSPECTRO");
  
  gInterpreter->AddIncludePath(JPAC+"/include/");
  gSystem->Load(JPAC+"/build/libJpacPhoto."+gSystem->GetSoExt());

 
  gSystem->Load(ELSPECTRO+"/lib/libelSpectro."+gSystem->GetSoExt());
  gInterpreter->AddIncludePath(ELSPECTRO+"/core");
  ROOT::Math::LorentzRotation ddd;
 
}
