{
  gSystem->Load("libEG");
 
  TString JPAC = gSystem->Getenv("JPACPHOTO");
  TString ELSPECTRO = gSystem->Getenv("ELSPECTRO");
  TString HEPMC3 = gSystem->Getenv("HEPMC3");
  
  gInterpreter->AddIncludePath(JPAC+"/include/");
  gSystem->Load(JPAC+"/build/libJpacPhoto."+gSystem->GetSoExt());

  gSystem->Load(ELSPECTRO+"/lib/libelSpectro."+gSystem->GetSoExt());
  gInterpreter->AddIncludePath(ELSPECTRO+"/core");

  gInterpreter->AddIncludePath(HEPMC3+"/hepmc3-install/include/");
  gSystem->Load(HEPMC3+"/hepmc3-install/lib64/libHepMC3."+gSystem->GetSoExt());

  ROOT::Math::LorentzRotation ddd;
 
}
