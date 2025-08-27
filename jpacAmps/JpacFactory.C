void JpacInclude(const TString& inc){
  auto ELSPECTRO = gSystem->Getenv("ELSPECTRO");  
  // gROOT->ProcessLine(Form("#include \"%s/jpacAmps/%s\"",ELSPECTRO,inc.Data()));
  gROOT->ProcessLine(Form(".L %s/jpacAmps/%s+",ELSPECTRO,inc.Data()));

}

void JpacFactory(){
  JpacInclude("Vector_mesons.h");
  JpacInclude("X_mesons.h");
  JpacInclude("Z_mesons.h");
  JpacInclude("Y_mesons.h");
  
}
