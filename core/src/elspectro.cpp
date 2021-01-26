#include <iostream>
#include <TRint.h>
#include <TEnv.h>
#include <TString.h>
#include <TSystem.h>


int main(int argc, char **argv) {


  //get command line options first check if makeall
  TString macroName;
  for(Int_t i=1;i<argc;i++){
    TString opt=argv[i];
    if((opt.Contains(".C"))) macroName=opt;
  }
  
  TRint  *app = new TRint("elSpectro", &argc, argv);
  // Run the TApplication (not needed if you only want to store the histograms.)
  app->ProcessLine(".x $ELSPECTRO/core/src/Load.C");

  app->ProcessLine(Form(".x %s",macroName.Data()));
 

  app->Terminate(0);
  
  return 0;
  
}
