#include "Interface.h"
#include "LorentzVector.h"
#include "DistTF1.h"
#include "FunctionsForElectronScattering.h"
#include "LundWriter.h"
#include "TwoBody_stu.h"


#include <TBenchmark.h>
#include <TH1.h>
#include <TH2.h>
#include <TFile.h>

// ---------------------------------------------------------------------------
// Diagnostic histograms
// ---------------------------------------------------------------------------

double minMass = 0.2;
double maxMass = 3;
TH1F hW("W","W",1000,1,5);
TH1F ht("t","t",1000,0,10);
TH1F hgE("gE","gE",1000,0,20);
TH1F hpP("pP","proton momentum",100,0,10);
TH1F hpTheta("pTh","proton #theta",100,0,180);
TH1F hnP("nP","neutron momentum",100,0,10);
TH1F hnTheta("nTh","neutron #theta",100,0,180);
TH1F heP("eP","e- momentum",100,0,10);
TH1F heTheta("eTh","pe #theta",100,0,180);
TH1F hpoP("poP","proton momentum",100,0,10);
TH1F hpoTheta("poTh","proton #theta",100,0,180);

void PhotoD_npDiLepton(double ebeamE=12,int nEvents = 10000) {

  //  using namespace elSpectro;
  //  elSpectro::Manager::Instance();
  
  auto bremPhoton = initial(22,0,11,
			    model(new Bremsstrahlung()),
			    new BremstrPhoton(ebeamE,0.1*ebeamE,ebeamE*0.999));
  //get 4-momentum of brems. photon
  auto photonP4=bremPhoton->GetInteracting4Vector();
  
  //deuteron target at rest
  auto dTarget= initial(45,0);
  auto dP4=dTarget->GetInteracting4Vector();
  
  auto npee = dynamic_cast<DecayModelDnpee*>(model(new DecayModelDnpee{}));
  auto production = photoprod( bremPhoton,dTarget,(new DecayModelW{0, npee}) );
 

  auto proton  = npee->GetProton();
  auto neutron  = npee->GetNeutron();
  auto ele  = npee->GetElectron();
  auto pos  = npee->GetPositron();
  

  
  // ---------------------------------------------------------------------------
  // Initialize LUND with 10000 events per file
  // ---------------------------------------------------------------------------
  
  //  writer(new LundWriter{Form("out_photo/gd_to_npDiLept_%d.dat",(int)ebeamE),10000});
  
  //initilase the generator, may take some time for making distribution tables 
  initGenerator();
  
  // ---------------------------------------------------------------------------
  // Generate events
  // ---------------------------------------------------------------------------
  
  gBenchmark->Start("e");

  for(int i=0;i<nEvents;i++){
    // for(int i=0;i<0;i++){
    if(i%1000==0) std::cout<<"event number "<<i<<std::endl;
    nextEvent();
    auto& parts = Manager::Instance().Particles().StableParticles();

    for(auto* pa : parts){
      cout<<"     particle "<<pa->Pdg()<<" "<<pa->P4().P()<<" "<<pa->Mass()<<endl;
    }
    cout<<" in "<<*photonP4+*dP4<<" out "<<proton->P4()+neutron->P4()+ele->P4()+pos->P4()<<endl;

    double W = (*photonP4+*dP4).M();
    
    //fill diagnostic histograms
    hW.Fill(W);
    hgE.Fill(photonP4->E());

    hpP.Fill(proton->P4().P());
    hpTheta.Fill(proton->P4().Theta()*TMath::RadToDeg());
    hnP.Fill(neutron->P4().P());
    hnTheta.Fill(neutron->P4().Theta()*TMath::RadToDeg());
    heP.Fill(ele->P4().P());
    heTheta.Fill(ele->P4().Theta()*TMath::RadToDeg());
    hpoP.Fill(pos->P4().P());
    hpoTheta.Fill(pos->P4().Theta()*TMath::RadToDeg());

 
  }
  gBenchmark->Stop("e");
  gBenchmark->Print("e");
 
  generator().Summary();
  // internally stored histograms for total ep cross section
  TH1D *hWdist = (TH1D*)gDirectory->FindObject("Wdist");
  TH1D *hGenWdist = (TH1D*)gDirectory->FindObject("genWdist");

  TFile *fout = TFile::Open(Form("out_Dnpee/eD_to_npee_%d.root",(int)ebeamE), "recreate");
  // total ep cross section inputs
  if(hWdist)hWdist->Write();
  if(hGenWdist)hGenWdist->Write();
 
  // generated event distributions
  hW.Write();
  hgE.Write();
  hpP.Write();
  hpTheta.Write();
  hnP.Write();
  hnTheta.Write();
  heP.Write();
  heTheta.Write();
  hpoP.Write();
  hpoTheta.Write();
  fout->Close();
 

  
}
