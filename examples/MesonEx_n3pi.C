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
TH1F hQ2("Q2","Q2",1000,0,5);
TH1D heE("eE","eE",1000,0,20);
TH1D heTh("eTh","eTh",1000,0,180);
TH1D hYTh("YTh","YTh",1000,0,180);
TH1F hW("W","W",1000,1,4);
TH1F ht("t","t",1000,0,10);
TH1F hgE("gE","gE",1000,0,20);
TH1F hgTh("gTh","gTh",1000,0,180);
TH1F hXM("MesonM","; X to #pi#pi#pi Mass (GeV)",1000,minMass,maxMass);
TH1F hM12("MesonM12",";  #pi+#pi+ Mass (GeV)",1000,minMass,maxMass);
TH1F hMm2("MesonMm2",";  #pi+#pi- Mass (GeV)",1000,minMass,maxMass);
TH1F hM1m("MesonM1m",";  #pi-#pi+ Mass (GeV)",1000,minMass,maxMass);
TH1F hP1("P1","Pi1",100,0,10);
TH1F hP2("P2","Pi2",100,0,10);
TH1F hPm("Pm","Pi-",100,0,10);

void MesonEx_n3pi(double ebeamE=10.4,int nEvents = 5e4) {

  using namespace elSpectro;
  elSpectro::Manager::Instance();
  
  elSpectro::LorentzVector elbeam(0,0,ebeamE,elSpectro::escat::E_el(ebeamE));
  elSpectro::LorentzVector prbeam(0,0,0,elSpectro::escat::M_pr());
 

  //flat pure phase sapce distribution =1 for particle id 9995
  mass_distribution(9995,new DistTF1{TF1("hh","1",0.,(elbeam+prbeam).M())});
  //add a Breit-Wigner resonance for particle id 9995
  //mass_distribution(9995,new DistTF1{TF1("hh","TMath::BreitWigner(x,1.5,0.1)",3.9,4.5)});
  auto X=static_cast<DecayingParticle*>(particle(9995,model(new PhaseSpaceDecay{{},{211,211,-211}})) );
  
  //decay of pGamma* tp X + n
  auto pGammaStarDecay = model(new DecayModelst{ {X},{2112} });
  
  // 
  //create mesonex electroproduction of X + neutron
  //TwoBody_stu{0.1, 0.9, 3 ,0,0} //0.1 strength  s distribution (flat angular dist.),  0.9 strength t distribution with slope b = 3
  mesonex( ebeamE ,  new DecayModelQ2W{0, pGammaStarDecay,new TwoBody_stu{0.1, 0.9, 3 , 0 , 0}} );

  //give limits to the detected electron 
  auto production=dynamic_cast<ElectronScattering*>(generator().Reaction());
  production->SetLimitTarRest_eThmin(1.5*TMath::DegToRad());
  production->SetLimitTarRest_eThmax(5.5*TMath::DegToRad());
  production->SetLimitTarRest_ePmin(0.4);
  production->SetLimitTarRest_ePmax(6);

  //get pointers to produced particles fror diagnostic histos
  auto pip1 = static_cast<DecayingParticle*>(X)->Model()->Products()[1];
  auto X2 = static_cast<DecayingParticle*>((X)->Model()->Products()[0]);

  auto pip2 = static_cast<DecayingParticle*>(X2)->Model()->Products()[0];
  auto pim = static_cast<DecayingParticle*>(X2)->Model()->Products()[1];
  
  auto neutron = pGammaStarDecay->Products()[1];
  auto electron = dynamic_cast<DecayModelQ2W*>(production->Model())->GetScatteredElectron();
  
  // ---------------------------------------------------------------------------
  // Initialize LUND
  // ---------------------------------------------------------------------------
  
  writer(new LundWriter{Form("out_mesonex/ep_to_nX3pi_%d.dat",(int)ebeamE)});
  
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

    //fill diagnostic histograms
    auto photon = elbeam - electron->P4();
    double Q2 = -photon.M2();
    double W = (photon+prbeam).M();
    double t = -1*(neutron->P4()-prbeam).M2();
    hQ2.Fill(Q2);
    hW.Fill(W);
    ht.Fill(t);
    hgE.Fill(photon.E());
    
    auto elec = electron->P4();
    heTh.Fill(elec.Theta() *TMath::RadToDeg());
    heE.Fill(elec.E());

    auto Xrec=pip1->P4()+pip2->P4()+pim->P4();
    hXM.Fill(Xrec.M());
    Xrec=pip1->P4()+pip2->P4();
    hM12.Fill(Xrec.M());
    Xrec=pip1->P4()+pim->P4();
    hM1m.Fill(Xrec.M());
    Xrec=pim->P4()+pip2->P4();
    hMm2.Fill(Xrec.M());

    hP1.Fill(pip1->P4().P());
    hP2.Fill(pip2->P4().P());
    hPm.Fill(pim->P4().P());
  }
  gBenchmark->Stop("e");
  gBenchmark->Print("e");
 
 
  // internally stored histograms for total ep cross section
  TH1D *hWdist = (TH1D*)gDirectory->FindObject("Wdist");
  TH1D *hGenWdist = (TH1D*)gDirectory->FindObject("genWdist");

  TFile *fout = TFile::Open(Form("out_mesonex/ep_to_nX3pi_%d_4.root",(int)ebeamE), "recreate");
  // total ep cross section inputs
  if(hWdist)hWdist->Write();
  if(hGenWdist)hGenWdist->Write();
 
  // generated event distributions
  hQ2.Write();
  hW.Write();
  ht.Write();
  hgE.Write();
  heTh.Write();
  heE.Write();
  hXM.Write();
  hM12.Write();
  hMm2.Write();
  hM1m.Write();
  fout->Close();
 
   generator().Summary();

}
