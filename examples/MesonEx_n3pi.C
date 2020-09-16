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
  
  
  mass_distribution(9995,new DistTF1{TF1("hh","1",0.,4)});
  //mass_distribution(9995,new DistTF1{TF1("hh","TMath::BreitWigner(x,4.22,0.05)",3.9,4.5)});
  auto X=static_cast<DecayingParticle*>(particle(9995,model(new PhaseSpaceDecay{{},{211,211,-211}})) );
  std::cout<<" X PRODUCTS "<<X->Model()->Products().size()<<endl;
  //decay of pGamma* tp X + n
  auto pGammaStarDecay = model(new DecayModelst{ {X},{2112} });
  
  // 
  //create mesonex electroproduction of X + neutron
  auto production=mesonex( ebeamE ,  new DecayModelQ2W{0, pGammaStarDecay });
  // generator().SetModelForMassPhaseSpace(pGammaStarDecay);
  
  
  elSpectro::LorentzVector elbeam(0,0,ebeamE,elSpectro::escat::E_el(ebeamE));
  elSpectro::LorentzVector prbeam(0,0,0,elSpectro::escat::M_pr());
  
  auto pip1 = static_cast<DecayingParticle*>(X)->Model()->Products()[1];
  auto X2 = static_cast<DecayingParticle*>((X)->Model()->Products()[0]);
  std::cout<<" X2 PRODUCTS "<<X2->Model()->Products().size()<<endl;

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
    if(i%100==0) std::cout<<"event number "<<i<<std::endl;
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
 
  // compute total ep cross section from internally stored histograms
  if(hGenWdist&&hWdist){
   if(hGenWdist->Integral()){
     hGenWdist->Scale(1.0/hGenWdist->Integral()); // normalize flux PDF
     double integrated_xsection = 0; // get sigma_ep from integral over W: f(W)*sigma_gp(W)
     for(int i=0; i<hWdist->GetNbinsX(); i++) {
       double W = hWdist->GetXaxis()->GetBinCenter(i+1);
       double WbinWidthScale = hWdist->GetBinWidth(i+1);
       double W_xsection = hWdist->GetBinContent(i+1);
    
       int genWbin = hGenWdist->GetXaxis()->FindBin(W);
       WbinWidthScale /= hGenWdist->GetXaxis()->GetBinWidth(genWbin);
       double W_fluxWeight = hGenWdist->GetBinContent(genWbin) * WbinWidthScale;
       integrated_xsection += W_xsection * W_fluxWeight;
     }
     cout<<"Total cross section ("<<(int)ebeamE<<","<<"): sigma_ep = "<<integrated_xsection<<" nb "<<endl;
   }
  }
}
