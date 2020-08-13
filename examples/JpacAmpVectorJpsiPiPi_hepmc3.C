#include "Interface.h"
#include "LorentzVector.h"
#include "DecayGammaN_Test.h"
#include "DistTF1.h"
#include "JpacModelQ2W.h"
#include "FunctionsForElectronScattering.h"
#include "HepMC3Writer.h"
#include "LundWriter.h"
#include "TwoBody_stu.h"

#include "reaction_kinematics.hpp"
#include "regge_trajectory.hpp"
#include "amplitudes/vector_exchange.hpp"
#include "amplitudes/pseudoscalar_exchange.hpp"
#include "amplitudes/pomeron_exchange.hpp"
#include "amplitudes/amplitude_sum.hpp"
#include "amplitudes/baryon_resonance.hpp"

#include <TBenchmark.h>
#include <TH1.h>
#include <TH2.h>
#include <TFile.h>

// ---------------------------------------------------------------------------
// Diagnostic histograms
// ---------------------------------------------------------------------------

double minMass = 3.0;
double maxMass = 4.6;
TH1F hQ2("Q2","Q2",1000,0,100);
TH1D heE("eE","eE",1000,0,20);
TH1D heTh("eTh","eTh",1000,0,180);
TH1F hW("W","W",1000,0,100);
TH1F ht("t","t",1000,0,10);
TH1F hgE("gE","gE",1000,0,20);
TH1F hMesonM("MesonM","; J/#psi#pi#pi Mass (GeV)",1000,minMass,maxMass);
TH1F hJpsiM("JpsiM","; e+e- Mass (GeV)",1000,2.9,3.3);
TH1F hRhoM("RhoM","; #pi+#pi- Mass (GeV)",1000,0.0,1.0);
TH2F hElePVsEta("ElePVsEta","; #eta; p (GeV)",200,-5,5,500,0,50);
TH2F hPosPVsEta("PosPVsEta","; #eta; p (GeV)",200,-5,5,500,0,50);
TH2F hPionpPVsEta("PionpPVsEta","; #eta; p (GeV)",200,-5,5,500,0,50);
TH2F hPionmPVsEta("PionmPVsEta","; #eta; p (GeV)",200,-5,5,500,0,50);
TH2F hRecoilPVsEta("RecoilPVsEta","; #eta; p (GeV)",200,0,10,1000,0,275);
TH2F hRecoilThetaVsP("RecoilThetaVsP","; p (GeV); #theta (mrad)",1000,0,275,200,0,200);
TH1F hRecoilPt("RecoilPt","; p_{T} (GeV)",200,0,5.0);


void JpacAmpVectorJpsiPiPi_hepmc3(double ebeamE = 5, double pbeamE = 41, int nEvents = 5e4) {

  using namespace jpacPhoto;
  using namespace elSpectro;
  elSpectro::Manager::Instance();
  
  // ---------------------------------------------------------------------------
  // AMPLITUDES
  // ---------------------------------------------------------------------------

  // ---------------------------------------------------------------------------
  // PSI(2S)
  // ---------------------------------------------------------------------------

  linear_trajectory alpha(+1, 0.941, 0.364, "pomeron");

  // higher mass kinematics
  reaction_kinematics * ptr2S = new reaction_kinematics(3.686, "#psi(2S)");

  // trajectory is the same but new kinematics for bigger phase-space
  pomeron_exchange pomeron_2S(ptr2S, &alpha, false, "100 x #psi(2S)");

  // Same t-slope but coupling is scaled by 1/4, multiplied by 10 for comparision
  pomeron_2S.set_params({10. * 0.379 / 4., 0.12});

  // ---------------------------------------------------------------------------
  // Y(4220)
  // ---------------------------------------------------------------------------

  // Best fit values from [2] from near threshold data
  linear_trajectory alpha19(+1, .94, 0.36, "pomeron (2019)");

  // Best fit values from [1] from high energy data
  linear_trajectory alpha16(+1, 1.1, 0.11, "pomeron (2016)");

  // Set up kinematics, determined entirely by vector meson mass
  reaction_kinematics * ptr = new reaction_kinematics(4.220, "Y(4220)");

  // Pomeron exchange fit to near threshold
  // helicity nonconserving
  pomeron_exchange Y(ptr, &alpha19, false, "2019 fit");
  Y.set_params({1.54 * .379, .12});

  // Pomeron exchange for higher energies
  // helicity conserving + alpha intercept > 1
  pomeron_exchange Y2(ptr, &alpha16, true, "2016 fit");
  Y2.set_params({1.54 * .159, 1.01});

  // Sum individual contributions together incoherently
  // amplitude_sum sum(ptr, {&Y, &Y2}, "Sum");
  amplitude_sum sum(ptr, {&Y2}, "Sum");
 
  // ---------------------------------------------------------------------------
  // elSpectro
  // ---------------------------------------------------------------------------
  
  //create a V decaying to J/psi pi+pi-
  auto jpsi=particle(443,model(new PhaseSpaceDecay({},{11,-11})));
  
  //mass_distribution(113,new DistTF1{TF1("hhRho","1",0.2,0.9)});
  mass_distribution(113,new DistTF1{TF1("hhRho","TMath::BreitWigner(x,0.775,0.151)",0.2,1.2)});
  auto rho=particle(113,model(new PhaseSpaceDecay({},{211,-211})));

  //mass_distribution(9995,new DistTF1{TF1("hh","1",3.9,4.5)});
  mass_distribution(9995,new DistTF1{TF1("hh","TMath::BreitWigner(x,4.22,0.05)",3.9,4.5)});
  auto v=particle(9995,model(new PhaseSpaceDecay{{jpsi,rho},{}}));
  
  //create eic electroproduction of X + proton
  auto jpac = new JpacModelQ2W{&sum, {v},{2212}, 1, 1, 0.5};
  auto production=eic( ebeamE, pbeamE, jpac );

  //just produce events with st given distribution and no jpac amplitude
  //auto jpac = new DecayModelQ2W{0, new PhaseSpaceDecay( {v},{2212}),new TwoBody_stu{0,1, 0.5 ,0,0}};
  //auto production=eic( ebeamE, pbeamE, jpac );

  elSpectro::LorentzVector elbeam(0,0,-1*ebeamE,elSpectro::escat::E_el(ebeamE));
  elSpectro::LorentzVector prbeam(0,0,pbeamE,elSpectro::escat::E_pr(pbeamE));
  
  auto pionp = static_cast<DecayingParticle*>(rho)->Model()->Products()[0];
  auto pionm = static_cast<DecayingParticle*>(rho)->Model()->Products()[1];
  auto posJ = static_cast<DecayingParticle*>(jpsi)->Model()->Products()[0];
  auto eleJ = static_cast<DecayingParticle*>(jpsi)->Model()->Products()[1];
  
  auto electron = jpac->GetScatteredElectron();
  auto proton = jpac->GetDecayBaryon();
  //Particle* proton =nullptr;
  //for(auto* pp:particles().StableParticles()){
  //  if(pp->Pdg()==2212){
  //    proton=pp;
  //    break;
  //  }
  // }
  // ---------------------------------------------------------------------------
  // Initialize HepMC3
  // ---------------------------------------------------------------------------
  
  writer(new HepMC3Writer{Form("out/vectorJpsiPiPi_hepmc3_elspectro_%d_%d.txt",(int)ebeamE,(int)pbeamE)});
  //writer(new LundWriter{Form("out/vectorJpsiPiPi_lund_elspectro_%d_%d.txt",(int)ebeamE,(int)pbeamE)});

  
  //initilase the generator, may take some time for making distribution tables 
  initGenerator();
  
   
  // ---------------------------------------------------------------------------
  // Generate events
  // ---------------------------------------------------------------------------
  
  gBenchmark->Start("e");

  for(int i=0;i<nEvents;i++){
    if(i%100==0) std::cout<<"event number "<<i<<std::endl;
    nextEvent();

    //fill diagnostic histograms
    auto photon = elbeam - electron->P4();
    double Q2 = -photon.M2();
    double W = (photon+prbeam).M();
    double t = -1*(proton->P4()-prbeam).M2();
    hQ2.Fill(Q2);
    hW.Fill(W);
    ht.Fill(t);
    hgE.Fill(photon.E());
    
    auto elec = electron->P4();
    heTh.Fill(elec.Theta() *TMath::RadToDeg());
    heE.Fill(elec.E());
    
    auto jpsiP4 = eleJ->P4() + posJ->P4();
    hJpsiM.Fill(jpsiP4.M());
    auto rhoP4 = pionp->P4() + pionm->P4();
    hRhoM.Fill(rhoP4.M());
    auto jpsi_rho=rhoP4+jpsiP4;
    hMesonM.Fill(jpsi_rho.M());
    
    hElePVsEta.Fill(eleJ->P4().Eta(), eleJ->P4().P());
    hPosPVsEta.Fill(posJ->P4().Eta(), posJ->P4().P());
    hPionpPVsEta.Fill(pionp->P4().Eta(), pionp->P4().P());
    hPionmPVsEta.Fill(pionm->P4().Eta(), pionm->P4().P());
    hRecoilPVsEta.Fill(proton->P4().Eta(), proton->P4().P());
    hRecoilThetaVsP.Fill(proton->P4().P(), proton->P4().Theta()*1000.);
    hRecoilPt.Fill(proton->P4().Pt());
  
    
  }
  gBenchmark->Stop("e");
  gBenchmark->Print("e");
 
 
  // internally stored histograms for total ep cross section
  TH1D *hWdist = (TH1D*)gDirectory->FindObject("Wdist");
  TH1D *hGenWdist = (TH1D*)gDirectory->FindObject("genWdist");

  TFile *fout = TFile::Open(Form("out/vectorJpsiPiPi_elspectro_%d_%d_diagnostic2.root",(int)ebeamE,(int)pbeamE), "recreate");
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
  hJpsiM.Write();
  hRhoM.Write();
  hMesonM.Write();
  hElePVsEta.Write();
  hPosPVsEta.Write();
  hPionpPVsEta.Write();
  hPionmPVsEta.Write();
  hRecoilPVsEta.Write();
  hRecoilThetaVsP.Write();
  hRecoilPt.Write();
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
     cout<<"Total cross section ("<<(int)ebeamE<<","<<(int)pbeamE<<"): sigma_ep = "<<integrated_xsection<<" nb "<<endl;
   }
  }
}
