//Just need jpacPhoto headers
#include "reaction_kinematics.hpp"
#include "regge_trajectory.hpp"
#include "amplitudes/vector_exchange.hpp"
#include "amplitudes/pseudoscalar_exchange.hpp"
#include "amplitudes/pomeron_exchange.hpp"
#include "amplitudes/amplitude_sum.hpp"
#include "amplitudes/baryon_resonance.hpp"

// ---------------------------------------------------------------------------
// Diagnostic histograms
// ---------------------------------------------------------------------------

double minMass = 3.6;
double maxMass = 4.0;
TH1F hQ2("Q2","Q2",1000,0,100);
TH1D heE("eE","eE",1000,0,20);
TH1D heTh("eTh","eTh",1000,0,180);
TH1F hW("W","W",1000,0,100);
TH1F ht("t","t",1000,0,10);
TH1F hgE("gE","gE",1000,0,20);
TH1F hMesonM("MesonM","; J/#psi#pi Mass (GeV)",1000,minMass,maxMass);
TH1F hJpsiM("JpsiM","; e+e- Mass (GeV)",1000,2.9,3.3);
TH1F hRhoM("RhoM","; #pi+#pi- Mass (GeV)",1000,0.2,1.2);
TH2F hElePVsEta("ElePVsEta","; #eta; p (GeV)",200,-5,5,500,0,50);
TH2F hPosPVsEta("PosPVsEta","; #eta; p (GeV)",200,-5,5,500,0,50);
TH2F hPionpPVsEta("PionpPVsEta","; #eta; p (GeV)",200,-5,5,500,0,50);
TH2F hPionmPVsEta("PionmPVsEta","; #eta; p (GeV)",200,-5,5,500,0,50);
TH2F hRecoilPVsEta("RecoilPVsEta","; #eta; p (GeV)",200,0,10,1000,0,275);
TH2F hRecoilThetaVsP("RecoilThetaVsP","; p (GeV); #theta (mrad)",1000,0,275,200,0,200);
TH1F hRecoilPt("RecoilPt","; p_{T} (GeV)",200,0,5.0);


void EIC_JPAC_X3872_Hists(double ebeamE = 5, double pbeamE = 41, int nEvents = 5e4) {

  LorentzVector elbeam(0,0,-1*ebeamE,escat::E_el(ebeamE));
  LorentzVector prbeam(0,0,pbeamE,escat::E_pr(pbeamE));

  // ---------------------------------------------------------------------------
  // AMPLITUDES
  // ---------------------------------------------------------------------------

  // Linear trajectory for the rho
  linear_trajectory alpha(-1, 0.5, 0.9, "EXD_linear");
  
  // Set up kinematics for the X(3872)
  reaction_kinematics * ptr = new reaction_kinematics(3.872, "X(3872)");

  // Initialize Reggeon amplitude with the above kinematics and regge_trajectory
  vector_exchange X3872rho(ptr, &alpha, "#rho");
  X3872rho.set_params({3.81E-3, 2.4, 14.6});
  
  vector_exchange X3872omega(ptr, &alpha, "#omega");
  X3872omega.set_params({9.51E-3, 16, 0.});
  
  // Sum individual contributions together incoherently
  amplitude_sum jpac_amp(ptr, {&X3872rho, &X3872omega}, "Sum");

  // ---------------------------------------------------------------------------
  // elSpectro
  // ---------------------------------------------------------------------------
  
  //create a X decaying to J/psi pi+pi-
  auto jpsi=particle(443,model(new PhaseSpaceDecay({},{11,-11})));
  //rho
  mass_distribution(113,new DistTF1{TF1("hhRho","TMath::BreitWigner(x,0.775,0.151)",0.2,0.7)});
  auto rho=particle(113,model(new PhaseSpaceDecay({},{211,-211})));
  //x : new resonance needs dummy PDG code = 9995
  //    and mass = 3.872 as well as mass distribution
  mass_distribution(9995,new DistTF1{TF1("hh","TMath::BreitWigner(x,3.872,0.001)",3.85,3.89)});
  auto x=particle(9995,3.872,model(new PhaseSpaceDecay{{jpsi,rho},{}}));

  //create eic electroproduction of X + proton
  auto pGammaStarDecay = new JpacModelst{&jpac_amp, {x},{2212} }; //photo-nucleon system
  auto photoprod = new DecayModelQ2W{0,pGammaStarDecay,new TwoBody_stu{0., 1.0, 2.5,0,0}};

  //combine beam, target and reaction products
  auto production=eic( ebeamE, pbeamE, photoprod );

  // ---------------------------------------------------------------------------
  // Initialize HepMC3
  // ---------------------------------------------------------------------------
  writer(new HepMC3Writer{Form("out/x3872_pac_%d_%d.txt",(int)ebeamE,(int)pbeamE)});

  
 // ---------------------------------------------------------------------------
  // Get event particles for filling histograms
  // ---------------------------------------------------------------------------
  auto pionp = static_cast<DecayingParticle*>(rho)->Model()->Products()[0];
  auto pionm = static_cast<DecayingParticle*>(rho)->Model()->Products()[1];
  auto posJ = static_cast<DecayingParticle*>(jpsi)->Model()->Products()[0];
  auto eleJ = static_cast<DecayingParticle*>(jpsi)->Model()->Products()[1];
  
  auto electron = photoprod->GetScatteredElectron();
  auto proton = photoprod->GetDecayBaryon();
 
  //initilase the generator, may take some time for making distribution tables 
  initGenerator();
  
  // ---------------------------------------------------------------------------
  // Generate events
  // ---------------------------------------------------------------------------
  
  gBenchmark->Start("e");//timer
  
  for(int i=0;i<nEvents;i++){
     nextEvent();
     if(i%1000==0) std::cout<<"event number "<<i<<std::endl;

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

  // ---------------------------------------------------------------------------
  // Report generator statistics, can be used for optimising
  // ---------------------------------------------------------------------------
  
  generator().Summary();

  
  // ---------------------------------------------------------------------------
  // Write diagnostic histograms
  // ---------------------------------------------------------------------------
  TFile *fout = TFile::Open(Form("out/x3872_elspectro_%d_%d_diagnostic.root",(int)ebeamE,(int)pbeamE), "recreate");
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

}

