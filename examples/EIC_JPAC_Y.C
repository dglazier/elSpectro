//Just need jpacPhoto headers
#include "reaction_kinematics.hpp"
#include "regge_trajectory.hpp"
#include "core/vector_exchange.hpp"
#include "core/pseudoscalar_exchange.hpp"
#include "core/pomeron_exchange.hpp"
#include "core/amplitude_sum.hpp"
#include "core/baryon_resonance.hpp"

#include "FunctionsForGenvector.h"

// ---------------------------------------------------------------------------
// Diagnostic histograms
// ---------------------------------------------------------------------------

TH1F hQ2("Q2","Q2",1000,0,100);
TH1D hgPhi("gPhi","gPhi",100,-180,180);
TH1D heE("eE","eE",1000,0,50);
TH2D heThPhi("eThPhi","eThPhi",50,-180,180,500,178,180);
TH1D heTh("eTh","eTh",1000,0,180);
TH1D hPTh("NTh","NTh",1000,0,180);
TH2D hPThPhi("PThPhi","PThPhi",50,-180,180,50,0,10);
TH1D hPPhi("NPh","NPh",1000,-180,180);
TH1D hePhi("ePh","ePh",90,-180,180);
TH1F hW("W","W",1000,0,100);
TH1F ht("t","t",1000,0,10);
TH1F hgE("gE","gE",1000,0,20);
TH1F hMesonM("MesonM","; J/#psi#pi Mass (GeV)",1000,3.5,5);
TH1F h2PiM("M2Pi",";#pi+#pi- Mass (GeV)",1000,0,2);
TH1F hJpsiM("JpsiM","; e+e- Mass (GeV)",1000,2.9,3.3);
TH2F hElePVsEta("ElePVsEta","; #eta; p (GeV)",200,-5,5,500,0,50);
TH2F hPosPVsEta("PosPVsEta","; #eta; p (GeV)",200,-5,5,500,0,50);
TH2F hPionPVsEta("PionPVsEta","; #eta; p (GeV)",200,-5,5,500,0,50);
TH2F hRecoilPVsEta("RecoilPVsEta","; #eta; p (GeV)",200,0,10,1000,0,275);
TH2F hRecoilThetaVsP("RecoilThetaVsP","; p (GeV); #theta (mrad)",1000,0,275,200,0,200);
TH1F hRecoilPt("RecoilPt","; p_{T} (GeV)",200,0,5.0);


void EIC_JPAC_Y(double ebeamE = 10, double pbeamE = 100, double nLumi=1E33, int nDays = 100) {

  Double_t crossingAngle=0.0; //0mrad
  //define e- beam, pdg =11
  auto elBeam = initial(11,ebeamE);
  auto elBeamP4=elBeam->GetInteracting4Vector();
  elBeam->SetAngleThetaPhi(TMath::Pi()-crossingAngle,0);

  //define pr beam, pdg =2212
  auto prBeam = initial(2212,pbeamE);
  auto prBeamP4=prBeam->GetInteracting4Vector();
  prBeam->SetAngleThetaPhi(crossingAngle,0);
  
  // ---------------------------------------------------------------------------
  // AMPLITUDES
  // ---------------------------------------------------------------------------

  // ---------------------------------------------------------------------------
  // Preliminaries
  // ---------------------------------------------------------------------------



  // Y(4260)
  auto kY=reaction_kinematics(M_Y4260);
  kY.set_meson_JP(1, -1);
  //  double R_Y = 1.55;//changed in JpacPhoto 11/10/2022
  double R_Y = 0.84;

    
  // Same but high-energy
  linear_trajectory alpha_traj_high{1, 1.15, 0.11, "HE"};
  double b_HE = 1.01;
  double A_HE = 0.16;
  pomeron_exchange Y_Amp_high{&kY, &alpha_traj_high, true, "#it{Y}(4260) high"};
  Y_Amp_high.set_params({A_HE * R_Y, b_HE});

  // Low - energy trajectory and couplings
  linear_trajectory alpha_traj_low{1, 0.94, 0.36, "LE"};
  double b_LE = 0.12;
  double A_LE = 0.38;
  pomeron_exchange Y_Amp_low{&kY, &alpha_traj_low, false, "#it{Y}(4260) low"};
  Y_Amp_low.set_params({A_LE * R_Y, b_LE});
    
  double low_s = 7*7.;//GeV
  double high_s = 20*20.;//GeV
  amplitude_blend jpac_amp{&kY,&Y_Amp_low,low_s,&Y_Amp_high,high_s, "#it{X} blend"};

 
  // -------------------------------------------------------------------
  // elSpectro
  // --------------------------------------------------------------------
  //create a V decaying to J/psi pi+pi-
  auto jpsi=particle(443,model(new PhaseSpaceDecay({},{11,-11})));

  //flat 2pi mass distribution
  mass_distribution(9996,new DistTF1{TF1("hhsigma","1",0.25,4)});
  auto sigma=particle(9996,model(new PhaseSpaceDecay({},{211,-211})));

  mass_distribution(9995,new DistTF1{TF1("hh","TMath::BreitWigner(x,4.22,0.05)",3.5,5)});
  auto Y=particle(9995,model(new PhaseSpaceDecay{{jpsi,sigma},{}}));
  Y->SetPdgMass(M_Y4260);
    
  //create eic electroproduction of X + proton
  auto pGammaStarDecay = JpacModelst{&jpac_amp, {Y},{2212} }; //photo-nucleon system
  //Decay g*p state, provide s channel and t-channel "shapes"
  //Note the amplitude will provide the actual t-distribution, this approximation speeds up sampling
  //TwoBody_stu{0., 1.0, 2.5} => 0% s-schannel, 100% t channel with slope 2.5 
  auto photoprod = DecayModelQ2W{0,&pGammaStarDecay,new TwoBody_stu{0., 1.0, 3}};
  //combine beam, target and reaction products
  //auto production=eic( ebeamE, pbeamE, &photoprod );
  auto production=eic( elBeam,prBeam,&photoprod);

  // ---------------------------------------------------------------------------
  // Initialize HepMC3
  // ---------------------------------------------------------------------------
  writer(new HepMC3Writer{Form("outCollision/jpac_Y_%d_%d.txt",(int)ebeamE,(int)pbeamE)});
  // exit(0);
  
  // ---------------------------------------------------------------------------
  //initilase the generator, may take some time for making distribution tables 
  // ---------------------------------------------------------------------------
  initGenerator();
  
  // ---------------------------------------------------------------------------
  //Set number of events via experimental luminosity and beamtime
  // ---------------------------------------------------------------------------
  production->SetCombinedBranchingFraction(0.06*0.01); //Jpsi->e+e- and 1% Y
  generator().SetNEvents_via_LuminosityTime(nLumi,24*60*60*nDays);
  //generator().SetNEvents(10000);
  // auto fastIntegral=production->IntegrateCrossSectionFast();
  //std::cout<<"       check fast cross section "<<fastIntegral<<std::endl;

  
  // exit(0);
  // ---------------------------------------------------------------------------
  // Get event particles for filling histograms
  // ---------------------------------------------------------------------------
  auto pionp = static_cast<DecayingParticle*>(sigma)->Model()->Product(0);
  auto pionm = static_cast<DecayingParticle*>(sigma)->Model()->Product(1);
  auto posJ = static_cast<DecayingParticle*>(jpsi)->Model()->Products()[0];
  auto eleJ = static_cast<DecayingParticle*>(jpsi)->Model()->Products()[1];

  auto electron = photoprod.GetScatteredElectron();
  auto neutron = photoprod.GetDecayBaryon();

  // ---------------------------------------------------------------------------
  // Generate events
  // ---------------------------------------------------------------------------
  
  gBenchmark->Start("generator");//timer

  while(finishedGenerator()==false){
    nextEvent();
    countGenEvent();
    if(generator().GetNDone()%1000==0) std::cout<<"event number "<<generator().GetNDone()<<std::endl;

    auto photon = *elBeamP4 - electron->P4();
   double Q2 = -photon.M2();
   double W = (photon+*prBeamP4).M();
   double t = -1*(neutron->P4()- *prBeamP4).M2();
   hQ2.Fill(Q2);
   hW.Fill(W);
   ht.Fill(t);
   hgE.Fill(photon.E());

   auto forScPhi= electron->P4();
    hgPhi.Fill(forScPhi.Phi()*TMath::RadToDeg());
   
     
   auto elec = electron->P4();
   heThPhi.Fill(forScPhi.Phi() *TMath::RadToDeg(),forScPhi.Theta() *TMath::RadToDeg());
   heTh.Fill(elec.Theta() *TMath::RadToDeg());
   hePhi.Fill(elec.Phi() *TMath::RadToDeg());
   heE.Fill(elec.E());
   
   auto jpsiP4 = eleJ->P4() + posJ->P4();
   hJpsiM.Fill(jpsiP4.M());
   auto sigmaP4 = pionp->P4() + pionm->P4();
   auto jpsi_pion=sigmaP4+jpsiP4;
   hMesonM.Fill(jpsi_pion.M());
   h2PiM.Fill(sigmaP4.M());

   hElePVsEta.Fill(eleJ->P4().Eta(), eleJ->P4().P());
   hPosPVsEta.Fill(posJ->P4().Eta(), posJ->P4().P());
   hPionPVsEta.Fill(pionp->P4().Eta(), pionp->P4().P());
   hRecoilPVsEta.Fill(neutron->P4().Eta(), neutron->P4().P());
   hRecoilThetaVsP.Fill(neutron->P4().P(), neutron->P4().Theta()*1000.);
   hRecoilPt.Fill(ROOT::Math::VectorUtil::Perp(prBeamP4->Vect(),neutron->P4().Vect()));
   
   hPTh.Fill(neutron->P4().Theta()*TMath::RadToDeg());
   hPThPhi.Fill(neutron->P4().Phi()*TMath::RadToDeg(),neutron->P4().Theta()*TMath::RadToDeg());
   hPPhi.Fill(neutron->P4().Phi()*TMath::RadToDeg());


  }
  
  gBenchmark->Stop("generator");
  gBenchmark->Print("generator");
  
  // ---------------------------------------------------------------------------
  // Report generator statistics, can be used for optimising
  // ---------------------------------------------------------------------------
  
  generator().Summary();

   
  // ---------------------------------------------------------------------------
  // Write diagnostic histograms
  // ---------------------------------------------------------------------------
  TH1D *hWdist = (TH1D*)gDirectory->FindObject("Wdist");
  TH1D *hGenWdist = (TH1D*)gDirectory->FindObject("genWdist");

  TFile *fout = TFile::Open(Form("out/jpac_y_%d_%d.txt",(int)ebeamE,(int)pbeamE), "recreate");
  // generated event distributions
  hQ2.Write();
  hW.Write();
  ht.Write();
  hgE.Write();
  hgPhi.Write();
  heTh.Write();
  hPTh.Write();
  hPPhi.Write();
  hePhi.Write();
  heE.Write();
  hJpsiM.Write();
  hMesonM.Write();
  h2PiM.Write();
  hElePVsEta.Write();
  hPosPVsEta.Write();
  hPionPVsEta.Write();
  hRecoilPVsEta.Write();
  hRecoilThetaVsP.Write();
  hRecoilPt.Write();
  hPThPhi.Write();
  heThPhi.Write();
  hWdist->Write();
  fout->Close();
  
}
