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

//Amplitude based on $JPACPHOTO/executables/XYZ_Plots/X3872_high.cpp

//To set luminosity and days change last 2 arguments
//e.g. for luminosoty 10^33 and 25 days, e- energy 100 and p energy 100
//with high energy paramterisation :
// 'EIC_JPAC_X3872.C("high",100,100,1E33,25)'
// To just run a fixed number of events leave last
// argument 0 and nLumi=number of events
// 'EIC_JPAC_X3872.C("high",100,100,1E4)'

void EIC_JPAC_X3872_Hists(string ampPar="high",double ebeamE = 5, double pbeamE = 41, double nLumi=100, int nDays = 0) {

  LorentzVector elbeam(0,0,-1*ebeamE,escat::E_el(ebeamE));
  LorentzVector prbeam(0,0,pbeamE,escat::E_pr(pbeamE));

  // ---------------------------------------------------------------------------
  // AMPLITUDES
  // ---------------------------------------------------------------------------

  // ---------------------------------------------------------------------------
  // Preliminaries
  // ---------------------------------------------------------------------------


  // X(3872)
  auto kX = reaction_kinematics{M_X3872};
  kX.set_JP(1, 1);

  // Nucleon couplings and cutoffs
  double gV_omega = 16., gT_omega = 0.;
  double LamOmega = 1.2;
  double gV_rho = 2.4, gT_rho = 14.6;
  double LamRho = 1.4;
  double gV_phi = -6.2, gT_phi = 2.1;
  double gV_psi = 1.6E-3, gT_psi = 0.;

  // Top couplings
  double gX_omega = 8.2E-3;
  double gX_rho = 3.6E-3;
    
  // Linear trajectory for the rho
  auto alpha = linear_trajectory{-1, 0.5, 0.9, "#rho - #omega"};

  // ---------------------------------------------------------------------------
  // High-Energy Amplitudes
  // ---------------------------------------------------------------------------
  //////////////////
  // X(3872)
  vector_exchange *X_omega{nullptr};
  if(ampPar=="high")X_omega= new vector_exchange(&kX, &alpha, "#omega");
  else if(ampPar=="low") X_omega=new vector_exchange(&kX, M_OMEGA, "#omega");
  else {cerr<<"invalid amplitude parameterisation "<<ampPar<<endl; exit(0);}
  X_omega->set_params({gX_omega, gV_omega, gT_omega});
  X_omega->set_formfactor(true, LamOmega);

  vector_exchange *X_rho{nullptr};
  if(ampPar=="high")X_rho= new vector_exchange(&kX, &alpha, "#rho");
  else if(ampPar=="low") X_rho=new vector_exchange(&kX, M_RHO, "#rho");
  else {cerr<<"invalid amplitude parameterisation "<<ampPar<<endl; exit(0);}
  X_rho->set_params({gX_rho, gV_rho, gT_rho});
  X_rho->set_formfactor(true, LamRho);

  std::vector<amplitude*> X_exchanges = {X_omega, X_rho};
  amplitude_sum jpac_amp(&kX, X_exchanges, "#it{X}(3872)");
 
  
  // ---------------------------------------------------------------------------
  // elSpectro
  // ---------------------------------------------------------------------------
  
  //create a X decaying to J/psi pi+pi-
  auto jpsi=particle(443,model(new PhaseSpaceDecay({},{11,-11})));
  //rho
  mass_distribution(113,new DistTF1{TF1("hhRho","TMath::BreitWigner(x,0.775,0.151)",0.2,0.7)});
  auto rho=particle(113,model(new PhaseSpaceDecay({},{211,-211})));
  //x
  mass_distribution(9995,new DistTF1{TF1("hh","TMath::BreitWigner(x,3.872,0.001)",3.85,3.89)});
  auto x=particle(9995,model(new PhaseSpaceDecay{{jpsi,rho},{}}));
  x->SetPdgMass(3.872);

  //create eic electroproduction of X + proton
  auto pGammaStarDecay = JpacModelst{&jpac_amp, {x},{2212} }; //photo-nucleon system
  //Decay g*p state, provide s channel and t-channel "shapes"
  //Note the amplitude will provide the actual t-distribution, this approximation speeds up sampling
  //TwoBody_stu{0., 1.0, 2.5} => 0% s-schannel, 100% t channel with slope 2.5 
  auto photoprod = DecayModelQ2W{0,&pGammaStarDecay,new TwoBody_stu{0., 1.0, 2.5}};

  //combine beam, target and reaction products
  auto production=eic( ebeamE, pbeamE, &photoprod );

  // ---------------------------------------------------------------------------
  // Initialize HepMC3
  // ---------------------------------------------------------------------------
  writer(new HepMC3Writer{Form("out/jpac_x3872_%s_%d_%d.txt",ampPar.data(),(int)ebeamE,(int)pbeamE)});
  
  
  // ---------------------------------------------------------------------------
  //initilase the generator, may take some time for making distribution tables 
  // ---------------------------------------------------------------------------
  initGenerator();
  
  // ---------------------------------------------------------------------------
  //Set number of events via experimental luminosity and beamtime
  // ---------------------------------------------------------------------------
  production->SetCombinedBranchingFraction(0.06); //Just Jpsi->e+e-
  generator().SetNEvents_via_LuminosityTime(nLumi,24*60*60*nDays);
  
  // ---------------------------------------------------------------------------
  // Get event particles for filling histograms
  // ---------------------------------------------------------------------------
  auto pionp = static_cast<DecayingParticle*>(rho)->Model()->Products()[0];
  auto pionm = static_cast<DecayingParticle*>(rho)->Model()->Products()[1];
  auto posJ = static_cast<DecayingParticle*>(jpsi)->Model()->Products()[0];
  auto eleJ = static_cast<DecayingParticle*>(jpsi)->Model()->Products()[1];
  
  auto electron = photoprod.GetScatteredElectron();
  auto proton = photoprod.GetDecayBaryon();
  // ---------------------------------------------------------------------------
  // Generate events
  // ---------------------------------------------------------------------------
  
  gBenchmark->Start("generator");//timer
  
  // for(int i=0;i<nEvents;i++){
  while(finishedGenerator()==false){
    nextEvent();
    countGenEvent();
    if(generator().GetNDone()%1000==0) std::cout<<"event number "<<generator().GetNDone()<<std::endl;

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
  
  gBenchmark->Stop("generator");
  gBenchmark->Print("generator");
  
  // ---------------------------------------------------------------------------
  // Report generator statistics, can be used for optimising
  // ---------------------------------------------------------------------------
  
  generator().Summary();

  delete X_rho;
  delete X_omega;

  // ---------------------------------------------------------------------------
  // Write diagnostic histograms
  // ---------------------------------------------------------------------------
  TFile *fout = TFile::Open(Form("out/jpac_x3872_%s_%d_%d_diagnostic.root",ampPar.data(),(int)ebeamE,(int)pbeamE), "recreate");
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

