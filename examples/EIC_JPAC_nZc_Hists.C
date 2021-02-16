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

TH1F hQ2("Q2","Q2",1000,0,100);
TH1D heE("eE","eE",1000,0,20);
TH1D heTh("eTh","eTh",1000,0,180);
TH1D hePhi("ePh","ePh",90,-180,180);
TH1F hW("W","W",1000,0,100);
TH1F ht("t","t",1000,0,10);
TH1F hgE("gE","gE",1000,0,20);
TH1F hMesonM("MesonM","; J/#psi#pi Mass (GeV)",1000,3.5,4.3);
TH1F hJpsiM("JpsiM","; e+e- Mass (GeV)",1000,2.9,3.3);
TH2F hElePVsEta("ElePVsEta","; #eta; p (GeV)",200,-5,5,500,0,50);
TH2F hPosPVsEta("PosPVsEta","; #eta; p (GeV)",200,-5,5,500,0,50);
TH2F hPionPVsEta("PionPVsEta","; #eta; p (GeV)",200,-5,5,500,0,50);
TH2F hRecoilPVsEta("RecoilPVsEta","; #eta; p (GeV)",200,0,10,1000,0,275);
TH2F hRecoilThetaVsP("RecoilThetaVsP","; p (GeV); #theta (mrad)",1000,0,275,200,0,200);
TH1F hRecoilPt("RecoilPt","; p_{T} (GeV)",200,0,5.0);


//Amplitude based on $JPACPHOTO/executables/XYZ_Plots/Z_low.cpp and Z_high.cpp

//To set luminosity and days change last 2 arguments
//e.g. for luminosoty 10^33 and 25 days, e- energy 100 and p energy 100
//with high energy paramterisation :
// 'EIC_JPAC_X3872.C("high",100,100,1E33,25)'
// To just run a fixed number of events leave last
// argument 0 and nLumi=number of events
// 'EIC_JPAC_X3872.C("high",100,100,1E4)'

void EIC_JPAC_nZc_Hists(string ampPar="high",double ebeamE = 5, double pbeamE = 41, double nLumi=100, int nDays = 0) {

  LorentzVector elbeam(0,0,-1*ebeamE,escat::E_el(ebeamE));
  LorentzVector prbeam(0,0,pbeamE,escat::E_pr(pbeamE));

  // ---------------------------------------------------------------------------
  // AMPLITUDES
  // ---------------------------------------------------------------------------
 // ---------------------------------------------------------------------------
  // Preliminaries
  // ---------------------------------------------------------------------------

  double g_NN = sqrt(4. * M_PI * 13.81); // Nucleon coupling same for all
  double LamPi = .9;  // 900 MeV cutoff for formfactor
  double bPi = 1. / (LamPi * LamPi);

  // Zc(3900)
  auto kZc = reaction_kinematics{M_ZC3900};
  kZc.set_JP(1, 1);

  double gc_Psi = 1.91; // psi coupling before VMD scaling
  double gc_Gamma = E * F_JPSI * gc_Psi / M_JPSI;
  std::vector<double> Zc_couplings = {gc_Gamma, g_NN};

  // Pion trajectory 
  int signature = +1;
  double alpha_prime = 0.7; // GeV^-2
  double alpha_0 =  - alpha_prime * M2_PION;
  auto alpha = linear_trajectory{signature, alpha_0, alpha_prime};

  pseudoscalar_exchange* ampZc{nullptr};
  // ---------------------------------------------------------------------------
  // low => Fixed-spin amplitudes
  // ---------------------------------------------------------------------------
  if(ampPar=="low") ampZc = new pseudoscalar_exchange{&kZc, M_PION, "#it{Z_{c}} (3900)^{+}"};
 // ---------------------------------------------------------------------------
  // high => Reggeized amplitudes
  // ---------------------------------------------------------------------------
  else if(ampPar=="high") ampZc = new pseudoscalar_exchange(&kZc, &alpha, "#it{Z_{c}}(3900)^{+}");
  else {cerr<<"invalid amplitude parameterisation "<<ampPar<<endl; exit(0);}

  ampZc->set_params(Zc_couplings);
  ampZc->set_formfactor(true, bPi);


  // ---------------------------------------------------------------------------
  // elSpectro
  // ---------------------------------------------------------------------------
  //create a Z decaying to J/psi pi+

  auto jpsi=particle(443,model(new PhaseSpaceDecay({},{11,-11})));

  //Zc
  double minMass = 3.5;
  double maxMass = 4.4;
  mass_distribution(9995,new DistTF1{TF1("hh",Form("TMath::BreitWigner(x,%lf,0.05)",M_ZC3900),minMass,maxMass)});
  auto Z=particle(9995,model(new PhaseSpaceDecay{{jpsi},{211}}));
  Z->SetPdgMass(M_ZC3900);

  //create eic electroproduction of X + proton
  auto pGammaStarDecay = JpacModelst{ampZc, {Z},{2112} }; //photo-nucleon system
  //auto pGammaStarDecay = DecayModelst{ {Z},{2212} }; //photo-nucleon system
  //Decay g*p state, provide s channel and t-channel "shapes"
  //Note the amplitude will provide the actual t-distribution, this approximation speeds up sampling
  //TwoBody_stu{0., 1.0, 2.5} => 0% s-schannel, 100% t channel with slope 2.5 
  auto photoprod = DecayModelQ2W{0,&pGammaStarDecay,new TwoBody_stu{0., 1.0, 2.5}};

  //combine beam, target and reaction products
  auto production=eic( ebeamE, pbeamE, &photoprod );

  // ---------------------------------------------------------------------------
  // Initialize HepMC3
  // ---------------------------------------------------------------------------
  writer(new HepMC3Writer{Form("out/jpac_Zc3900_%s_%d_%d.txt",ampPar.data(),(int)ebeamE,(int)pbeamE)});
  
  
  // ---------------------------------------------------------------------------
  //initilase the generator, may take some time for making distribution tables 
  // ---------------------------------------------------------------------------
  initGenerator();
  
  // ---------------------------------------------------------------------------
  //Set number of events via experimental luminosity and beamtime
  // ---------------------------------------------------------------------------
  production->SetCombinedBranchingFraction(0.06); //Just Jpsi->e+e-
  generator().SetNEvents_via_LuminosityTime(nLumi,24*60*60*nDays);

  auto fastIntegral=production->IntegrateCrossSectionFast();
  std::cout<<"       check fast cross section "<<fastIntegral<<std::endl;

  // ---------------------------------------------------------------------------
  // Get event particles for filling histograms
  // ---------------------------------------------------------------------------
  auto pionp = static_cast<DecayingParticle*>(Z)->Model()->Products()[1];
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

      auto photon = elbeam - electron->P4();
   double Q2 = -photon.M2();
   double W = (photon+prbeam).M();
   double t = -1*(neutron->P4()-prbeam).M2();
   hQ2.Fill(Q2);
   hW.Fill(W);
   ht.Fill(t);
   hgE.Fill(photon.E());
   
   //LorentzVector totalIn = elbeam + prbeam;
   //LorentzVector totalOut = electron->P4() + neutron->P4() + eleJ->P4() + posJ->P4() + pionp->P4();
   //cout<<totalIn.E()<<" "<<totalOut.E()<<endl;

   auto elec = electron->P4();
   heTh.Fill(elec.Theta() *TMath::RadToDeg());
   hePhi.Fill(elec.Phi() *TMath::RadToDeg());
   heE.Fill(elec.E());
   
   auto jpsiP4 = eleJ->P4() + posJ->P4();
   hJpsiM.Fill(jpsiP4.M());
   auto jpsi_pion=pionp->P4()+jpsiP4;
   hMesonM.Fill(jpsi_pion.M());
   
   hElePVsEta.Fill(eleJ->P4().Eta(), eleJ->P4().P());
   hPosPVsEta.Fill(posJ->P4().Eta(), posJ->P4().P());
   hPionPVsEta.Fill(pionp->P4().Eta(), pionp->P4().P());
   hRecoilPVsEta.Fill(neutron->P4().Eta(), neutron->P4().P());
   hRecoilThetaVsP.Fill(neutron->P4().P(), neutron->P4().Theta()*1000.);
   hRecoilPt.Fill(neutron->P4().Pt());


  }
  
  gBenchmark->Stop("generator");
  gBenchmark->Print("generator");
  
  // ---------------------------------------------------------------------------
  // Report generator statistics, can be used for optimising
  // ---------------------------------------------------------------------------
  
  generator().Summary();

  delete ampZc;
  
  // ---------------------------------------------------------------------------
  // Write diagnostic histograms
  // ---------------------------------------------------------------------------

  TFile *fout = TFile::Open(Form("out/jpac_Zc3900_%s_%d_%d_diagnostic.root",ampPar.data(),(int)ebeamE,(int)pbeamE), "recreate");
  // generated event distributions
  hQ2.Write();
  hW.Write();
  ht.Write();
  hgE.Write();
  heTh.Write();
  heE.Write();
  hJpsiM.Write();
  hMesonM.Write();
  hElePVsEta.Write();
  hPosPVsEta.Write();
  hPionPVsEta.Write();
  hRecoilPVsEta.Write();
  hRecoilThetaVsP.Write();
  hRecoilPt.Write();
  fout->Close();

}
