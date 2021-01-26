//Just need jpacPhoto headers
#include "reaction_kinematics.hpp"
#include "regge_trajectory.hpp"
#include "amplitudes/vector_exchange.hpp"
#include "amplitudes/pseudoscalar_exchange.hpp"
#include "amplitudes/pomeron_exchange.hpp"
#include "amplitudes/amplitude_sum.hpp"
#include "amplitudes/baryon_resonance.hpp"


//Amplitude based on $JPACPHOTO/executables/XYZ_Plots/Z_low.cpp and Z_high.cpp

//To set luminosity and days change last 2 arguments
//e.g. for luminosoty 10^33 and 25 days, e- energy 100 and p energy 100
//with high energy paramterisation :
// 'EIC_JPAC_X3872.C("high",100,100,1E33,25)'
// To just run a fixed number of events leave last
// argument 0 and nLumi=number of events
// 'EIC_JPAC_X3872.C("high",100,100,1E4)'

void EIC_JPAC_nZc(string ampPar="high",double ebeamE = 5, double pbeamE = 41, double nLumi=100, int nDays = 0) {

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

  //or can just do generator().SetNEvents(1E6);
  auto fastIntegral=production->IntegrateCrossSectionFast();
  std::cout<<"       check fast cross section "<<fastIntegral<<std::endl;
  // ---------------------------------------------------------------------------
  // Generate events
  // ---------------------------------------------------------------------------
  
  gBenchmark->Start("generator");//timer

  while(finishedGenerator()==false){
    nextEvent();
    countGenEvent();
    if(generator().GetNDone()%1000==0) std::cout<<"event number "<<generator().GetNDone()<<std::endl;
  }
  
  gBenchmark->Stop("generator");
  gBenchmark->Print("generator");
  
  // ---------------------------------------------------------------------------
  // Report generator statistics, can be used for optimising
  // ---------------------------------------------------------------------------
  
  generator().Summary();

  delete ampZc;

}
