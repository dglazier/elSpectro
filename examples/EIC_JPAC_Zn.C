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

void EIC_JPAC_Zn(string ampPar="high",double ebeamE = 5, double pbeamE = 41, double nLumi=100, int nDays = 0) {

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
  double mZc = 3.8884; // GeV
  reaction_kinematics * kZc = new reaction_kinematics(mZc);
  kZc->set_JP(1, 1);

  double gc_Psi = 1.91; // psi coupling before VMD scaling
  double gc_Gamma = E * F_JPSI * gc_Psi / M_JPSI;
  std::vector<double> Zc_couplings = {gc_Gamma, g_NN};

  // Zb(10610)
  double mZb = 10.6072;
  reaction_kinematics * kZb = new reaction_kinematics(mZb);
  kZb->set_JP(1, 1);

  double gb_Ups1 = 0.49, gb_Ups2 = 3.30, gb_Ups3 = 9.22;
  double gb_Gamma = E * (F_UPSILON1S * gb_Ups1 / M_UPSILON1S 
                       + F_UPSILON2S * gb_Ups2 / M_UPSILON2S
                       + F_UPSILON3S * gb_Ups3 / M_UPSILON3S);  
  std::vector<double> Zb_couplings = {gb_Gamma, g_NN};

  
  // Zb(10650)
  double mZbp = 10.6522;
  reaction_kinematics * kZbp = new reaction_kinematics(mZbp);
  kZbp->set_JP(1, 1);

  double gbp_Ups1 = 0.21, gbp_Ups2 = 1.47, gbp_Ups3 = 4.8;
  double gbp_Gamma = E * (F_UPSILON1S * gbp_Ups1 / M_UPSILON1S 
                       +  F_UPSILON2S * gbp_Ups2 / M_UPSILON2S
                       +  F_UPSILON3S * gbp_Ups3 / M_UPSILON3S);  
  std::vector<double> Zbp_couplings = {gbp_Gamma, g_NN};
  
  // ---------------------------------------------------------------------------
  // Fixed-spin amplitudes
  // ---------------------------------------------------------------------------

  pseudoscalar_exchange Zc_fixedspin(kZc, M_PION, "#it{Z_{c}} (3900)^{+}");
  Zc_fixedspin.set_params(Zc_couplings);
  Zc_fixedspin.set_formfactor(true, bPi);

  pseudoscalar_exchange Zb_fixedspin(kZb, M_PION,  "#it{Z_{b}} (10610)^{+}");
  Zb_fixedspin.set_params(Zb_couplings);
  Zb_fixedspin.set_formfactor(true, bPi);

  pseudoscalar_exchange Zbp_fixedspin(kZbp, M_PION, "#it{Z'_{b}} (10650)^{+}");
  Zbp_fixedspin.set_params(Zbp_couplings);
  Zbp_fixedspin.set_formfactor(true, bPi);


  
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

  delete X_rho;
  delete X_omega;

}
