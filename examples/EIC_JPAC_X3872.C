//Just need jpacPhoto headers
#include "reaction_kinematics.hpp"
#include "regge_trajectory.hpp"
#include "amplitudes/vector_exchange.hpp"
#include "amplitudes/pseudoscalar_exchange.hpp"
#include "amplitudes/pomeron_exchange.hpp"
#include "amplitudes/amplitude_sum.hpp"
#include "amplitudes/baryon_resonance.hpp"


void EIC_JPAC_X3872(double ebeamE = 5, double pbeamE = 41, int nEvents = 5e4) {

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
  //x
  mass_distribution(9995,new DistTF1{TF1("hh","TMath::BreitWigner(x,3.872,0.001)",3.85,3.89)});
  auto x=particle(9995,model(new PhaseSpaceDecay{{jpsi,rho},{}}));
  x->SetPdgMass(3.872);

  //create eic electroproduction of X + proton
  auto pGammaStarDecay = new JpacModelst{&jpac_amp, {x},{2212} }; //photo-nucleon system
  auto photoprod = new DecayModelQ2W{0,pGammaStarDecay,new TwoBody_stu{0., 1.0, 2.5,0,0}};

  //combine beam, target and reaction products
  auto production=eic( ebeamE, pbeamE, photoprod );

  // ---------------------------------------------------------------------------
  // Initialize HepMC3
  // ---------------------------------------------------------------------------
  writer(new HepMC3Writer{Form("out/x3872_pac_%d_%d.txt",(int)ebeamE,(int)pbeamE)});
  
  
  //initilase the generator, may take some time for making distribution tables 
  initGenerator();
  
  // ---------------------------------------------------------------------------
  // Generate events
  // ---------------------------------------------------------------------------
  
  gBenchmark->Start("e");//timer
  
  for(int i=0;i<nEvents;i++){
     nextEvent();
     if(i%1000==0) std::cout<<"event number "<<i<<std::endl;
  }

  gBenchmark->Stop("e");
  gBenchmark->Print("e");

  // ---------------------------------------------------------------------------
  // Report generator statistics, can be used for optimising
  // ---------------------------------------------------------------------------
  
  generator().Summary();

}

