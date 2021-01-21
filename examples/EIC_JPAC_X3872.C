//Just need jpacPhoto headers
#include "reaction_kinematics.hpp"
#include "regge_trajectory.hpp"
#include "amplitudes/vector_exchange.hpp"
#include "amplitudes/pseudoscalar_exchange.hpp"
#include "amplitudes/pomeron_exchange.hpp"
#include "amplitudes/amplitude_sum.hpp"
#include "amplitudes/baryon_resonance.hpp"


//Amplitude based on $JPACPHOTO/executables/XYZ_Plots/X3872_high.cpp
double Frixione(Double_t *x,Double_t *p);

void EIC_JPAC_X3872(string ampPar="high",double ebeamE = 5, double pbeamE = 41, double lumi=1E33, int nDays = 25) {

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
  generator().SetNEvents_via_LuminosityTime(lumi,24*60*60*nDays);
  //or can just do generator().SetNEvents(1E6);
  auto fastIntegral=production->IntegrateCrossSectionFast();
  std::cout<<" check fast cross section "<<fastIntegral<<std::endl;
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
