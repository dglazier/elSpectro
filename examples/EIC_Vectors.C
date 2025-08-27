#include "BuildReaction.h"
#include "ParticleFactory.h"
#include "ElectronScattering.h"
#include "JpacAmplitude.h"
#include "Manager.h"
#include "HepMC3Writer.h"
#include "Vector_mesons.h"
#include "Jpsi.h"


//void EIC_Vectors(double ebeamE = 26.7, double pbeamE = 820, double nLumi=2E33, double nDays = 1./24/60/10*2000) {
void EIC_Vectors(double ebeamE = 10, double pbeamE = 100, double nLumi=6.1E33, double nDays = 0.1/24/60) {
  gRandom->SetSeed(1111);
  //create a generator 
  auto generator = elSpectro::Manager{};
  elSpectro::particles::ParticleFactory::Init();
  auto& particleFactory = elSpectro::particles::ParticleFactory::Instance();
  // ---------------------------------------------------------------------------
  // TWO-BODY PARTICLES
  // ---------------------------------------------------------------------------
   //int proton_id = 2212; //baryon stable
  //auto proton = elSpectro::particles::ParticleFactory::Instance().Get("proton"); //meson decays
 
  // ---------------------------------------------------------------------------
  // AMPLITUDES
  // ---------------------------------------------------------------------------
  auto omega = particleFactory.CreateDecayingParticle("omega"); //meson decays
  jpacAmps::omega omega_amp;
  auto jpac_omega = omega_amp.Get();
  auto two_body_omega = elSpectro::JpacTwoBody{jpac_omega,{ omega },{2212} };

  auto rho=particleFactory.CreateDecayingParticle("rho0");
  jpacAmps::rho rho_amp;
  auto jpac_rho = rho_amp.Get();
  auto two_body_rho = elSpectro::JpacTwoBody{jpac_rho,{ rho },{2212} }; 
 
  auto phi=particleFactory.CreateDecayingParticle("phi");
  jpacAmps::phi phi_amp;
  auto jpac_phi = phi_amp.Get();
  auto two_body_phi = elSpectro::JpacTwoBody{jpac_phi,{ phi },{2212} }; 
 
  auto Jpsi = particleFactory.CreateDecayingParticle("J/psi"); //meson decays
  auto two_body_jpsi = elSpectro::JpacTwoBody{jpacAmps::Jpsi_P().Get(),{ Jpsi },{2212} };

  // ---------------------------------------------------------------------------
  // PRODUCTION PROCESS
  // ---------------------------------------------------------------------------
 
  //photo-nucleon system decaying to meson and baryon

  
  auto elBeam =  elSpectro::CollidingParticle{11,ebeamE};
  auto elBeamP4=elBeam.GetInteracting4Vector();
  elBeam.SetAngleThetaPhi(TMath::Pi(),0);

  //define pr beam, pdg =2212
  auto prBeam =  elSpectro::CollidingParticle{2212,pbeamE};
  auto prBeamP4=prBeam.GetInteracting4Vector();


  ////////*****************************
  ////Need to have FormationQ2W without taking model and decayer
  ////Then add decays
  //  auto formation =  elSpectro::FormationQ2W{1.,elSpectro::CloneModel(two_body_jpsi)};
  auto formation =  elSpectro::FormationQ2W{1.,elSpectro::CloneModel(two_body_rho)};
  formation.AddHadronicChannel(1,elSpectro::CloneModel(two_body_omega));
  formation.AddHadronicChannel(1,elSpectro::CloneModel(two_body_phi));
  formation.AddHadronicChannel(1,elSpectro::CloneModel(two_body_jpsi));

  formation.setThreshold(2.0);
  
  auto allProduction = new elSpectro::ElectronScattering(elBeam,prBeam, elSpectro::CloneModel( formation ) );

  // allProduction->SetLimit_Q2min(2);

 
  //create eic electroproduction of X + proton
  //allProduction now owned by generator
  generator.Reaction( allProduction);

  // generator.Reaction( phoPro::Build_ep_Collision(ebeamE,pbeamE,two_body_omega));
 
  // ---------------------------------------------------------------------------
  // Initialize HepMC3
  // ---------------------------------------------------------------------------
  generator.SetWriter(new HepMC3Writer{Form("/home/dglazier/dump/jpac_rho_%d_%d_testing.txt",(int)ebeamE,(int)pbeamE)});
  
  // ---------------------------------------------------------------------------
  //initilase the generator, may take some time for making distribution tables 
  // ---------------------------------------------------------------------------
  generator.InitGeneration();
  
  // ---------------------------------------------------------------------------
  //Set number of events via experimental luminosity and beamtime
  // ---------------------------------------------------------------------------
  //omega_prod->SetCombinedBranchingFraction(omega.BranchRatio()); //Just Jpsi->e+e-
  generator.SetNEvents_via_LuminosityTime(nLumi,24*60*60*nDays);
  generator.SetNEvents(10000);
  //auto fastIntegral=rho_prod->IntegrateCrossSectionFast();
  // std::cout<<"       check fast cross section "<<fastIntegral<<std::endl;

  
  // ---------------------------------------------------------------------------
  //Generate events
  // ---------------------------------------------------------------------------
  std::cout<<" ***************** start processing "<<endl;
  generator.processAll();

 
 }
