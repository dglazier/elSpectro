#include "ParticleFactory.h"
#include "ElectronScattering.h"
#include "JpacAmplitude.h"
#include "Manager.h"
#include "HepMC3Writer.h"
#include "Vector_mesons.h"

#include "BuildReaction.h"

//void EIC_rho(double ebeamE = 10, double pbeamE = 0, double nLumi=6.1E33, double nDays = 1./24/60/60) {
void EIC_rho(double ebeamE = 18, double pbeamE = 275, double nLumi=6.1E33, double nDays = 1./24/60/60/10) {

  //create a generator 
  auto generator = elSpectro::Manager{};
  elSpectro::particles::ParticleFactory::Init();
  auto& particleFactory = elSpectro::particles::ParticleFactory::Instance();
  // ---------------------------------------------------------------------------
  // TWO-BODY PARTICLES
  // ---------------------------------------------------------------------------
  //phoPro::RhoMeson rho; //meson decays
  auto rho=particleFactory.CreateDecayingParticle("rho0");
  
  int proton_id = 2212; //baryon stable
  
  // ---------------------------------------------------------------------------
  // AMPLITUDES
  // ---------------------------------------------------------------------------
  jpacAmps::rho rho_amp;
  auto jpac_amp = rho_amp.Get();
  // auto jpac_amp = rho_amp.Amp_Pomeron();

  
  // ---------------------------------------------------------------------------
  // PRODUCTION PROCESS
  // ---------------------------------------------------------------------------
 
  //photo-nucleon system decaying to meson and baryon
  cout<<"rho "<<rho.Pdg()<<endl;
  auto two_body = elSpectro::JpacTwoBody{jpac_amp,{ rho },{proton_id} }; 
  //create eic electroproduction of X + proton
  auto elscat =   phoPro::Build_ep_Collision(ebeamE,pbeamE,two_body) ;
  elscat->SetLimit_Ymin(0.05);
  generator.Reaction( elscat );

  
  // ---------------------------------------------------------------------------
  // Initialize HepMC3
  // ---------------------------------------------------------------------------
  generator.SetWriter( new  elSpectro::HepMC3Writer{Form("/home/dglazier/elspectro_out/eic_rho_%d_%d.hepmc",(int)ebeamE,(int)pbeamE)});
  
  // ---------------------------------------------------------------------------
  //initilase the generator, may take some time for making distribution tables 
  // ---------------------------------------------------------------------------
  generator.InitGeneration();
  
  // ---------------------------------------------------------------------------
  //Set number of events via experimental luminosity and beamtime
  // ---------------------------------------------------------------------------
  //rho_prod->SetCombinedBranchingFraction(rho.BranchRatio()); 
  generator.SetNEvents_via_LuminosityTime(nLumi,24*60*60*nDays);
  // generator.SetNEvents(1000);
  //auto fastIntegral=generator.Reaction()->IntegrateCrossSectionFast();
  //std::cout<<"       check fast cross section "<<fastIntegral<<std::endl;

  
  // ---------------------------------------------------------------------------
  //Generate events
  // ---------------------------------------------------------------------------
  generator.processAll();
  


 }
