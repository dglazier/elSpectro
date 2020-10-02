#include "Interface.h"
#include "LorentzVector.h"
#include "DecayGammaN_Test.h"
#include "DistTF1.h"
#include "JpacModelQ2W.h"
#include "FunctionsForElectronScattering.h"

#include "reaction_kinematics.hpp"
#include "regge_trajectory.hpp"
#include "amplitudes/pseudoscalar_exchange.hpp"
#include "amplitudes/pomeron_exchange.hpp"
#include "amplitudes/amplitude_sum.hpp"

#include <TBenchmark.h>
#include <TH1.h>

//declare hists here so can draw interactive
  TH1F hQ2("Q2","Q2",1000,0,10);
  TH1D heE("eE","eE",1000,0,100);
  TH1D heTh("eTh","eTh",1000,0,180);
  TH1F hW("W","W",1000,0,50);
  TH1F hWhad("Whad","W",1000,0,50);
  TH1F hgE("gE","gE",1000,-100,100);
  TH1F hRhoM("RhoM","#rho Mass",1000,0.2,1.2);
  TH1F h2RhoM("Rho2M","#rho Mass",1000,0.2,1.2);
  TH1F hRhoE("RhoE","#rho Energy",100,0,40);
  TH1F hRhoTh("RhoTh","#rho #theta",1000,0,180);
  TH1F hPE("PE","proton Energy",100,0,100);
  TH1F hPTh("PTh","proton #theta",1000,0,180);
  TH1F ht("t","t",1000,00,50);
  TH1F hu("u","u",1000,00,50);

void JpacAmpRho(){
  using namespace jpacPhoto;
  using namespace elSpectro;
  elSpectro::Manager::Instance();

  // Kinematics for 780 MeV vector
  double rhoMass=0.77549000;
  double rhoWidth=0.15100000;
  auto* ptrRho = new jpacPhoto::reaction_kinematics(rhoMass, "rho0");
  // Amplitudes
  // For regge amps need pion trajectory
  linear_trajectory alpha(1, -0.7*mPi*mPi, 0.7, "pionic trajectory");
  pseudoscalar_exchange rho_pi(ptrRho, mPi, "#pi");
  pseudoscalar_exchange rho_regge(ptrRho, &alpha, "Reggeon exchange");
  // Couplings for rho  width
  rho_pi.set_params({5, sqrt(rhoWidth*1000.*M_PI*14.4)});
  rho_regge.set_params({4, sqrt(rhoWidth*1000.*M_PI*14.4)});

  // Best fit values from [1] from high energy
  linear_trajectory alpha2016(1, 1.1, 0.11, "pomeron");
  pomeron_exchange rho_pomeron(ptrRho, &alpha2016, false, "pomeron");
  rho_pomeron.set_params({30,1.});

 
  amplitude_sum sum(ptrRho, {&rho_pi, &rho_pomeron}, "Sum");

  //////////////////////////////////////////////////////////////////////
  //elSpectro generator
  
  //create a rho decaying to pi+,pi-
  mass_distribution(113,new DistTF1{
      TF1("hh",Form("TMath::BreitWigner(x,%f,%f)",rhoMass,rhoWidth),0.2,3)}
    );
  
  auto rho=particle(113,model(new PhaseSpaceDecay{{},{211,-211}}));

  //construct elSpectro Jpac decay model {jpac_amp, {decaying particles},and/or{pdg codes}, s_strength, t_strength, t_slope} where the latter 3 arguments define the decay angle distribution envelope
  auto jpac = new JpacModelQ2W{&sum, {rho},{2212}, 1, 0,1 };
   
  //give generator decay model for g*N system
  auto production=eic( 10 , 100, jpac  ); 
  elSpectro::LorentzVector elbeam(0,0,-10,elSpectro::escat::E_el(10));
  elSpectro::LorentzVector prbeam(0,0,100,elSpectro::escat::E_pr(100));
 
  auto pionp = static_cast<DecayingParticle*>(rho)->Model()->Products()[0];
  auto pionm = static_cast<DecayingParticle*>(rho)->Model()->Products()[1];

  auto electron = jpac->GetScatteredElectron();
  auto proton = jpac->GetDecayBaryon();
  
  gBenchmark->Start("e");

  //initiase the generator, may take some time for making distribution tables 
  production->InitGen();
  for(int i=0;i<10000;i++){
    if(i%100==0) std::cout<<"event number "<<i<<std::endl;
    production->GenerateProducts();
    
    auto photon = elbeam - electron->P4();
    hQ2.Fill(-photon.M2());
 
    hW.Fill((photon+prbeam).M());
    
    hgE.Fill(photon.E());
  
    heTh.Fill(electron->P4().Theta() *TMath::RadToDeg());
    heE.Fill(electron->P4().E());
    
    auto rho_pions=pionp->P4()+pionm->P4();
    hRhoM.Fill(rho_pions.M());
    hRhoTh.Fill(rho_pions.Theta()*TMath::RadToDeg());
    hRhoE.Fill(rho_pions.E());

    hPTh.Fill(proton->P4().Theta()*TMath::RadToDeg());
    hPE.Fill(proton->P4().E());

    ht.Fill(-(prbeam-proton->P4()).M2());
    hu.Fill((photon-proton->P4()).M2());
    
    auto rhoProton=proton->P4()+rho_pions;
    hWhad.Fill(rhoProton.M());
    if(rhoProton.M()<100&&rhoProton.M()>0)h2RhoM.Fill(rho_pions.M());
 
    auto sum1 = rhoProton+electron->P4();
    auto sum2 = elbeam+prbeam;

    //if(sum1!=sum2)cout<<"energy moemtnum not conserved! "<<sum1<<" "<<sum2<<endl;
  }
  gBenchmark->Stop("e");
  gBenchmark->Print("e");


  

}
