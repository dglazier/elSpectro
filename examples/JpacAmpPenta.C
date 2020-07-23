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
#include "amplitudes/baryon_resonance.hpp"

#include <TBenchmark.h>
#include <TH1.h>
#include <TH2.h>

//declare hists here so can draw interactive
  TH1F hQ2("Q2","Q2",1000,0,10);
  TH1D heE("eE","eE",1000,0,100);
  TH1D heTh("eTh","eTh",1000,0,180);
  TH1F hW("W","W",1000,0,25);
  TH1F hWhad("Whad","W",1000,0,50);
  TH1F hgE("gE","gE",1000,-100,100);
  TH1F hJpsiM("JpsiM","#rho Mass",1000,0.2,4.2);
  TH1F hJpsiE("JpsiE","#rho Energy",100,0,40);
  TH1F hJpsiTh("JpsiTh","#rho #theta",1000,0,180);
  TH1F hPE("PE","proton Energy",100,0,100);
  TH1F hPTh("PTh","proton #theta",1000,0,180);
  TH1F ht("t","t",1000,00,50);
  TH2F hWt("Wt","Wt",100,00,25,100,0,10);
  TH1F hu("u","u",1000,00,50);

void JpacAmpPenta(){
  using namespace jpacPhoto;
  using namespace elSpectro;
  elSpectro::Manager::Instance();
  // ---------------------------------------------------------------------------
  // AMPLITUDES
  // ---------------------------------------------------------------------------

  // Set up Kinematics
  reaction_kinematics * ptr = new reaction_kinematics(mJpsi, "jpsi");
  // ---------------------------------------------------------------------------
  // S - CHANNEL

  // Two different pentaquarks
  // masses and widths from 2015 LHCb paper [2]
  baryon_resonance P_c4450(ptr, 3, -1, 4.45, 0.040, "P_{c}(4450)");
  P_c4450.set_params({1, .7071}); // 2% branching fraction and equal photocouplings
  //  P_c4450.set_params({0.01, .7071}); // 2% branching fraction and equal photocouplings

  baryon_resonance P_c4380(ptr, 5, +1, 4.38, 0.205, "P_{c}(4380)");
  P_c4380.set_params({0.2, .7071}); // 2% branching fraction and equal photocouplings
  //  P_c4380.set_params({0.01, .7071}); // 2% branching fraction and equal photocouplings

  // ---------------------------------------------------------------------------
  // T - CHANNEL

  // Set up pomeron trajectory
  // Best fit values from [1]
  linear_trajectory alpha(+1, 0.941, 0.364, "pomeron");

  // Create amplitude with kinematics and trajectory
  pomeron_exchange background(ptr, &alpha, "Background");

  // normalization and t-slope
  background.set_params({0.379, 0.12});

  // ---------------------------------------------------------------------------
  // SUM
  // ---------------------------------------------------------------------------
  // Incoherent sum of the s and t channels
  amplitude_sum sum5q(ptr, {&background, &P_c4450}, "5q Sum");
  amplitude_sum sum10q(ptr, {&background, &P_c4450, &P_c4380}, "10q Sum");

  //////////////////////////////////////////////////////////////////////
  //elSpectro generator
    
  auto jpsi=particle(443,model(new PhaseSpaceDecay({},{11,-11})));
 
  //construct elSpectro Jpac decay model {jpac_amp, {decaying particles},and/or{pdg codes}, s_strength, t_strength, t_slope} where the latter 3 arguments define the decay angle distribution envelope
  //here decay to jpsi+proton with flat (s-channel envelope)
  auto jpac = new JpacModelQ2W{&sum5q, {jpsi},{2212}, 1, 1, 1 }; 
   
  //give generator decay model for g*N system
  auto production=eic( 5 , 40, jpac  ); 
  elSpectro::LorentzVector elbeam(0,0,-5,elSpectro::escat::E_el(5));
  elSpectro::LorentzVector prbeam(0,0,40,elSpectro::escat::E_pr(40));
 
  auto elJ = static_cast<DecayingParticle*>(jpsi)->Model()->Products()[0];
  auto posJ = static_cast<DecayingParticle*>(jpsi)->Model()->Products()[1];

  auto electron = jpac->GetScatteredElectron();
  auto proton = jpac->GetDecayBaryon();
  
  gBenchmark->Start("e");

  //initiase the generator, may take some time for making distribution tables 
  production->InitGen();
 
  for(int i=0;i<100000;i++){
    if(i%100==0) std::cout<<"event number "<<i<<std::endl;
    production->GenerateProducts();
    
    auto photon = elbeam - electron->P4();
    hQ2.Fill(-photon.M2());
 
    hW.Fill((photon+prbeam).M());
    
    hgE.Fill(photon.E());
  
    heTh.Fill(electron->P4().Theta() *TMath::RadToDeg());
    heE.Fill(electron->P4().E());
    
    auto Jpsi=elJ->P4()+posJ->P4();
    
    hJpsiM.Fill(Jpsi.M());
    hJpsiTh.Fill(Jpsi.Theta()*TMath::RadToDeg());
    hJpsiE.Fill(Jpsi.E());

    hPTh.Fill(proton->P4().Theta()*TMath::RadToDeg());
    hPE.Fill(proton->P4().E());

    ht.Fill(-(prbeam-proton->P4()).M2());
    hWt.Fill((photon+prbeam).M(),-(prbeam-proton->P4()).M2());
    hu.Fill((photon-proton->P4()).M2());
    
    auto jpsiProton=proton->P4()+Jpsi;
    hWhad.Fill(jpsiProton.M());
  
    auto sum1 = jpsiProton+electron->P4();
    auto sum2 = elbeam+prbeam;

    //if(sum1!=sum2)cout<<"energy moemtnum not conserved! "<<sum1<<" "<<sum2<<endl;
  }
  gBenchmark->Stop("e");
  gBenchmark->Print("e");


  

}
