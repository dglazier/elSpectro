#include "Interface.h"
#include "LorentzVector.h"
#include "DistTF1.h"
#include "FunctionsForElectronScattering.h"
#include "LundWriter.h"
#include "TwoBody_stu.h"

#include "reaction_kinematics.hpp"
#include "regge_trajectory.hpp"
#include "amplitudes/vector_exchange.hpp"
#include "amplitudes/pseudoscalar_exchange.hpp"
#include "amplitudes/pomeron_exchange.hpp"
#include "amplitudes/amplitude_sum.hpp"
#include "amplitudes/baryon_resonance.hpp"

#include <TBenchmark.h>
#include <TH1.h>
#include <TH2.h>
#include <TFile.h>

jpacPhoto::amplitude_sum Amplitude();

// ---------------------------------------------------------------------------
// Diagnostic histograms
// ---------------------------------------------------------------------------

double minMass = 0.2;
double maxMass = 4;
TH1F hQ2("Q2","Q2",1000,0,5);
TH2F hWQ2("WQ2","WQ2",100,0,5,100,0,10);
TH2F hWt("Wt","Wt",100,3.8,4.6,100,0,9);
TH1D heE("eE","eE",1000,0,20);
TH1D heTh("eTh","eTh",1000,0,18);
TH1D hYTh("YTh","YTh",1000,0,180);
TH1F hW("W","W",1000,1,5);
TH1F ht("t","t",1000,0,10);
TH1F hgE("gE","gE",1000,0,20);
TH1F hgTh("gTh","gTh",1000,0,180);
TH1F hXM("MesonM","; X to #pi#pi#pi Mass (GeV)",1000,minMass,maxMass);
TH1F hP1("P1","Pi1",100,0,10);
TH1F hP2("P2","Pi2",100,0,10);
TH1F hPm("Pm","Pi-",100,0,10);
TH1F  hVectorCosTh("CosThGJ","cos(#theta_{GJ})",100,-1,1);
TH1F  hVectorPhi("PhiGJ","#phi_{GJ}",100,-180,180);
TH1F  hScatPhi("Phi_Y","#phi_{Y}",100,-180,180);
TH2F  hVectorvScatPhi("PhiGJPhi_Y","#phi_{GJ} v #phi_{Y}",50,-180,180,50,-180,180);
TH2F  hVectorPhvVectorTh("ThPhiGJ","#phi_{GJ} v #theta_{GJ}",50,-1,1,50,-180,180);
TH2F  hVectorPhvVectorTh2("ThPhiGJ2","#phi_{GJ} v #theta_{GJ}",50,-1,1,50,-180,180);

void MesonEx_JpsiPenta(double ebeamE=10.4,int nEvents = 5e1) {
  gRandom->SetSeed(0);
  using namespace elSpectro;
  elSpectro::Manager::Instance();
  
  //mass_distribution(9995,new DistTF1{TF1("hh","1",0.,4)});
  elSpectro::LorentzVector elbeam(0,0,ebeamE,elSpectro::escat::E_el(ebeamE));
  elSpectro::LorentzVector prbeam(0,0,0,elSpectro::escat::M_pr());
  
  //create a V decaying to J/psi pi+pi-
  // mass_distribution(443,new DistTF1{TF1("hhjpsi","TMath::BreitWigner(x,3.0969160,0.0000929)",3,3.2)});

  /* //THIS PART FOR JUST PHASE SPACE
  auto jpsi=particle(443,model(new PhaseSpaceDecay({},{11,-11})));
   
  //decay of gamma* + p  to p + J/psi
  //depends on s and t
  auto pGammaStarDecay = static_cast<DecayModelst*>(model(new DecayModelst{ {jpsi},{2212} }));
  // 
  //create mesonex electroproduction of X + proton
  //TwoBody_stu{0.1, 0.9, 3 ,0,0} //0.1 strength  s distribution (flat angular dist.),  0.9 strength t distribution with slope b = 3
  mesonex( ebeamE ,  new DecayModelQ2W{0, pGammaStarDecay,new TwoBody_stu{0.1, 0.9, 2 , 0 , 0} });
  */
  
  //OR JPAC amplitude versionmodel(new VectorSDMEDecay{{},{11,-11}})
  auto jpsi=particle(443,model(new VectorSDMEDecay({},{11,-11})));
  auto jpac_amp = Amplitude();
  auto pGammaStarDecay = new JpacModelst{&jpac_amp, {jpsi},{2212} };
  auto photoprod = new DecayModelQ2W{0,pGammaStarDecay,new TwoBodyFlat{}};
  //  auto photoprod = new DecayModelQ2W{0,pGammaStarDecay,new TwoBody_stu{0.1, 0.9, 0.5,0,0}};
  //auto photoprod = new DecayModelQ2W{0,pGammaStarDecay,new TwoBody_stu{0.5, 0.0, 1,0,0}};
  mesonex( ebeamE ,  photoprod);

  


  auto production=dynamic_cast<ElectronScattering*>(generator().Reaction());
  // production->SetLimitTarRest_eThmin(1.5*TMath::DegToRad());
  //production->SetLimitTarRest_eThmax(5.5*TMath::DegToRad());
  // production->SetLimitTarRest_ePmin(0.4);
  //production->SetLimitTarRest_ePmax(6);
  //production->SetLimit_Q2max(0.01);
  // production->GiveZVertexDist( new DistConst(-3) );
  production->GiveZVertexDist( new DistUniform(-3-2.5,-3+2.5) );
  

  auto Jel = static_cast<DecayingParticle*>(jpsi)->Model()->Products()[1];
  auto Jpo = static_cast<DecayingParticle*>(jpsi)->Model()->Products()[0];
  
  auto proton = pGammaStarDecay->GetBaryon();
  auto electron = dynamic_cast<DecayModelQ2W*>(production->Model())->GetScatteredElectron();
  // ---------------------------------------------------------------------------
  // Initialize LUND
  // ---------------------------------------------------------------------------
  
  writer(new LundWriter{Form("out_mesonex/ep_to_pJpsi_%d.dat",(int)ebeamE)});
  
  //initilase the generator, may take some time for making distribution tables 
  initGenerator();
  
  // ---------------------------------------------------------------------------
  //Set number of events via experimental luminosity and beamtime
  // ---------------------------------------------------------------------------
  production->SetCombinedBranchingFraction(0.06); //Just Jpsi->e+e-
  generator().SetNEvents_via_LuminosityTime(1E35,24*60*60*80);
  //or can just do generator().SetNEvents(1E6);
  auto fastIntegral=production->IntegrateCrossSectionFast();
  std::cout<<" check fast cross section "<<fastIntegral<<std::endl;

  // ---------------------------------------------------------------------------
  // Generate events
  // ---------------------------------------------------------------------------
  
  gBenchmark->Start("e");

  while(finishedGenerator()==false){
    nextEvent();
    countGenEvent();
    if(generator().GetNDone()%1000==0) std::cout<<"event number "<<generator().GetNDone()<<std::endl;
    
    auto photon = elbeam - electron->P4();
    double W = (photon+prbeam).M();
    
    //fill diagnostic histograms
    double Q2 = -photon.M2();
    hWQ2.Fill(W,Q2);
    hQ2.Fill(Q2);
    hW.Fill(W);
    hgE.Fill(photon.E());
    
    auto elec = electron->P4();
    heTh.Fill(elec.Theta() *TMath::RadToDeg());
    heE.Fill(elec.E());

    auto Xrec=Jel->P4()+Jpo->P4();
    hXM.Fill(Xrec.M());
 
    hP1.Fill(Jel->P4().P());
    hPm.Fill(Jpo->P4().P());

    double t = -1*(proton->P4()-prbeam).M2();// + kine::t0(W,0,prbeam.M(),Xrec.M(),prbeam.M());
    ht.Fill(t);
    hWt.Fill(W,t);

        //SDME related angles
    MomentumVector elScatAngles;
    auto gStarN=(photon+prbeam);
    kine::electroCMDecay(&gStarN,&elbeam,&photon,&Xrec,&elScatAngles);
    hScatPhi.Fill(elScatAngles.Phi()*TMath::RadToDeg());
    MomentumVector vectorAngles;
    kine::mesonDecayGJ(&photon,&Xrec,&(proton->P4()),&Jel->P4(),&vectorAngles);
    hVectorCosTh.Fill(TMath::Cos(vectorAngles.Theta()));
    hVectorPhi.Fill(vectorAngles.Phi()*TMath::RadToDeg());
    hVectorvScatPhi.Fill(elScatAngles.Phi()*TMath::RadToDeg(),vectorAngles.Phi()*TMath::RadToDeg());
    hVectorPhvVectorTh.Fill(TMath::Cos(vectorAngles.Theta()),vectorAngles.Phi()*TMath::RadToDeg());
    if(W>4.4&&W<4.5)
    hVectorPhvVectorTh2.Fill(TMath::Cos(vectorAngles.Theta()),vectorAngles.Phi()*TMath::RadToDeg());
 
  }
  gBenchmark->Stop("e");
  gBenchmark->Print("e");
 
  generator().Summary();
  // internally stored histograms for total ep cross section
  TH1D *hWdist = (TH1D*)gDirectory->FindObject("Wdistlow");
  TH1D *hGenWdist = (TH1D*)gDirectory->FindObject("genWdist");

  TFile *fout = TFile::Open(Form("out_mesonex/ep_to_pJpsi_%d.root",(int)ebeamE), "recreate");
  // total ep cross section inputs
  if(hWdist)hWdist->Write();
  if(hGenWdist)hGenWdist->Write();
 
  // generated event distributions
  hWQ2.Write();
  hWt.Write();
  hQ2.Write();
  hW.Write();
  ht.Write();
  hgE.Write();
  heTh.Write();
  heE.Write();
  hXM.Write();
  fout->Close();
 

  
}

jpacPhoto::amplitude_sum Amplitude(){
  using namespace jpacPhoto;

  reaction_kinematics * ptr = new reaction_kinematics(3.0969160);//Jpsi mass
  ptr->set_JP(1, -1);
  
  // ---------------------------------------------------------------------------
  // S - CHANNEL

  // Two different pentaquarks
  // masses and widths from 2015 LHCb paper [2]
  auto P_c4450 =new baryon_resonance(ptr, 3, -1, 4.45, 0.040, "P_{c}(4450)");
  P_c4450->set_params({0.01, .7071}); // 2% branching fraction and equal photocouplings

  auto P_c4380 = new baryon_resonance(ptr, 5, +1, 4.38, 0.205, "P_{c}(4380)");
  P_c4380->set_params({0.01, .7071}); // 2% branching fraction and equal photocouplings

  // ---------------------------------------------------------------------------
  // T - CHANNEL

  // Set up pomeron trajectory
  // Best fit values from [1]
   auto alpha =new linear_trajectory(+1, 0.941, 0.364, "pomeron");

  // Create amplitude with kinematics and trajectory
  auto  background = new pomeron_exchange(ptr, alpha, false, "Background");

  // normalization and t-slope
  background->set_params({0.379, 0.12});

  // ---------------------------------------------------------------------------
  // SUM
  // ---------------------------------------------------------------------------
  // Incoherent sum of the s and t channels
  // amplitude_sum sum5q(ptr, {background}, "5q Sum");
    amplitude_sum sum5q(ptr, {background, P_c4450}, "5q Sum");
  //  amplitude_sum sum10q(ptr, {background, P_c4450, P_c4380}, "10q Sum");

  return sum5q;
}
