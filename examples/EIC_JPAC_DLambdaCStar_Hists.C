//Just need jpacPhoto headers
#include "reaction_kinematics.hpp"
#include "regge_trajectory.hpp"
#include "amplitudes/vector_exchange.hpp"
#include "amplitudes/dirac_exchange.hpp"
#include "amplitudes/amplitude_sum.hpp"

#include "FunctionsForGenvector.h"

// ---------------------------------------------------------------------------
// Diagnostic histograms
// ---------------------------------------------------------------------------

TH1F hQ2("Q2","Q2",1000,0,100);
TH1D hgPhi("gPhi","gPhi",100,-180,180);
TH1D heE("eE","eE",1000,0,50);
TH2D heThPhi("eThPhi","eThPhi",50,-180,180,500,178,180);
TH1D heTh("eTh","eTh",1000,0,180);
TH1D hPTh("NTh","NTh",1000,0,180);
TH2D hPThPhi("PThPhi","PThPhi",50,-180,180,50,0,10);
TH1D hPPhi("NPh","NPh",1000,-180,180);
TH1D hePhi("ePh","ePh",90,-180,180);
TH1F hW("W","W",1000,0,40);
TH1F ht("t","t",1000,0,10);
TH1F hgE("gE","gE",1000,0,20);

TH1D hDMass("DMass","DMass",100,1.5,2.5);
TH1D hLMass("LMass","LMass",100,2,3);
TH1D hLStarMass("LStarMass","LMass",100,2.5,3);
TH1D hLStarMassCut("LStarMassCut","LMassCut",100,2.5,3);
TH1D hDKTh("DKTh","DKTh",1000,0,180);
TH1D hDKPh("DKPh","DKPh",1000,-180,180);
TH1D hDKP("DKP","DKP",1000,0,50);
TH1D hDpiTh("DpiTh","DpiTh",1000,0,180);
TH1D hDpiPh("DpiPh","DpiPh",1000,-180,180);
TH1D hDpiP("DKpiP","DpiP",1000,0,50);
TH1D hLKTh("LKTh","LKTh",1000,0,180);
TH1D hLKPh("LKPh","LKPh",1000,-180,180);
TH1D hLKP("LKP","LKP",1000,0,50);
TH1D hLpiTh("LpiTh","LpiTh",1000,0,180);
TH1D hLpiPh("LpiPh","LpiPh",1000,-180,180);
TH1D hLpiP("LKpiP","LpiP",1000,0,50);
TH1D hLprTh("LprTh","LprTh",1000,0,180);
TH1D hLprPh("LprPh","LprPh",1000,-180,180);
TH1D hLprP("LKprP","LprP",1000,0,50);

//based on $JPACPHOTO/jpacPhoto/executables/open_charm/open_charm_photoproduction.cpp
void EIC_JPAC_DLambdaCStar_Hists(double ebeamE = 5, double pbeamE = 41, double nLumi=1E33, int nDays = 11) {

  Double_t crossingAngle=0.0; //30mrad
  //define e- beam, pdg =11
  auto elBeam = initial(11,ebeamE);
  auto elBeamP4=elBeam->GetInteracting4Vector();
  elBeam->SetAngleThetaPhi(TMath::Pi()-crossingAngle,0);

  //define pr beam, pdg =2212
  auto prBeam = initial(2212,pbeamE);
  auto prBeamP4=prBeam->GetInteracting4Vector();
  prBeam->SetAngleThetaPhi(crossingAngle,0);
  
  // ---------------------------------------------------------------------------
  // AMPLITUDES
  // ---------------------------------------------------------------------------

   // Form factor parameter
    double eta = 1.;
    double lambdaQCD = 0.25;
    bool du_result = false;

    // ---------------------------------------------------------------------------
    // D phototproduction
    // ---------------------------------------------------------------------------

    // Set up Kinematics for Dbar LambdaC in final state
    reaction_kinematics kD(M_D, M_LAMBDAC, M_PROTON);
    kD.set_JP(PSEUDO_SCALAR);

    vector_exchange d_dstarEx(&kD, M_DSTAR, "D^{*} exchange");
    d_dstarEx.set_params({0.134, -13.2, 0.});
    d_dstarEx.set_formfactor(2, M_DSTAR + eta * lambdaQCD);
    d_dstarEx.set_debug(1);

    dirac_exchange d_lamcEx(&kD, M_LAMBDAC, "#Lambda_{c} exchange");
    d_lamcEx.set_params({sqrt(4.* PI * ALPHA), -4.3});
    d_lamcEx.set_formfactor(2, M_LAMBDAC + eta * lambdaQCD);

    amplitude_sum d_sum(&kD,  {&d_dstarEx, &d_lamcEx}, "Sum");



   // ---------------------------------------------------------------------------
  // elSpectro
  // ---------------------------------------------------------------------------
  //create a D decaying to K-pi+ (3.95% branch)
  auto D=particle(411,model(new PhaseSpaceDecay({},{-321,211})));

  //Create Lambda_c+ decaying to pK-pi+ (6.28% branch)
  auto Lc=particle(4122,model(new PhaseSpaceDecay{{},{2212,-321,211}}));

  //Add 5 LC* state mass distributions
  mass_distribution(9995,new DistTF1{TF1("hh",Form("TMath::BreitWigner(x,2.595,0.003)+TMath::BreitWigner(x,2.625,0.001)+TMath::BreitWigner(x,2.86,0.068)+TMath::BreitWigner(x,2.88,0.006)+TMath::BreitWigner(x,2.94,0.020)"),2.500,3.00)});
  auto LcStar=particle(9995,model(new PhaseSpaceDecay{{Lc},{211,-211}}));

  //create eic electroproduction of X + proton
  auto pGammaStarDecay = JpacModelst{&d_sum, {D,LcStar},{} }; //photo-nucleon system
  auto photoprod = DecayModelQ2W{0,&pGammaStarDecay,new TwoBody_stu{1., 0, 1.}};

  //combine beam, target and reaction products
  auto production=eic( elBeam,prBeam,&photoprod );

  std::cout<<"Electron lab theta "<<elBeamP4->Theta()*TMath::RadToDeg()<<endl;
  // ---------------------------------------------------------------------------
  // Initialize HepMC3
  // ---------------------------------------------------------------------------
  writer(new HepMC3Writer{Form("outCollision/jpac_DLcStar_%d_%d.txt",(int)ebeamE,(int)pbeamE)});
  // exit(0);
  
  // ---------------------------------------------------------------------------
  //initilase the generator, may take some time for making distribution tables 
  // ---------------------------------------------------------------------------
  initGenerator();
  //generator().SuppressPhaseSpace(0.0001);
  // ---------------------------------------------------------------------------
  //Set number of events via experimental luminosity and beamtime
  // ---------------------------------------------------------------------------
  production->SetCombinedBranchingFraction(0.0395*0.0628*0.1); //D Lambda_c +10% for lower production??..??
  generator().SetNEvents_via_LuminosityTime(nLumi,24*60*60*nDays);
  // generator().SetNEvents(20000);
  // auto fastIntegral=production->IntegrateCrossSectionFast();
  //std::cout<<"       check fast cross section "<<fastIntegral<<std::endl;

  
  // exit(0);
  // ---------------------------------------------------------------------------
  // Get event particles for filling histograms
  // ---------------------------------------------------------------------------
  auto D_K = static_cast<DecayingParticle*>(D)->Model()->Product(0);
  auto D_pi = static_cast<DecayingParticle*>(D)->Model()->Product(1);
  auto L_prK = static_cast<DecayingParticle*>(Lc)->Model()->Product(0);
  auto L_pi = static_cast<DecayingParticle*>(Lc)->Model()->Product(2);
  auto L_pr = static_cast<const DecayingParticle*>(L_prK)->Model()->Product(0);
  auto L_K = static_cast<const DecayingParticle*>(L_prK)->Model()->Product(1);
  
  auto LStar_Lcpip = static_cast<DecayingParticle*>(LcStar)->Model()->Product(0);

  
  auto LStar_pip = static_cast<const DecayingParticle*>(LStar_Lcpip)->Model()->Product(1);
  auto LStar_pim = static_cast<const DecayingParticle*>(LcStar)->Model()->Product(2);

  auto electron = photoprod.GetScatteredElectron();
 
  // ---------------------------------------------------------------------------
  // Generate events
  // ---------------------------------------------------------------------------
  
  gBenchmark->Start("generator");//timer

  while(finishedGenerator()==false){
    nextEvent();
    countGenEvent();
    if(generator().GetNDone()%10==0) std::cout<<"event number "<<generator().GetNDone()<<std::endl;
    //  cout<<LStar_pim->Pdg()<<" "<<LStar_pip->Pdg()<<endl;
    auto photon = *elBeamP4 - electron->P4();
    double Q2 = -photon.M2();
    double W = (photon+*prBeamP4).M();
    hQ2.Fill(Q2);
    hW.Fill(W);
    hgE.Fill(photon.E());

    auto P4D = D_K->P4()+D_pi->P4();
    hDMass.Fill(P4D.M());

    auto P4L = L_pr->P4()+L_K->P4()+L_pi->P4();
    double t = -1*(L_pr->P4()- *prBeamP4).M2();
    ht.Fill(t);
    // cout<<"pr "<<L_pr->P4()<<" K "<<L_K->P4()<<" pi "<<L_pi->P4()<<endl;
    // cout<<"L "<<P4L<<" "<<P4L.M()<<" "<<Lc->P4()<<" "<<Lc->P4().M()<<endl;
    // cout<<"D "<<P4D<<" "<<P4D.M()<<" "<<D->P4()<<" "<<D->P4().M()<<endl;
    hLMass.Fill(P4L.M());
    auto P4LStar = P4L+LStar_pim->P4()+LStar_pip->P4();
    hLStarMass.Fill(P4LStar.M());
    if(L_pr->P4().Theta()>4*TMath::DegToRad()
       &&L_K->P4().Theta()>4*TMath::DegToRad()
       &&L_pi->P4().Theta()>4*TMath::DegToRad()
       &&LStar_pim->P4().Theta()>4*TMath::DegToRad()
       &&LStar_pip->P4().Theta()>4*TMath::DegToRad()
       )
      hLStarMassCut.Fill(P4LStar.M());
     // auto Total1 =  *elBeamP4+*prBeamP4;
    // auto Total2 =  P4D+P4L+electron->P4();
    // ////  auto Total3 =  D->P4()+Lc->P4()+electron->P4();
    // cout<<"Total "<<Total1<< " "<<Total2<< " "<<Total3<<endl;
    
    hDKP.Fill(D_K->P4().P());
    hDKTh.Fill(D_K->P4().Theta()*TMath::RadToDeg());
    hDKPh.Fill(D_K->P4().Phi()*TMath::RadToDeg());
    hDpiP.Fill(D_pi->P4().P());
    hDpiTh.Fill(D_pi->P4().Theta()*TMath::RadToDeg());
    hDpiPh.Fill(D_pi->P4().Phi()*TMath::RadToDeg());
    
    hLKP.Fill(L_K->P4().P());
    hLKTh.Fill(L_K->P4().Theta()*TMath::RadToDeg());
    hLKPh.Fill(L_K->P4().Phi()*TMath::RadToDeg());
    hLpiP.Fill(L_pi->P4().P());
    hLpiTh.Fill(L_pi->P4().Theta()*TMath::RadToDeg());
    hLpiPh.Fill(L_pi->P4().Phi()*TMath::RadToDeg());
    hLprP.Fill(L_pr->P4().P());
    hLprTh.Fill(L_pr->P4().Theta()*TMath::RadToDeg());
    hLprPh.Fill(L_pr->P4().Phi()*TMath::RadToDeg());

  }
  
  gBenchmark->Stop("generator");
  gBenchmark->Print("generator");
  
  // ---------------------------------------------------------------------------
  // Report generator statistics, can be used for optimising
  // ---------------------------------------------------------------------------
  
  generator().Summary();

  // delete ampZc;
  
  // ---------------------------------------------------------------------------
  // Write diagnostic histograms
  // ---------------------------------------------------------------------------
  TH1D *hWdist = new TH1D(*(TH1D*)gDirectory->FindObject("Wdist"));
 
  TFile *fout = TFile::Open(Form("outCollision/jpac_LcStarD.root"), "recreate");
  // generated event distributions
  hQ2.Write();
  hW.Write();
  ht.Write();
  hDKP.Write();
  hDKTh.Write();
  hDKPh.Write();
  hDpiP.Write();
  hDpiTh.Write();
  hDpiPh.Write();
    
  hLKP.Write();
  hLKTh.Write();
  hLKPh.Write();
  hLpiP.Write();
  hLpiTh.Write();
  hLpiPh.Write();
  hLprP.Write();
  hLprTh.Write();
  hLprPh.Write();
  hDMass.Write();
  hLMass.Write();
  hLStarMass.Write();
  fout->Close();
  
}
