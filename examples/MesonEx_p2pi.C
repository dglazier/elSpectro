#include "Interface.h"
#include "LorentzVector.h"
#include "DistTF1.h"
#include "FunctionsForElectronScattering.h"
#include "LundWriter.h"
#include "TwoBody_stu.h"


#include <TBenchmark.h>
#include <TH1.h>
#include <TH2.h>
#include <TFile.h>

// ---------------------------------------------------------------------------
// Diagnostic histograms
// ---------------------------------------------------------------------------

double minMass = 0.2;
double maxMass = 3;
TH1F hQ2("Q2","Q2",1000,0,5);
TH2F hWQ2("WQ2","WQ2",100,0,5,100,0,10);
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

void MesonEx_p2pi(double ebeamE=10.4,int nEvents = 5e4) {

  using namespace elSpectro;
  elSpectro::Manager::Instance();
  
  //mass_distribution(9995,new DistTF1{TF1("hh","1",0.,4)});
   //define e- beam, pdg =11 momentum = _beamP
  auto elBeam = initial(11,ebeamE);
  auto elbeam=elBeam->GetInteracting4Vector();
  //proton target at rest
  auto prTarget= initial(2212,0);
  auto prbeam=prTarget->GetInteracting4Vector();

  //  elSpectro::LorentzVector elbeam(0,0,ebeamE,elSpectro::escat::E_el(ebeamE));
  //  elSpectro::LorentzVector prbeam(0,0,0,elSpectro::escat::M_pr());
  
   //flat pure phase sapce distribution =1 for particle id 9995
  //mass_distribution(9995,new DistTF1{TF1("hh","1",0.,(elbeam+prbeam).M())});
  //add a Breit-Wigner resonance for particle id 9995
  mass_distribution(9995,new DistTF1{TF1("hh","0.9*TMath::BreitWigner(x,0.78,0.149) + 0.1*TMath::BreitWigner(x,1.27,0.187)",0.,(*elbeam+*prbeam).M())});

  //produced meson decaying to pi+ pi- with mass distribution 9995
  auto X=static_cast<DecayingParticle*>(particle(9995,model(new PhaseSpaceDecay{{},{211,-211}})) );
  
  //decay of gamma* + p  to p + X
  //depends on s and t
  auto pGammaStarDecay = static_cast<DecayModelst*>(model(new DecayModelst{ {X},{2212} }));
  
  mesonex( elBeam,prTarget ,  new DecayModelQ2W{0, pGammaStarDecay,new TwoBody_stu{0.1, 0.9, 4 , 0 , 0}} );


  // 
  //create mesonex electroproduction of X + proton
  //TwoBody_stu{0.1, 0.9, 3 ,0,0} //0.1 strength  s distribution (flat angular dist.),  0.9 strength t distribution with slope b = 3
  // mesonex( ebeamE ,  new DecayModelQ2W{0, pGammaStarDecay,new TwoBody_stu{0, 1, 5 , 0 , 0} });
  // mesonex( ebeamE ,  new DecayModelQ2W{0, pGammaStarDecay,new TwoBody_stu{0.1, 0.9, 3 , 0 , 0} });
  auto production=dynamic_cast<ElectronScattering*>(generator().Reaction());
  production->SetLimitTarRest_eThmin(2*TMath::DegToRad());
  production->SetLimitTarRest_eThmax(5*TMath::DegToRad());
  production->SetLimitTarRest_ePmin(0.5);
  production->SetLimitTarRest_ePmax(6);

  

  auto pip1 = static_cast<DecayingParticle*>(X)->Model()->Product(1);
  auto pim = static_cast<DecayingParticle*>(X)->Model()->Product(0);
  
  auto proton = pGammaStarDecay->GetBaryon();
  auto electron = dynamic_cast<DecayModelQ2W*>(production->Model())->GetScatteredElectron();
  // ---------------------------------------------------------------------------
  // Initialize LUND with 10000 events per file
  // ---------------------------------------------------------------------------
  
  writer(new LundWriter{Form("out_mesonex/ep_to_pX2pi_%d.dat",(int)ebeamE),10000});
  
  //initilase the generator, may take some time for making distribution tables 
  initGenerator();
  
  // ---------------------------------------------------------------------------
  // Generate events
  // ---------------------------------------------------------------------------
  
  gBenchmark->Start("e");

  for(int i=0;i<nEvents;i++){
    // for(int i=0;i<0;i++){
    if(i%1000==0) std::cout<<"event number "<<i<<std::endl;
    nextEvent();
    
    auto photon = *elbeam - electron->P4();
    double W = (photon+*prbeam).M();
    
    //fill diagnostic histograms
    double Q2 = -photon.M2();
    hWQ2.Fill(W,Q2);
    hQ2.Fill(Q2);
    hW.Fill(W);
    hgE.Fill(photon.E());
    
    auto elec = electron->P4();
    heTh.Fill(elec.Theta() *TMath::RadToDeg());
    heE.Fill(elec.E());

    auto Xrec=pip1->P4()+pim->P4();
    hXM.Fill(Xrec.M());
 
    hP1.Fill(pip1->P4().P());
    hPm.Fill(pim->P4().P());

    double t = -1*(proton->P4()-*prbeam).M2();// + kine::t0(W,0,prbeam.M(),Xrec.M(),prbeam.M());
    ht.Fill(t);

  }
  gBenchmark->Stop("e");
  gBenchmark->Print("e");
 
  generator().Summary();
  // internally stored histograms for total ep cross section
  TH1D *hWdist = (TH1D*)gDirectory->FindObject("Wdist");
  TH1D *hGenWdist = (TH1D*)gDirectory->FindObject("genWdist");

  TFile *fout = TFile::Open(Form("out_mesonex/ep_to_pX2pi_%d_z.root",(int)ebeamE), "recreate");
  // total ep cross section inputs
  if(hWdist)hWdist->Write();
  if(hGenWdist)hGenWdist->Write();
 
  // generated event distributions
  hWQ2.Write();
  hQ2.Write();
  hW.Write();
  ht.Write();
  hgE.Write();
  heTh.Write();
  heE.Write();
  hXM.Write();
  fout->Close();
 

  
}
