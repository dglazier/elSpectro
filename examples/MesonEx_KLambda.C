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
TH1D heE("eE","eE",1000,0,20);
TH1D heTh("eTh","eTh",1000,0,180);
TH1D hYTh("YTh","YTh",1000,0,180);
TH1F hW("W","W",1000,1,4);
TH1F hV("Vertex","Vertex",100,0,20);
TH1F ht("t","t",1000,0,10);
TH1F hgE("gE","gE",1000,0,20);
TH1F hgTh("gTh","gTh",1000,0,180);
TH1F hLM("MesonM","; #Lambda to #pip Mass (GeV)",1000,1,1.2);
TH1F hPKp("PKp","K+",100,0,10);
TH1F hPKm("PKm","K-",100,0,10);
TH1F hPp("Pp","proton",100,0,10);

void MesonEx_KLambda(double ebeamE=10.4,int nEvents = 5e4) {

  using namespace elSpectro;
  elSpectro::Manager::Instance();
  
   //define e- beam, pdg =11 momentum = _beamP
  auto elBeam = initial(11,ebeamE);
  auto elbeam=elBeam->GetInteracting4Vector();
  //proton target at rest
  auto prTarget= initial(2212,0);
  auto prbeam=prTarget->GetInteracting4Vector();


  auto Lambda=static_cast<DecayingParticle*>(particle(3122,model(new PhaseSpaceDecay{{},{2212,-211}})) );
  
  //decay of pGamma* tp Lambda + K
  auto pGammaStarDecay = static_cast<DecayModelst*>(model(new DecayModelst{ {Lambda},{321} }));
  
  // 
  //create mesonex electroproduction of X + neutron
  //TwoBody_stu{0.1, 0.9, 3 ,0,0} //0.1 strength  s distribution (flat angular dist.),  0.9 strength t distribution with slope b = 3
  mesonex( elBeam,prTarget ,  new DecayModelQ2W{0, pGammaStarDecay,new TwoBody_stu{0.1, 0.9, 3 , 0 , 0}} );

  //give limits to the detected electron 
  auto production=dynamic_cast<ElectronScattering*>(generator().Reaction());
  production->SetLimitTarRest_eThmin(1.5*TMath::DegToRad());
  production->SetLimitTarRest_eThmax(5.5*TMath::DegToRad());
  production->SetLimitTarRest_ePmin(0.4);
  production->SetLimitTarRest_ePmax(6);

  //get pointers to produced particles fror diagnostic histos
  auto Km = Lambda->Model()->Product(1);
  auto proton = Lambda->Model()->Product(0);
  
  
  auto Kp = pGammaStarDecay->GetMeson();
  auto electron = dynamic_cast<DecayModelQ2W*>(production->Model())->GetScatteredElectron();
  
  // ---------------------------------------------------------------------------
  // Initialize LUND
  // ---------------------------------------------------------------------------
  
  writer(new LundWriter{Form("out_mesonex/ep_to_LambdaK_%d.dat",(int)ebeamE)});
  
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

    //fill diagnostic histograms
    auto photon = *elbeam - electron->P4();
    double Q2 = -photon.M2();
    double W = (photon+*prbeam).M();
    double t = -1*(proton->P4()-*prbeam).M2();
    hQ2.Fill(Q2);
    hW.Fill(W);
    ht.Fill(t);
    hgE.Fill(photon.E());
    
    auto elec = electron->P4();
    heTh.Fill(elec.Theta() *TMath::RadToDeg());
    heE.Fill(elec.E());

    auto Lrec=proton->P4()+Km->P4();
    hLM.Fill(Lrec.M());
 
    hPKp.Fill(Kp->P4().P());
    hPKm.Fill(Km->P4().P());
    hPp.Fill(proton->P4().P());

    hV.Fill(proton->VertexPosition()->P()/10);
  }
  gBenchmark->Stop("e");
  gBenchmark->Print("e");
 
 
  // internally stored histograms for total ep cross section
  TH1D *hWdist = (TH1D*)gDirectory->FindObject("Wdist");
  TH1D *hGenWdist = (TH1D*)gDirectory->FindObject("genWdist");

  TFile *fout = TFile::Open(Form("out_mesonex/ep_to_LambdaK_%d.root",(int)ebeamE), "recreate");
  // total ep cross section inputs
  if(hWdist)hWdist->Write();
  if(hGenWdist)hGenWdist->Write();
 
  // generated event distributions
  hQ2.Write();
  hW.Write();
  ht.Write();
  hgE.Write();
  heTh.Write();
  heE.Write();
  hLM.Write();
  hPKp.Write();
  hPKm.Write();
  hPp.Write();
  fout->Close();
 
   generator().Summary();

}
