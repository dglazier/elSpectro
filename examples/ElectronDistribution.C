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


double Frixione(Double_t *x,Double_t *p);

void ElectronDistribution(double ebeamE=10.4,int nEvents = 1E7) {

  using namespace elSpectro;
  elSpectro::Manager::Instance();
  

  generator().Reaction(new ElectronScattering(ebeamE,0,0,0,new ScatteredElectron_xy(ebeamE , escat::M_pr(),escat::M_pr() ),new PhaseSpaceDecay{{},{-2211,11}}));
  auto production=dynamic_cast<ElectronScattering*>(generator().Reaction());

  elSpectro::LorentzVector elbeam(0,0,ebeamE,elSpectro::escat::E_el(ebeamE));
  elSpectro::LorentzVector prbeam(0,0,0,elSpectro::escat::M_pr());
  
  auto electron =production->Model()->Products()[1];

  //  production->SetLimitTarRest_eThmin(10*TMath::DegToRad());
  // production->SetLimitTarRest_eThmax(25*TMath::DegToRad());
  //production->SetLimitTarRest_ePmin(1);
  production->SetLimitTarRest_ePmax(9);
  // production->SetLimit_Q2min(1);
  // production->SetLimit_Q2max(4);
  // production->SetLimit_Xmin(0.2);
  // production->SetLimit_Xmax(0.5);
  
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
  }
  gBenchmark->Stop("e");
  gBenchmark->Print("e");

  auto yGen=static_cast<TH1D*>(gDirectory->Get("ydist"));
  
  auto frixione = new TF1{"frixione",&Frixione,0.,1,1};
  frixione->SetParameter(0,ebeamE);
  frixione->SetNpx(yGen->GetNbinsX());
  frixione->SetLineColor(1);

  yGen->Scale(frixione->Eval(0.5)/yGen->Interpolate(0.5));

  yGen->Draw();
  frixione->Draw("same");

}
double Frixione(Double_t *x,Double_t *p){
  auto e0=p[0];//only parameter is electron beam energy
  auto y=x[0]; //e' energy has to be x-axis

  const double alpha = 1./137.;
  const double PI = TMath::Pi();
  
  // double y = Eg/Eb;
  double me = 0.00051;
  double Mp = 0.9383;
  double q2_max = -1 * me*me*y*y/(1 - y) ;
  //double q2_max = -1 ;
  double eg=e0*y;
  double q2_min = - 2*eg*Mp;//add this
  //double q2_min = - 4;//add this
  // if(q2_min<-4)q2_min = - 4;
   /////////////////////double q2_min = - 4*e0*e0*y;

  auto flux = alpha/(2*PI) * (2*me*me*y*(1/q2_max - 1/q2_min) + (1 + (1 - y)*(1-y))/y * log(q2_min/q2_max));
  
  if(flux==TMath::Infinity()) return 0;
  if(TMath::IsNaN(flux)) return 0;

  if(flux<0)return 0;
  
  return flux;
}
