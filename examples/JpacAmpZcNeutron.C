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
#include <TFile.h>

#include "HepMC3/Units.h"
#include "HepMC3/GenEvent.h"
#include "HepMC3/GenParticle.h"
#include "HepMC3/GenVertex.h"
#include "HepMC3/Print.h"
#include "HepMC3/GenRunInfo.h"
#include "HepMC3/WriterAscii.h"

// ---------------------------------------------------------------------------
// Diagnostic histograms
// ---------------------------------------------------------------------------

double minMass = 3.5;
double maxMass = 4.3;
TH1F hQ2("Q2","Q2",1000,0,100);
TH1D heE("eE","eE",1000,0,20);
TH1D heTh("eTh","eTh",1000,0,180);
TH1F hW("W","W",1000,0,100);
TH1F ht("t","t",1000,0,10);
TH1F hgE("gE","gE",1000,0,20);
TH1F hMesonM("MesonM","; J/#psi#pi Mass (GeV)",1000,minMass,maxMass);
TH1F hJpsiM("JpsiM","; e+e- Mass (GeV)",1000,2.9,3.3);
TH2F hElePVsEta("ElePVsEta","; #eta; p (GeV)",200,-5,5,500,0,50);
TH2F hPosPVsEta("PosPVsEta","; #eta; p (GeV)",200,-5,5,500,0,50);
TH2F hPionPVsEta("PionPVsEta","; #eta; p (GeV)",200,-5,5,500,0,50);
TH2F hRecoilPVsEta("RecoilPVsEta","; #eta; p (GeV)",200,0,10,1000,0,275);
TH2F hRecoilThetaVsP("RecoilThetaVsP","; p (GeV); #theta (mrad)",1000,0,275,200,0,200);
TH1F hRecoilPt("RecoilPt","; p_{T} (GeV)",200,0,5.0);

void JpacAmpZcNeutron(double ebeamE = 5, double pbeamE = 41, int nEvents = 5e4) {

 using namespace jpacPhoto;
 using namespace elSpectro;
 using namespace HepMC3;
 elSpectro::Manager::Instance();

 // ---------------------------------------------------------------------------
 // AMPLITUDES
 // ---------------------------------------------------------------------------
 
 // For regge amps need pion trajectory
 linear_trajectory alpha(1, -0.7*mPi*mPi, 0.7, "pionic trajectory");

 // ---------------------------------------------------------------------------
 // ZC(3900)
 // ---------------------------------------------------------------------------
 
 // Kinematics for 3900 MeV vector
 reaction_kinematics * ptr = new reaction_kinematics(3.9, "Z_{c}^{+}(3900)");
 
 // Amplitudes
 pseudoscalar_exchange Z3900(ptr, mPi, "#pi exchange");
 pseudoscalar_exchange Z3900R(ptr, &alpha, "Reggeon exchange");
 
 // Couplings for 4 MeV width
 Z3900.set_params({0.67 * 3.90, sqrt(4.*M_PI*14.4)});
 Z3900R.set_params({0.67 * 3.90, sqrt(4.*M_PI*14.4)});

 //amplitude_sum sum(ptr, {&Z3900, &Z3900R}, "Sum");
 amplitude_sum sum(ptr, {&Z3900}, "Sum");

 // ---------------------------------------------------------------------------
 // elSpectro
 // ---------------------------------------------------------------------------

 //create a Zc decaying to J/psi pi+
 mass_distribution(9995,new DistTF1{TF1("hh","TMath::BreitWigner(x,3.9,0.05)",minMass,maxMass)});
 auto jpsi=particle(443,model(new PhaseSpaceDecay({},{11,-11})));
 auto zc=particle(9995,model(new PhaseSpaceDecay{{jpsi},{211}}));
 
 //create eic electroproduction of Zc + neutron
 auto jpac = new JpacModelQ2W{&sum, {zc},{2112}, 0, 1, 2.0};
 auto production=eic( ebeamE, pbeamE, jpac );
 elSpectro::LorentzVector elbeam(0,0,-1*ebeamE,elSpectro::escat::E_el(ebeamE));
 elSpectro::LorentzVector prbeam(0,0,pbeamE,elSpectro::escat::E_pr(pbeamE));

 auto pionp = static_cast<DecayingParticle*>(zc)->Model()->Products()[1];
 auto posJ = static_cast<DecayingParticle*>(jpsi)->Model()->Products()[0];
 auto eleJ = static_cast<DecayingParticle*>(jpsi)->Model()->Products()[1];

 auto electron = jpac->GetScatteredElectron();
 auto neutron = jpac->GetDecayBaryon();

 gBenchmark->Start("e");

 //initialize the generator, may take some time for making distribution tables 
 production->InitGen();

 // compute total ep cross section from JPAC sigma_gp(s) and f(s)
 // need 2 histograms... 
 // f(s) genWdist
 // sigma(s) Wdist

 // ---------------------------------------------------------------------------
 // Initialize HepMC3
 // ---------------------------------------------------------------------------
 
 shared_ptr<GenRunInfo> run = std::make_shared<GenRunInfo>();
 Writer *ascii_io = new WriterAscii(Form("out/zc_piexch_elspectro_%d_%d.txt",(int)ebeamE,(int)pbeamE),run);
 
 // ---------------------------------------------------------------------------
 // Generate events
 // ---------------------------------------------------------------------------
 
 gBenchmark->Start("e");
 for(int i=0;i<nEvents;i++){
   if(i%100==0) std::cout<<"event number "<<i<<std::endl;
   production->GenerateProducts();
   
   auto photon = elbeam - electron->P4();
   double Q2 = -photon.M2();
   double W = (photon+prbeam).M();
   double t = -1*(neutron->P4()-prbeam).M2();
   hQ2.Fill(Q2);
   hW.Fill(W);
   ht.Fill(t);
   hgE.Fill(photon.E());
   
   //LorentzVector totalIn = elbeam + prbeam;
   //LorentzVector totalOut = electron->P4() + neutron->P4() + eleJ->P4() + posJ->P4() + pionp->P4();
   //cout<<totalIn.E()<<" "<<totalOut.E()<<endl;

   auto elec = electron->P4();
   heTh.Fill(elec.Theta() *TMath::RadToDeg());
   heE.Fill(elec.E());
   
   auto jpsiP4 = eleJ->P4() + posJ->P4();
   hJpsiM.Fill(jpsiP4.M());
   auto jpsi_pion=pionp->P4()+jpsiP4;
   hMesonM.Fill(jpsi_pion.M());
   
   hElePVsEta.Fill(eleJ->P4().Eta(), eleJ->P4().P());
   hPosPVsEta.Fill(posJ->P4().Eta(), posJ->P4().P());
   hPionPVsEta.Fill(pionp->P4().Eta(), pionp->P4().P());
   hRecoilPVsEta.Fill(neutron->P4().Eta(), neutron->P4().P());
   hRecoilThetaVsP.Fill(neutron->P4().P(), neutron->P4().Theta()*1000.);
   hRecoilPt.Fill(neutron->P4().Pt());

   // ---------------------------------------------------------------------------
   // HepMC Writer beam particles
   // ---------------------------------------------------------------------------
   GenEvent *evt = new GenEvent(Units::GEV,Units::MM);
   const int nParticles = 7;
   GenParticlePtr particlePtr[nParticles];
   GenVertexPtr v1 = std::make_shared<GenVertex>();
   particlePtr[0] = std::make_shared<GenParticle>( FourVector( elbeam.X(), elbeam.Y(), elbeam.Z(), elbeam.E() ), 11, 3 );
   particlePtr[1] = std::make_shared<GenParticle>( FourVector( prbeam.X(), prbeam.Y(), prbeam.Z(), prbeam.E() ), 2212, 3 );
   for(int i=0; i<2; i++) v1->add_particle_in (particlePtr[i]);
   
   // ---------------------------------------------------------------------------
   // HepMC Writer final state particles (need to generalize)
   // ---------------------------------------------------------------------------
   particlePtr[2] = std::make_shared<GenParticle>( FourVector( electron->P4().X(), electron->P4().Y(), electron->P4().Z(), electron->P4().E() ), electron->Pdg(), 1 );
   particlePtr[3] = std::make_shared<GenParticle>( FourVector( neutron->P4().X(), neutron->P4().Y(), neutron->P4().Z(), neutron->P4().E() ), neutron->Pdg(), 1 );
   particlePtr[4] = std::make_shared<GenParticle>( FourVector( pionp->P4().X(), pionp->P4().Y(), pionp->P4().Z(), pionp->P4().E() ), pionp->Pdg(), 1 );
   particlePtr[5] = std::make_shared<GenParticle>( FourVector( posJ->P4().X(), posJ->P4().Y(), posJ->P4().Z(), posJ->P4().E() ), posJ->Pdg(), 1 );
   particlePtr[6] = std::make_shared<GenParticle>( FourVector( eleJ->P4().X(), eleJ->P4().Y(), eleJ->P4().Z(), eleJ->P4().E() ), eleJ->Pdg(), 1 );
   for(int i=2; i<nParticles; i++) v1->add_particle_out (particlePtr[i]);

   evt->add_vertex(v1);
   
   evt->add_attribute("Q2",make_shared<DoubleAttribute>(Q2));
   evt->add_attribute("W",make_shared<DoubleAttribute>(W));
   ascii_io->write_event(*evt);

   delete evt;
 }
 gBenchmark->Stop("e");
 gBenchmark->Print("e");
 
 delete ascii_io;

 // internally stored histograms for total ep cross section
 TH1D *hWdist = (TH1D*)gDirectory->FindObjectAny("Wdist"); // sigma_gp(W) from JPAC
 TH1D *hGenWdist = (TH1D*)gDirectory->FindObjectAny("genWdist"); // f(W) photon flux distribution

 TFile *fout = TFile::Open(Form("out/zc_elspectro_%d_%d_diagnostic.root",(int)ebeamE,(int)pbeamE), "recreate");
 // total ep cross section inputs
 hWdist->Write();
 hGenWdist->Write();
 // generated event distributions
 hQ2.Write();
 hW.Write();
 ht.Write();
 hgE.Write();
 heTh.Write();
 heE.Write();
 hJpsiM.Write();
 hMesonM.Write();
 hElePVsEta.Write();
 hPosPVsEta.Write();
 hPionPVsEta.Write();
 hRecoilPVsEta.Write();
 hRecoilThetaVsP.Write();
 hRecoilPt.Write();
 fout->Close();

 // compute total ep cross section from internally stored histograms
 hGenWdist->Scale(1.0/hGenWdist->Integral()); // normalize flux PDF
 double integrated_xsection = 0; // get sigma_ep from integral over W: f(W)*sigma_gp(W)
 for(int i=0; i<hWdist->GetNbinsX(); i++) {
   double W = hWdist->GetXaxis()->GetBinCenter(i+1);
   double WbinWidthScale = hWdist->GetBinWidth(i+1);
   double W_xsection = hWdist->GetBinContent(i+1);
   
   int genWbin = hGenWdist->GetXaxis()->FindBin(W);
   WbinWidthScale /= hGenWdist->GetXaxis()->GetBinWidth(genWbin);
   double W_fluxWeight = hGenWdist->GetBinContent(genWbin) * WbinWidthScale;
   integrated_xsection += W_xsection * W_fluxWeight;
 }
 cout<<"Total cross section ("<<(int)ebeamE<<","<<(int)pbeamE<<"): sigma_ep = "<<integrated_xsection<<" nb "<<endl;

}
