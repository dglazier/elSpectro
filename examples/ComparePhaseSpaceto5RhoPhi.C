#include "FunctionsForKinematics.h"

void ComparePhaseSpaceto5RhoPhi() {

  
  TLorentzVector target(0.0, 0.0, 0.0, 0.938);
  TLorentzVector beam(0.0, 0.0, 10.4, 10.4);
  TLorentzVector W = beam + target;

 
  
  Double_t masses[5] = { 0.938, 0.139, 0.139,0.494,0.494} ;
  TGenPhaseSpace event;
  event.SetDecay(W, 5, masses);

   TH1F* h1rho = new TH1F("hrho","hrho", 100,0,1.5);
   TH1F* h1phi = new TH1F("hphi","hphi", 100,0.9,1.5);
   TH1F* h1X = new TH1F("hX","hrho", 100,0,50);
  
   
   auto bwRho = TF1("hh","TMath::BreitWigner(x,0.78,0.1)",0.3,3.5);
   auto bwPhi = TF1("hh","TMath::BreitWigner(x,1.19,0.05)",1,1.4);
 
   auto bwmaxRho=bwRho.GetMaximum();
   auto bwmaxPhi=bwPhi.GetMaximum();
   
   gBenchmark->Start("tgenphasespace");
   Long64_t counter=0;
   Long64_t Nevents=1000;
   
   TLorentzVector *pProton = event.GetDecay(0);
   
   TLorentzVector *pPip    = event.GetDecay(1);
   TLorentzVector *pPim    = event.GetDecay(2);
   TLorentzVector *pKp    = event.GetDecay(3);
   TLorentzVector *pKm    = event.GetDecay(4);
 
   for (Int_t n=0;;n++) {
   
     Double_t weight = event.Generate();
     
     TLorentzVector pRho = *pPip + *pPim;
     TLorentzVector pPhi = *pKp + *pKm;
     
   
     if(weight > gRandom->Uniform() ){//phase space weight
      
       if(bwRho.Eval(pRho.M()) > bwmaxRho*gRandom->Uniform() && pRho.M()>0.2){ //weight with rho
	 if(bwPhi.Eval(pPhi.M()) > bwmaxPhi*gRandom->Uniform() && pPhi.M()>1){ //weight with phi
  
	   h1rho->Fill(pRho.M());
	   h1phi->Fill(pPhi.M());
	   h1X->Fill((pPhi+pRho).M());
	   counter++;
	   if(counter>Nevents) break;
	 }
       }
     }
   }
   gBenchmark->Stop("tgenphasespace");
   gBenchmark->Print("tgenphasespace");
   TCanvas* can =new TCanvas();
   can->Divide(3,1);
   can->cd(1);
   h1rho->Draw();
   can->cd(2);
   h1phi->Draw();
   can->cd(3);
   h1X->Draw();

   
   using namespace elSpectro;
   elSpectro::Manager::Instance();
   
   TH1F* h1rho_el = new TH1F("hrho_el","hrho", 100,0,1.5);
   h1rho_el->SetLineColor(2);
   TH1F* h1phi_el = new TH1F("hphi_el","hphi", 100,0.9,1.5);
   h1phi_el->SetLineColor(2);
   TH1F* h1X_el = new TH1F("hX_el","hX", 100,0,50);
   h1X_el->SetLineColor(2);
 
   mass_distribution(113,new DistTF1{TF1("Mrho","TMath::BreitWigner(x,0.78,0.1)",0.2,3.5)});
   auto rho=dynamic_cast<DecayingParticle*>( particle(113,model(new PhaseSpaceDecay{{},{211,-211}})));
   
   mass_distribution(333,new DistTF1{TF1("Mphi","TMath::BreitWigner(x,1.19,0.05)",0.9,1.8)});
   auto phi=dynamic_cast<DecayingParticle*>( particle(333,model(new PhaseSpaceDecay{{},{321,-321}})));

   mass_distribution(9995,new DistTF1{TF1("MX","1",1.2,14)});//phase space mass
   auto X=dynamic_cast<DecayingParticle*>( particle(9995,model(new PhaseSpaceDecay{{rho,phi},{}})));

   auto pX=dynamic_cast<DecayingParticle*>( particle(-2211,model(new PhaseSpaceDecay{{X},{2212}})));
   generator().SetModelForMassPhaseSpace(pX->Model());
   
   pX->SetXYZT(W.X(),W.Y(),W.Z(),W.T());

   // pX->Print();
   
   auto pionp = rho->Model()->Products()[0];
   auto pionm = rho->Model()->Products()[1];
   auto Kp = phi->Model()->Products()[0];
   auto Km = phi->Model()->Products()[1];
  
   auto proton = pX->Model()->Products()[1];

   
   gBenchmark->Start("elspectro");
   int extra=1000;
 
   for (Int_t n=0;n<Nevents*extra;n++) {
     pX->GenerateProducts();

     h1rho_el->Fill(rho->P4().M());
     h1phi_el->Fill(phi->P4().M());
     h1X_el->Fill(X->P4().M());
      

   }
   h1rho_el->Scale(1./extra);
   h1phi_el->Scale(1./extra);
   h1X_el->Scale(1./extra);
   gBenchmark->Stop("elspectro");
   gBenchmark->Print("elspectro");

   
   can->cd(1);
   h1rho_el->DrawCopy("hist same");
   can->cd(2);
   h1phi_el->DrawCopy("hist same");
   can->cd(3);
   h1X_el->DrawCopy("hist same");
   cout<<"NOTE elSpectro ran "<<extra <<" times as many events"<<endl;
    
}
