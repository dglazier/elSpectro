void ComparePhaseSpaceto4PiRho() {

  
  TLorentzVector target(0.0, 0.0, 0.0, 0.938272);
  TLorentzVector beam(0.0, 0.0, 20.4, 20.4);
  TLorentzVector W = beam + target;

 
   //(Momentum, Energy units are Gev/C, GeV)
  Double_t masses[7] = { 0.938272, 0.13957, 0.13957,0.13957, 0.13957,0.13957, 0.13957};//0.13957, 0.13957} ;

   TGenPhaseSpace event;
   
   event.SetDecay(W, 7, masses);

   TH2F *h2 = new TH2F("h2","h2", 50,1.1,20.8, 50,1.1,20.8);
   TH2F *h2weight = new TH2F("h2weight","h2weight",100,0,10.2,100,0,0.5);
   TH1F* h1 = new TH1F("h1","h1", 100,0,5);
   TH1F* hMeson = new TH1F("hMeson","Meson Mass", 100,0,10);
   TH1F* h1un = new TH1F("h1un","h1un", 100,0,5);
   h1un->SetLineColor(3);
   
   auto bw = TF1("hh","TMath::BreitWigner(x,0.78,0.1)",0.0,5.5);
   //auto bw = TF1("hh","1",0.0,5.5);
   //auto bw = TF1("hh","1",0.3,10.5);
   auto bwmax=bw.GetMaximum();
   gBenchmark->Start("tgenphasespace");
   Long64_t counter=0;
   Long64_t Nevents=500;
   
   for (Int_t n=0;;n++) {
   
     Double_t weight = event.Generate();

     TLorentzVector *pProton = event.GetDecay(0);
     
     TLorentzVector *pPip    = event.GetDecay(1);
     TLorentzVector *pPim    = event.GetDecay(2);
     
     TLorentzVector pPPip = *pProton + *pPip;
     TLorentzVector pPPim = *pProton + *pPim;
     TLorentzVector pRho = *pPip + *pPim;

     if(weight > gRandom->Uniform() ){
       if(bw.Eval(pRho.M()) > bwmax*gRandom->Uniform()){
	 h2->Fill(pPPip.M2() ,pPPim.M2());
	 h1->Fill(pRho.M());
	 hMeson->Fill((W-*pProton).M());
	 if(counter++>Nevents) break;
       }
     }
     h1un->Fill(pRho.M(),0.008);
 
     h2weight->Fill(pRho.M(),weight);
   }
   gBenchmark->Stop("tgenphasespace");
   gBenchmark->Print("tgenphasespace");
   TCanvas* can =new TCanvas();
   can->Divide(2,1);
   can->cd(1);
   h1->Draw();
   h1un->DrawCopy("same l");
   can->cd(2);
   h2->Draw();

   new TCanvas();
   h2weight->Draw("col1");

   TCanvas* canMeson = new TCanvas();
   hMeson->Draw("");
   
   TH2F *hel2 = new TH2F("hel2","h2", 50,1.1,20.8, 50,1.1,20.8);
   hel2->SetMarkerColor(2);
   TH1F* hel1 = new TH1F("hel1","h1", 100,0,5);
   hel1->SetLineColor(2);
   TH1F* hMeson2 = new TH1F("hMeson2","Meson Mass", 100,0,10);
   hMeson2->SetLineColor(2);
   
   using namespace elSpectro;
   elSpectro::Manager::Instance();

   elSpectro::LorentzVector elW{W.X(),W.Y(),W.Z(),W.T()};
   
   mass_distribution(113,new DistTF1{TF1("hh","TMath::BreitWigner(x,0.78,0.1)",0.0,5.5)});
   mass_distribution(9995,new DistTF1{TF1("hh","1",0.0,15.32)});
   //  mass_distribution(113,new DistTF1{TF1("hh","1",0.0,10.5)});
   auto rho=dynamic_cast<DecayingParticle*>( particle(113,model(new PhaseSpaceDecay{{},{211,-211}})));
   //auto rho2=dynamic_cast<DecayingParticle*>( particle(113,model(new PhaseSpaceDecay{{},{2212,-211}})));
   //auto meson=dynamic_cast<DecayingParticle*>( particle(9995,model(new PhaseSpaceDecay{{rho,rho2},{}})));

   //Note the following different methods will give pure
   //phase space distributions
   //if 113 mass distribution is flat

   //auto prho=dynamic_cast<DecayingParticle*>( particle(-2211,model(new PhaseSpaceDecay{{meson},{211}}))); //this OK
   auto prho=dynamic_cast<DecayingParticle*>( particle(-2211,model(new PhaseSpaceDecay{{rho},{211,-211,211,-211,2212}}))); //this OK
   //   auto prho=dynamic_cast<DecayingParticle*>( particle(-2211,model(new PhaseSpaceDecay{{},{2212,211,-211,211,-211}})));//this OK
   generator().SetModelForMassPhaseSpace(prho->Model());
   
   prho->SetXYZT(W.X(),W.Y(),W.Z(),W.T());


   auto pionp = rho->Model()->Product(0);
   auto pionm = rho->Model()->Product(1);
   
   Particle* proton =nullptr;
   //Particle* pionp =nullptr;
   //Particle* pionm =nullptr;
    for(auto& pp:particles().StableParticles()){
     if(pp->Pdg()==2212){
       proton=pp;
     }
     /*if(pp->Pdg()==211){
       pionp=pp;
     }
     if(pp->Pdg()==-211){
       pionm=pp;
       }*/
   }
    
    prho->PostInit(nullptr);
    gBenchmark->Start("elspectro");

   double scale = 100;
   
   for (Int_t n=0;n<Nevents*scale;n++) {
     prho->GenerateProducts();
     auto p4rho =pionp->P4()+pionm->P4();
      hel1->Fill(p4rho.M());
      
      auto p4ppip= pionp->P4()+proton->P4();
      auto p4ppim= pionm->P4()+proton->P4();
      hel2->Fill(p4ppip.M2() ,p4ppim.M2());
      hMeson2->Fill((elW-proton->P4()).M());

   }
   
   gBenchmark->Stop("elspectro");
   gBenchmark->Print("elspectro");
   can->cd(1);
   hel1->Scale(1./scale);
   hel1->DrawCopy("same");
   h1un->Draw("same");
   can->cd(2);
   hel2->Scale(1./scale);
   hel2->Draw("same");

   canMeson->cd();
   hMeson2->Scale(1./scale);
   hMeson2->DrawCopy("same");

   new TCanvas();
   hel1->Sumw2();
   hel1->Divide(h1);
   hel1->Draw();
   
}
