void ComparePhaseSpaceto2() {

 
   //(Momentum, Energy units are Gev/C, GeV)
   Double_t masses[2] = {0.139, 0.139} ;

   TGenPhaseSpace event;
   TLorentzVector W(0,0,0.5,sqrt(0.775*0.775+0.5*0.5));
   
   event.SetDecay(W, 2, masses);

   TH1F*  h2= new TH1F("RhoM","#rho Mass",1000,0.2,1.2);
   TH1F*  h2p= new TH1F("pip","#pi+ momentum",1000,0,1);
   
   gBenchmark->Start("tgenphasespace");
   for (Int_t n=0;n<1E7;n++) {
   
     Double_t weight = event.Generate();

     
     TLorentzVector *pPip    = event.GetDecay(0);
     TLorentzVector *pPim    = event.GetDecay(1);
     
     TLorentzVector pRho = *pPim + *pPip;
     
     h2->Fill(pRho.M());
     h2p->Fill(pPim->P());
     
   }
   gBenchmark->Stop("tgenphasespace");
   gBenchmark->Print("tgenphasespace");


   
   using namespace elSpectro;
   elSpectro::Manager::Instance();

    auto rho=dynamic_cast<DecayingParticle*>( particle(113,model(new PhaseSpaceDecay{{},{211,-211}})));
   generator().SetModelForMassPhaseSpace(rho->Model());
   
   rho->SetXYZ(0,0,0.5);

   auto pionp = rho->Model()->Product(0);
   auto pionm = rho->Model()->Product(1);

   TH1F* hel=new TH1F("RhoMelS","#rho Mass",1000,0.2,1.2);
   hel->SetLineColor(2);
   TH1F*  help= new TH1F("pipelS","#pi+ momentum",1000,0,1);
   help->SetLineColor(2);

   gBenchmark->Start("elspectro");
 
   for (Int_t n=0;n<1E7;n++) {
     rho->GenerateProducts();
     auto rho_pions=pionp->P4()+pionm->P4();
     hel->Fill(rho_pions.M());
     help->Fill(pionp->P4().P());
   }
   gBenchmark->Stop("elspectro");
   gBenchmark->Print("elspectro");

   hel->Draw();
   h2->Draw("same");
   new TCanvas;
   help->Draw();
   h2p->Draw("same");
 
}
