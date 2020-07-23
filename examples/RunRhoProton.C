{

 using namespace elSpectro;
 elSpectro::Manager::Instance();

 //create a rho decaying to pi+,pi-
 mass_distribution(113,new DistTF1{TF1("hh","TMath::BreitWigner(x,0.78,0.1)",0.2,2)});
 auto rho=particle(113,model(new PhaseSpaceDecay{{},{211,-211}}));

 //create eic electroproduction of rho + proton
 auto pGammaStarDecay = model(new PhaseSpaceDecay{ {rho},{2212} });
 
 auto production=eic( 10 , 100, new DecayModelQ2W{1.2, pGammaStarDecay }  ); //threshold of W = 1.2
 generator().SetModelForMassPhaseSpace(pGammaStarDecay);
 
 auto products = production->Model()->Products();
 
 auto pionp = static_cast<DecayingParticle*>(rho)->Model()->Products()[0];
 auto pionm = static_cast<DecayingParticle*>(rho)->Model()->Products()[1];

 Particle* proton =nullptr;
 for(auto* pp:particles().StableParticles()){
   if(pp->Pdg()==2212){
     proton=pp;
     break;
   }
 }
 
 TH1F hQ2("Q2","Q2",1000,0,200);
 TH1D heE("eE","eE",1000,0,100);
 TH1D heTh("eTh","eTh",1000,0,180);
 TH1F hW("W","W",1000,0,50);
 TH1F hWhad("Whad","W",1000,0,50);
 TH1F hgE("gE","gE",1000,0,100);
 TH1F hRhoM("RhoM","#rho Mass",1000,0.2,1.2);
 TH1F h2RhoM("Rho2M","#rho Mass",1000,0.2,1.2);
 TH1F hRhoE("RhoE","#rho Energy",100,0,40);
 TH1F hRhoTh("RhoTh","#rho #theta",1000,0,180);
 TH1F hPE("PE","proton Energy",100,0,100);
 TH1F hPTh("PTh","proton #theta",1000,0,180);
 
 gBenchmark->Start("e");
 production->InitGen();
 for(int i=0;i<1E6;i++){
    production->GenerateProducts();
  
    auto photonNuc = products[0]->P4();
    hW.Fill(photonNuc.M());
    photonNuc.SetE(photonNuc.T()- 0.93827208816);
    hQ2.Fill(-photonNuc.M2());
    hgE.Fill(photonNuc.T()-photonNuc.M());
  
    auto elec = products[1]->P4();
    heTh.Fill(elec.Theta() *TMath::RadToDeg());
    heE.Fill(elec.E());
    
    auto rho_pions=pionp->P4()+pionm->P4();
    hRhoM.Fill(rho_pions.M());
    hRhoTh.Fill(rho_pions.Theta()*TMath::RadToDeg());
    hRhoE.Fill(rho_pions.E());

    hPTh.Fill(proton->P4().Theta()*TMath::RadToDeg());
    hPE.Fill(proton->P4().E());
  
    auto rhoProton=proton->P4()+rho_pions;
    hWhad.Fill(rhoProton.M());
    if(rhoProton.M()<100&&rhoProton.M()>0)h2RhoM.Fill(rho_pions.M());
    
    }
  gBenchmark->Stop("e");
  gBenchmark->Print("e");

}
