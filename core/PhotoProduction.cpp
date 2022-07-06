#include "PhotoProduction.h"
#include "FunctionsForGenvector.h"
#include "Manager.h"
#include "Interface.h" //for generator
#include "DecayModelW.h"
#include "DontDecay.h"
#include "BremstrPhoton.h"
#include "Bremsstrahlung.h"
#include "Distribution.h"
#include "ScatteredElectron_xy.h"
#include <TDatabasePDG.h>
#include <TH1F.h>
#include <TFile.h>
#include <TBenchmark.h>

#include <Math/Functor.h>
#include <RooFunctorBinding.h>
#include <RooRealVar.h>
#include <RooArgList.h>

namespace elSpectro{
  int PhotoProduction::NintegralsPhotoProduction=0;

 /////////////////////////////////////////////////////////////////////
  PhotoProduction::PhotoProduction(CollidingParticle *photon,CollidingParticle* target,  DecayModel* model):
    ProductionProcess{photon,target,model},
    _beamPhot{22},
    _beamNucl{target->GetInteractingPdg()}
  {
    //For convenience keep our own
    //pointers to photon and target
    _photonptr = photon;
    _targetptr = target;
    
    SetNominalBeamCondtion();
  }

  /////////////////////////////////////////////////////////////////////
  void PhotoProduction::SetNominalBeamCondtion(){

    auto bremPhot =dynamic_cast<const BremstrPhoton*>(_photonptr->Decayer());
    auto bremModel =dynamic_cast<Bremsstrahlung*>(_photonptr->Model());
    bremModel->Products()[0]->SetXYZT(0,0,bremPhot->GetBeamEnergy(),bremPhot->GetBeamEnergy());
    
    _reactionInfo._ebeam = _photonptr->GetNominal4Vector();
 
    _beamPhot.SetP4(*(_photonptr->GetInteracting4Vector()));
    
    _beamNucl.SetP4(*(_targetptr->GetInteracting4Vector()));
    _massIon=_beamNucl.PdgMass();


    //For nominal beam conditions set nucleon rest frame vectors
    //in case info is needed in PostInit stage
    //Boost into ion rest frame
    auto prBoost=_beamNucl.P4().BoostToCM();
    _nuclRestNucl=LorentzVector(0,0,0,_beamNucl.Mass());
    _nuclRestPhot= boost(_beamPhot.P4(),prBoost);

    //set inital lab particles
    //this can be written to output file
    std::cout<<"PhotoProduction vertex "<<_photonptr->VertexPosition()<<std::endl;
    AddInitialParticlePtr(_photonptr);
    AddInitialParticlePtr(_targetptr);
 

    std::cout<<" PhotoProduction::SetNominalBeamCondtion() photon lab "<<_beamPhot.P4()<<std::endl;
   std::cout<<" PhotoProduction::SetNominalBeamCondtion() tar lab "<<_beamNucl.P4()<<std::endl;
   std::cout<<" PhotoProduction::SetNominalBeamCondtion() photon prest "<<_nuclRestPhot<<std::endl;
   std::cout<<" PhotoProduction::SetNominalBeamCondtion() tar prest "<<_nuclRestNucl<<std::endl;

   }
 
  /////////////////////////////////////////////////////////////////////////
  ////Set up the decay, adjust for given and kinematic thresholds etc
  void PhotoProduction::InitGen(){
    std::cout<<"PhotoProduction InitGen "<<std::endl;
    //pass on lorentzvectors in nucleon rest frame
    //This is the internal frame for the generator
    _reactionInfo._target=&_nuclRestNucl;
    _reactionInfo._photon =&_nuclRestPhot;
 
    auto& unproducts=Model()->UnstableProducts();
    //if(unproducts.empty()==true) return;
    
    if(unproducts.size()!=1) {
      std::cerr<<"PhotoProduction::InitGen need a W model with just a gamma + N decay product"<<std::endl;
    }

    double minMass=_massIon;
    if(unproducts.empty()==false){
      _gammaN = unproducts[0];//should only be gamma*N decaying product
      minMass=_gammaN->MinimumMassPossible();
    }
    if(auto ModelW=dynamic_cast<DecayModelW*>(Model())){
      auto thresh=ModelW->getThreshold();
      if(minMass<thresh)minMass=thresh;
    }
    
    if(Decayer()==nullptr){
      //Need to give ebeam (in ion rest), mass of ion, W threshold
      auto tempDecayer=new DontDecay();
      //tempDecayer->SetModel(Model());
      SetDecayer(tempDecayer); //give it to a sink
    }
    mutableDecayer()->PostInit(dynamic_cast<ReactionInfo*>(&_reactionInfo));

    auto decayer= dynamic_cast<DontDecay* >(mutableDecayer());
    
    if(decayer!=nullptr){

      //Do momemntum =>y, W limits last
      _Wmin=  minMass;
      auto brem=dynamic_cast<BremstrPhoton*>(_photonptr);
      if(brem!=nullptr){
	// W^2 - M^2  = 2M(Eg) 
	minMass=TMath::Sqrt(2*_massIon*brem->GetMinEnergy() + _massIon*_massIon);
	std::cout<<"PhotoProduction::InitGen()  minMass from min brem "<<minMass<<" from Emin = "<<brem->GetMinEnergy() <<std::endl;
     }
      
     }
    std::cout<<"PhotoProduction::InitGen() final minimum W "<<_Wmin<< " and energy "<<_Emin<<std::endl;
    
    if(_gammaN!=nullptr){
      _gammaN->SetMinMass(_Wmin);
       if(auto WModel=dynamic_cast<DecayModelW*>(Model())){
	 WModel->setThreshold(_Wmin);
     }
       generator().SetModelForMassPhaseSpace(_gammaN->Model());
    }

    
    std::cout<<"PhotoProduction::InitGen posit init brem"<<std::endl;
    ProductionProcess::PostInit(dynamic_cast<ReactionInfo*>(&_reactionInfo));

    //now polphotonvector has been set need to recall bremPostInit to set it
    auto brem=dynamic_cast<BremstrPhoton*>(_photonptr);
    if(brem!=nullptr){ //now we can get photoPolVector
 	brem->PostInit(dynamic_cast<ReactionInfo*>(&_reactionInfo));
    }
  
  }
  //////////////////////////////////////////////////////////////////////////
  ///Use Frixione + sigma(W) to integrate cross section over x , y and t
  double PhotoProduction::IntegrateCrossSectionFast(){
    std::cout<<"PhotoProduction::IntegrateCrossSectionFast() "<<std::endl;

    gBenchmark->Start("IntegrateCrossSectionFast");
    auto collision=MakeCollision();
    
    auto gammaNModel =dynamic_cast<DecayModelst*>(_gammaN->Model());
    auto threshold=gammaNModel->GetMeson()->PdgMass()+gammaNModel->GetBaryon()->PdgMass();
    TH1D* hWdist=new  TH1D("sdisthigh","sdisthigh",100,threshold,collision.M());
    gammaNModel->HistIntegratedXSection( *hWdist);
    
    TH1D* hEdist= new TH1D("egammadist","egammadist",100,threshold,collision.M());
    auto bremPhot =dynamic_cast<const BremstrPhoton*>(_photonptr->Decayer());
    std::cout<<"PhotoProduction::IntegrateCrossSectionFast() bremsstrahlung distribution:" <<std::endl;
    DistTF1* Edist = bremPhot->GetPhotonEdist();
    (Edist->GetTF1()).Print();
    Double_t maxE = bremPhot->GetBeamEnergy();
    std::cout<<"PhotoProduction::IntegrateCrossSectionFast() maxE = "<< maxE<<std::endl;
    for(int i=0;i<1000000;i++){
      Double_t energy = (Edist->SampleSingle())*maxE;
      LorentzVector b(0,0,energy,energy);
      LorentzVector t(0,0,0,_beamNucl.Mass());
      // LorentzVector t(0,0,0,0.938);
      Double_t W = (b+t).M();
      if(W<hEdist->GetXaxis()->GetXmin()){
        i--;
        continue;
      }
      hEdist->Fill(W);
    }

    double integrated_xsection = 0; // get sigma_ep from integral over W: f(W)*sigma_gp(W)
    for(int i=0; i<hWdist->GetNbinsX(); i++) {
      double W = hWdist->GetXaxis()->GetBinCenter(i+1);
      double WbinWidthScale = hWdist->GetBinWidth(i+1);
      double W_xsection = hWdist->GetBinContent(i+1);
      
      double W_fluxWeight =  hEdist->GetBinContent(i+1)/1000000 * hEdist->GetNbinsX();
      // double W_fluxWeight =  1;
      integrated_xsection += W_xsection * W_fluxWeight* WbinWidthScale;
    }
    gBenchmark->Stop("IntegrateCrossSectionFast");
    gBenchmark->Print("IntegrateCrossSectionFast");
 
    return integrated_xsection;

    // return 0.0;
   }

 
  LorentzVector PhotoProduction::MakeCollision(){
    //Generate collision 4-momentum
    if(_photonptr!=nullptr){
      _photonptr->GenerateComponents();
      _beamPhot.SetP4(*(_photonptr->GetInteracting4Vector()));
    }
    if(_targetptr!=nullptr){
      _targetptr->GenerateComponents();
      _beamNucl.SetP4(*(_targetptr->GetInteracting4Vector()));
    }

    //First, Eventually want to sample from beam divergence distributions
    LorentzVector collision = _beamPhot.P4() + _beamNucl.P4();
    //Boost into ion rest frame
    auto prBoost=_beamNucl.P4().BoostToCM();
    collision=boost(collision,prBoost);
    SetBoostToLab(-prBoost);
    _nuclRestNucl=LorentzVector(0,0,0,_beamNucl.Mass());
    _nuclRestPhot= boost(_beamPhot.P4(),prBoost);
    //set decay parent for gamma + p
    SetXYZT(collision.X(),collision.Y(),collision.Z(),collision.T());
    //   std::cout<<"PhotoProduction::MakeCollision() "<<collision<<std::endl;
    return collision;
  }
  //////////////////////////////////////////////////////////////////////////
  ///
  double PhotoProduction::IntegrateCrossSection(){
    std::cout<<"PhotoProduction::IntegrateCrossSection() needs impmented "<<std::endl;
 
    return 0.0;
  }
/////////////////////////////////////////////////////////////////////////
  DecayStatus  PhotoProduction::GenerateProducts(){

    auto collision=MakeCollision();
    
    //proceed through decay chain
    while(DecayingParticle::GenerateProducts()!=DecayStatus::Decayed){
      _nsamples++;
      collision=MakeCollision();
    }//DecayModelW
    
     
    //Boost all stable particles back to lab
    auto prBoost=_beamNucl.P4().BoostToCM();
    Manager::Instance().Particles().BoostStable(-prBoost);
   
    return DecayStatus::Decayed;
  }

}

