#include "ElectronScattering.h"
#include "FunctionsForGenvector.h"
#include "Manager.h"
#include "Interface.h" //for generator
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
  int ElectronScattering::NintegralsElectronScattering=0;

  /////////////////////////////////////////////////////////////////////
  ElectronScattering::ElectronScattering(double ep,double ionp, DecayModel* model, int ionpdg):
    _pElectron{ep},
    _pIon{ionp},
    _angleElectron{TMath::Pi()},
    _angleIon{0},
    _pdgIon{ionpdg},
    _beamElec{11},
    _beamNucl{ionpdg},
    ProductionProcess{0,nullptr,model}
  {
      
      SetBeamCondtion();
  }
  /////////////////////////////////////////////////////////////////////
  ElectronScattering::ElectronScattering(double ep,double ionp,
		     double eangle,double ionangle,  DecayModel* model, int ionpdg):
    _pElectron{ep},
    _pIon{ionp},
    _angleElectron{eangle},
    _angleIon{ionangle},
    _pdgIon{ionpdg},
    _beamElec{11},
    _beamNucl{ionpdg},
    ProductionProcess{0,nullptr,model}
  {
      SetBeamCondtion();
  }
  /////////////////////////////////////////////////////////////////////
  ElectronScattering::ElectronScattering(CollidingParticle *electron,CollidingParticle* target,  DecayModel* model):
    ProductionProcess{electron,target,model},
    _beamElec{11},
    _beamNucl{target->GetInteractingPdg()}
  {
    //For convenience keep our own
    //pointers to electron and target
    _electronptr =electron;
    _targetptr = target;
    
      SetNominalBeamCondtion();
  }

  /////////////////////////////////////////////////////////////////////
  void ElectronScattering::SetNominalBeamCondtion(){
    
    
    _beamElec.SetP4(*(_electronptr->GetInteracting4Vector()));
    
    _beamNucl.SetP4(*(_targetptr->GetInteracting4Vector()));
    _massIon=_beamNucl.PdgMass();

    /*  //not sure why this was there
    //For decaying
    SetXYZT(_beamElec.P4().X(),_beamElec.P4().Y(),
	    _beamElec.P4().Z(),_beamElec.P4().T());
    */

    //For nominal beam conditions set nucleon rest frame vectors
    //in case info is needed in PostInit stage
    //Boost into ion rest frame
    auto prBoost=_beamNucl.P4().BoostToCM();
    _nuclRestNucl=LorentzVector(0,0,0,_beamNucl.Mass());
    _nuclRestElec= boost(_beamElec.P4(),prBoost);

    //set inital lab particles
    //this can be written to output file
    std::cout<<"ElectronScattering vertex "<<_electronptr->VertexPosition()<<std::endl;
    AddInitialParticlePtr(_electronptr);
    AddInitialParticlePtr(_targetptr);
  //set inital lab particles
    //AddInitialParticlePtr(&_beamElec);
    //AddInitialParticlePtr(&_beamNucl);


    std::cout<<" ElectronScattering::SetNominalBeamCondtion() e- lab "<<_beamElec.P4()<<std::endl;
   std::cout<<" ElectronScattering::SetNominalBeamCondtion() tar lab "<<_beamNucl.P4()<<std::endl;
   std::cout<<" ElectronScattering::SetNominalBeamCondtion() e- prest "<<_nuclRestElec<<std::endl;
   std::cout<<" ElectronScattering::SetNominalBeamCondtion() tar prest "<<_nuclRestNucl<<std::endl;

   }
  /////////////////////////////////////////////////////////////////////
  void ElectronScattering::SetBeamCondtion(){
    
    _massIon=TDatabasePDG::Instance()->GetParticle(_pdgIon)->Mass();

    
    _beamElec.SetXYZT(0,0,_pElectron,
		      escat::E_el(_pElectron));
    auto p4=_beamElec.P4();
    genvector::LorentzRotateY(p4,_angleElectron);
    _beamElec.SetP4(p4);
    std::cout<<"ElectronScattering::SetBeamCondtion() Electron "<< _beamElec.P4()<<std::endl;
    
    
    _beamNucl.SetXYZT(0,0,_pIon,
		      TMath::Sqrt(_pIon*_pIon + _massIon*_massIon));
    p4=_beamNucl.P4();
    genvector::LorentzRotateY(p4,_angleIon);
    _beamNucl.SetP4(p4);
    std::cout<<"ElectronScattering::SetBeamCondtion() Nucl "<< _beamNucl.P4()<<std::endl;
 
    //For decaying
    SetXYZT(_beamElec.P4().X(),_beamElec.P4().Y(),
	    _beamElec.P4().Z(),_beamElec.P4().T());


    //For nominal beam conditions set nucleon rest frame vectors
    //in case info is needed in PostInit stage
    //Boost into ion rest frame
    auto prBoost=_beamNucl.P4().BoostToCM();
    _nuclRestNucl=LorentzVector(0,0,0,_beamNucl.Mass());
    _nuclRestElec= boost(_beamElec.P4(),prBoost);

    //set inital lab particles
    AddInitialParticlePtr(&_beamElec);
    AddInitialParticlePtr(&_beamNucl);


   }
  /////////////////////////////////////////////////////////////////////////
  void ElectronScattering::InitGen(){
    std::cout<<"Electron Scattering InitGen "<<std::endl;
    //pass on lorentzvectors in nucleon rest frame
    //This is the internal frame for the generator
    _reactionInfo._target=&_nuclRestNucl;
    _reactionInfo._ebeam =&_nuclRestElec;
    
     
 
    auto& unproducts=Model()->UnstableProducts();
    //if(unproducts.empty()==true) return;
    
    if(unproducts.size()!=1) {
      std::cerr<<"ElectronScattering::InitGen need a Q2W model with just a gamma*N decay product"<<std::endl;
    }

    double minMass=_massIon;
    if(unproducts.empty()==false){
      _gStarN = unproducts[0];//should only be gamma*N decaying product
      std::cout<<"Electron Scattering min mass "<<_gStarN->MinimumMassPossible()<<std::endl;
      minMass=_gStarN->MinimumMassPossible();
    }
    if(auto Q2WModel=dynamic_cast<DecayModelQ2W*>(Model())){
      auto thresh=Q2WModel->getThreshold();
      if(minMass<thresh)minMass=thresh;
    }
    
    //default scatteredelectron_xy, now have all parameters
    if(Decayer()==nullptr){
      //Need to give ebeam (in ion rest), mass of ion, W threshold
      auto tempDecayer=new ScatteredElectron_xy(_nuclRestElec.P(), _massIon, minMass);
      tempDecayer->SetModel(Model());
      SetDecayer(tempDecayer); //give it to a sink
      
    }
    mutableDecayer()->PostInit(dynamic_cast<ReactionInfo*>(&_reactionInfo));

    auto decayer= dynamic_cast<ScatteredElectron_xy* >(mutableDecayer());
    
    if(decayer!=nullptr){
      //Set any thresholds and ranges
      if(_Q2min!=0)  decayer->Dist().SetQ2min(_Q2min);
      if(_Q2max!=0)  decayer->Dist().SetQ2max(_Q2max);
      if(_Xmin!=0)  decayer->Dist().SetXmin(_Xmin);
      if(_Xmax!=0)  decayer->Dist().SetXmax(_Xmax);
      if(_Ymin!=0)  decayer->Dist().SetYmin(_Ymin);
      if(_Ymax!=0)  decayer->Dist().SetYmax(_Ymax);
      if(_eThmin!=0)  decayer->Dist().SetThmin(_eThmin);
      if(_eThmax!=0)  decayer->Dist().SetThmax(_eThmax);

      //Do momemntum =>y, W limits last
      _Wmin=  minMass;
 
      if(_ePmax!=0||_Ymin!=0) { //convert to y limit
	//Find lowest allowed y
	double y =0;
	if(_ePmax!=0) y = (_nuclRestElec.E()-escat::E_el(_ePmax))/_nuclRestElec.E();
	if(_Ymin!=0){
	  if(y<_Ymin){
	    y=_Ymin;
	    _ePmax=_nuclRestElec.E() - _nuclRestElec.E()*y;
	  }
	  // W^2 - M^2 + Q2 = 2M(Eg) = 2M*ebeam*y
	  minMass=TMath::Sqrt(2*_massIon*_nuclRestElec.E()*y + _massIon*_massIon);
	}
      	
	std::cout<<" settting Ymin "<< y <<" "<<_Ymin<<" "<<_ePmax<<std::endl;
	decayer->Dist().SetYmin(y);

	//we can limit the W threshold if we have a given Q2max value
	//This can greatly speed up sampling when the max value is higher
	//below the Q2 allowed minimum W.
	double  Q2PQ2max=0;
	double  ThPQ2max=0;
	if(_Q2max!=0) Q2PQ2max=(_Q2max);
	if(_eThmax!=0) ThPQ2max=escat::Q2_cosThy(_nuclRestElec.E(),TMath::Cos(_eThmax),y);
	double useQ2max=0;

	if(Q2PQ2max ==0 && ThPQ2max) useQ2max = ThPQ2max ;//need lowest allowed max
	else  useQ2max = Q2PQ2max;
	
	if(useQ2max!=0){
	  auto W2min= 2*_massIon * (_nuclRestElec.E() - escat::E_el(_ePmax) )
	    + _massIon*_massIon - useQ2max;
	  if( W2min>minMass*minMass){
	    _Wmin=TMath::Sqrt(W2min);
	    minMass=_Wmin;
	  }
	}
	
      }
       if(_ePmin!=0) { //convert to y limit
	if( (_nuclRestElec.E() < escat::E_el(_ePmin)) ){
	  std::cerr<<"ElectronScattering::InitGen() Error, requested Minimum electron momentum (in proton rest frame) higher than beam energy "<< escat::E_el(_ePmin)<<" > "<<_nuclRestElec.E()<<std::endl;
	  exit(0);
	}
	decayer->Dist().SetYmax((_nuclRestElec.E()-escat::E_el(_ePmin))/_nuclRestElec.E());
      }
      //Finally set reaction threshold and this calcualte ylimits
      std::cout<<" setting lowest W as "<<minMass<<std::endl;
      decayer->Dist().SetWThreshold(minMass);
      _Wmin=minMass;
      if(dynamic_cast<DistVirtPhotFlux_xy*>(&(decayer->Dist())))_Wmin=dynamic_cast<DistVirtPhotFlux_xy*>((&decayer->Dist()))->GetWMin();//can be effected by angle limits etc
    }
    
    std::cout<<"ElectronScattering::InitGen() final minimum W "<<_Wmin<<std::endl;
    if(_gStarN!=nullptr){
      _gStarN->SetMinMass(_Wmin);
       if(auto Q2WModel=dynamic_cast<DecayModelQ2W*>(Model())){
	 Q2WModel->setThreshold(_Wmin);
     }
      generator().SetModelForMassPhaseSpace(_gStarN->Model());
    }

    
    ProductionProcess::PostInit(dynamic_cast<ReactionInfo*>(&_reactionInfo));
    //now scattered electron should be set
    //_scattered= _reactionInfo._scattered;

  
  }
  //////////////////////////////////////////////////////////////////////////
  ///Use Frixione + sigma(W) to integrate cross section over x , y and t
  double ElectronScattering::IntegrateCrossSectionFast(){
    gBenchmark->Start("IntegrateCrossSectionFast");
    auto collision=MakeCollision();

    auto Q2WModel =dynamic_cast<DecayModelQ2W*>(Model());
    Q2WModel->ZeroPhoton();//need Q2=0 for P1CM
    
    auto gStarModel =dynamic_cast<DecayModelst*>(_gStarN->Model());
    auto threshold=gStarModel->GetMeson()->PdgMass()+gStarModel->GetBaryon()->PdgMass();
    TH1D* hWdist=new  TH1D("sdisthigh","sdisthigh",100,threshold,collision.M());
    gStarModel->HistIntegratedXSection( *hWdist);
    
    double integrated_xsection = 0; // get sigma_ep from integral over W: f(W)*sigma_gp(W)

    
    for(int i=0; i<hWdist->GetNbinsX(); i++) {
      double W = hWdist->GetXaxis()->GetBinCenter(i+1);
      double WbinWidthScale = hWdist->GetBinWidth(i+1);
      double W_xsection = hWdist->GetBinContent(i+1);
      
      double y = (W*W-escat::M2_pr())/_nuclRestElec.E()/2/escat::M_pr();
      double W_fluxWeight = escat::Frixione(_nuclRestElec.E(),y) * W /_nuclRestElec.E() /escat::M_pr();
      integrated_xsection += W_xsection * W_fluxWeight* WbinWidthScale;
    }
    gBenchmark->Stop("IntegrateCrossSectionFast");
    gBenchmark->Print("IntegrateCrossSectionFast");
 
    return integrated_xsection;
  }

 
  LorentzVector ElectronScattering::MakeCollision(){
    //Generate collision 4-momentum
    if(_electronptr!=nullptr){
      _electronptr->GenerateComponents();
      _beamElec.SetP4(*(_electronptr->GetInteracting4Vector()));
    }
    if(_targetptr!=nullptr){
      _targetptr->GenerateComponents();
      _beamNucl.SetP4(*(_targetptr->GetInteracting4Vector()));
    }

    //First, Eventually want to sample from beam divergence distributions
    LorentzVector collision = _beamElec.P4() + _beamNucl.P4();
    //Boost into ion rest frame
    auto prBoost=_beamNucl.P4().BoostToCM();
    collision=boost(collision,prBoost);
    _nuclRestNucl=LorentzVector(0,0,0,_beamNucl.Mass());
    _nuclRestElec= boost(_beamElec.P4(),prBoost);
    //set decay parent for e -> e'g*
    SetXYZT(collision.X(),collision.Y(),collision.Z(),collision.T());
    //  std::cout<<"ElectronScattering::MakeCollision() "<<collision<<std::endl;
    return collision;
  }
  //////////////////////////////////////////////////////////////////////////
  ///Use RooFit integrator to integrate cross section over x , y and t
  double ElectronScattering::IntegrateCrossSection(){
    
    auto collision=MakeCollision();
 

    auto photonFlux= dynamic_cast<ScatteredElectron_xy* >(mutableDecayer());

    // auto xvar = RooRealVar(Form("xIntegral%lf_%lf",photonFlux->Dist().GetMaxLnX(),photonFlux->Dist().GetMaxLnX()),"xIntegral",TMath::Exp(photonFlux->Dist().GetMinLnX()),TMath::Exp(photonFlux->Dist().GetMinLnX()),TMath::Exp(photonFlux->Dist().GetMaxLnX()),"");
    //auto yvar = RooRealVar(Form("yIntegral%lf_%lf",photonFlux->Dist().GetMaxLnY(),photonFlux->Dist().GetMaxLnY()),"yIntegral",TMath::Exp(photonFlux->Dist().GetMinLnY()),TMath::Exp(photonFlux->Dist().GetMinLnY()),TMath::Exp(photonFlux->Dist().GetMaxLnY()),"");

    auto xvar = RooRealVar(Form("xIntegral%lf_%lf",photonFlux->Dist().GetMaxLnX(),photonFlux->Dist().GetMaxLnX()),"xIntegral",(photonFlux->Dist().GetMinLnX()),(photonFlux->Dist().GetMinLnX()),(photonFlux->Dist().GetMaxLnX()),"");
    auto yvar = RooRealVar(Form("yIntegral%lf_%lf",photonFlux->Dist().GetMaxLnY(),photonFlux->Dist().GetMaxLnY()),"yIntegral",(photonFlux->Dist().GetMinLnY()),(photonFlux->Dist().GetMinLnY()),(photonFlux->Dist().GetMaxLnY()),"");


    
    // auto cthvar = RooRealVar("CosThIntegral","CosThIntegral",-0.,-0.99999,0.99999,"");
    //auto cthvar = RooRealVar("CosThIntegral","CosThIntegral",0.99,-0.99999,0.99999,"");
    auto cthvar = RooRealVar("CosThIntegral","CosThIntegral",0.99,-1,1,"");
    
    auto gStarModel =dynamic_cast<DecayModelst*>(_gStarN->Model());
    auto Q2WModel =dynamic_cast<DecayModelQ2W*>(Model());
    
    photonFlux->Dist().SetWThresholdVal(gStarModel->GetMeson()->PdgMass()+gStarModel->GetBaryon()->PdgMass());

    auto Eel=_nuclRestElec.E();

    double_t threshW= gStarModel->GetMeson()->PdgMass()+gStarModel->GetBaryon()->PdgMass();
    Double_t maxVal=0;
    auto fXYcosth = [this,&photonFlux,&gStarModel,&Q2WModel,&Eel,&threshW,&maxVal](const double *x)
      {
	if(x[0]==0) return 0.; //x
	if(x[1]==0) return 0.; //y
	auto val = photonFlux->Dist().Eval(x);
	if(TMath::IsNaN(val)) return 0.;
	if(val==0) return 0.;
	//calculate scatered electron at x and y 
	photonFlux->GenerateGivenXandY(P4(),Model()->Products(),TMath::Exp(x[0]),TMath::Exp(x[1]));
	//calculate virtual photon
	Q2WModel->Intensity();
	//get value of dsigma(s)/dcosth cross section at x,y,costh
	Double_t dsigma_costh=gStarModel->dsigma_costh(x[2]);
	val*=dsigma_costh;
	//additional (not real photo) Q2dependence of cross section
	if(TMath::IsNaN(val)) return 0.;
	if(val<0) return 0.;
	val*=Q2WModel->Q2H1Rho();
	return val;
      };
  
    auto wrapPdf=ROOT::Math::Functor( fXYcosth , 3);

    //Append integral number to name to prevent RooFit cahce if not wanted
    //Note call SetCacheIntegrals() to use cahced values
    TString pdfname(Form("ElScatterIntegral%d",NintegralsElectronScattering));
    auto pdf = RooFunctorPdfBinding(pdfname, "ElScatterIntegral", wrapPdf, RooArgList(xvar,yvar,cthvar));
    if(_cacheIntegrals==0) NintegralsElectronScattering++;//work around RooFit agressive caching!
    // pdf->Print();
    
    auto roovars= RooArgSet(xvar,yvar,cthvar);
     
 
    gBenchmark->Start("RooFitIntegral");

    auto RFintegral=pdf.getNorm(roovars);
  
    gBenchmark->Stop("RooFitIntegral");
    gBenchmark->Print("RooFitIntegral");
     
    std::cout<<" ElectronScattering::IntegrateCrossSection()  "<<RFintegral<<" nb "<<std::endl<<" giving a photon flux weighted average photoproduction cross section of "<<RFintegral/photonFlux->Dist().Integral()<<" nb"<<std::endl;
    std::cout<<" W range "<<_Wmin<<" - "<< collision.M() <<" =  "<< ( collision.M()- _Wmin)<<std::endl;
    //xvar.Print();
    //yvar.Print();
    //cthvar.Print();
    
    photonFlux->Dist().SetWThresholdVal(Q2WModel->getThreshold());
   
    return RFintegral;
  }
/////////////////////////////////////////////////////////////////////////
  DecayStatus  ElectronScattering::GenerateProducts(){

    auto collision=MakeCollision();
    
    //proceed through decay chain
    while(DecayingParticle::GenerateProducts()!=DecayStatus::Decayed){
      _nsamples++;
      collision=MakeCollision();
      //std::cout<<"ElectronScattering::GenerateProducts() "<<_nsamples<<std::endl;
    }//DecayModelQ2W
    
     
    //Boost all stable particles back to lab
    auto prBoost=_beamNucl.P4().BoostToCM();
    Manager::Instance().Particles().BoostStable(-prBoost);
    //Manager::Instance().Particles().BoostToFrame(-prBoost,collision);
   
    return DecayStatus::Decayed;
  }

}

