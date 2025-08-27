#include "ElectronScattering.h"
#include "FormationQ2W.h"
#include "TwoBodyProduction.h"
#include "FunctionsForGenvector.h"
#include "FunctionsForKinematics.h"
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
  // ElectronScattering::ElectronScattering(double ep,double ionp, DecayModel* model, int ionpdg):
  //   _pElectron{ep},
  //   _pIon{ionp},
  //   _angleElectron{TMath::Pi()},
  //   _angleIon{0},
  //   _pdgIon{ionpdg},
  //   _beamElec{11},
  //   _beamNucl{ionpdg},
  //   ProductionProcess{0,nullptr,model}
  // {
      
  //     SetBeamCondtion();
  // }
  // /////////////////////////////////////////////////////////////////////
  // ElectronScattering::ElectronScattering(double ep,double ionp,
  // 		     double eangle,double ionangle,  DecayModel* model, int ionpdg):
  //   _pElectron{ep},
  //   _pIon{ionp},
  //   _angleElectron{eangle},
  //   _angleIon{ionangle},
  //   _pdgIon{ionpdg},
  //   _beamElec{11},
  //   _beamNucl{ionpdg},
  //   ProductionProcess{0,nullptr,model}
  // {
  //     SetBeamCondtion();
  // }
  /////////////////////////////////////////////////////////////////////
  ElectronScattering::ElectronScattering(const CollidingParticle& electron,const CollidingParticle& target,decaymodel_ptr model):
    ProductionProcess{electron,target,model},
    _beamElec{11},
    _beamNucl{target.GetInteractingPdg()}
  {
    //For convenience keep our own
    //pointers to electron and target
    _electronptr =Incident1(); //order given to ProductionProcess
    _targetptr = Incident2();
    
    SetNominalBeamCondtion();
  }

  /////////////////////////////////////////////////////////////////////
  void ElectronScattering::SetNominalBeamCondtion(){
    
    
    _beamElec.SetP4(_electronptr->GetNominal4Vector());
    
    _beamNucl.SetP4(_targetptr->GetNominal4Vector());
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
    std::cout<<"ElectronScattering vertex "<<_electronptr->VertexPosition()<<" "<<_electronptr<< " "<<_targetptr<<std::endl;
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
    std::cout<<"ElectronScattering::SetBeamCondtion() ptrs "<<&_beamElec<<" "<<&_beamNucl<<std::endl;

   }
  /////////////////////////////////////////////////////////////////////////
  void ElectronScattering::InitGen(){
    std::cout<<"Electron Scattering InitGen 1"<<std::endl;
    //pass on lorentzvectors in nucleon rest frame
    //This is the internal frame for the generator
    _reactionInfo._target=_nuclRestNucl;
    _reactionInfo._ebeam =_nuclRestElec;

    //DecayModelQ2W should be initialised with single product gamma*N
    //auto& unproducts=Model()->Products();
    std::cout<<"ElectronScattering::InitGen() Model "<<Model()->GetName()<<std::endl;
									   // auto* formation = dynamic_cast<FormationQ2W*>(Model());
    // if(unproducts.size()!=1) {
    //   std::cerr<<"ElectronScattering::InitGen need a Q2W model with just a gamma*N decay product"<<std::endl;
									   //}
    std::cout<<"Electron Scattering InitGen 2 "<<Q2WModel()<<" "<<std::endl;
    double minMass=_massIon;
    //if(unproducts.empty()==false){
    auto& gStarN = Q2WModel()->GetGammaN();
    std::cout<<"Electron Scattering min mass "<<gStarN.MinimumMassPossible()<<std::endl;
    minMass=gStarN.MinimumMassPossible();
    //}
  
    //    if(auto Q2WModel=dynamic_cast<FormationQ2W*>(Model())){
    auto thresh=Q2WModel()->getThreshold();
    if(minMass<thresh)minMass=thresh;
    // }
    std::cout<<"Electron Scattering InitGen 2"<<std::endl;
  
    //default scatteredelectron_xy, now have all parameters
    // if(Decayer()==nullptr){
    //Need to give ebeam (in ion rest), mass of ion, W threshold
    // auto tempDecayer=new ScatteredElectron_xy(_nuclRestElec.P(), _massIon, minMass);
    //tempDecayer->SetModel(Model());
    //SetDecayer0(tempDecayer); //give it to a sink
    // SetDecayer0(new ScatteredElectron_xy(_nuclRestElec.P(), _massIon, minMass));
    SetDecayer0(std::make_shared<ScatteredElectron_xy>(_nuclRestElec.P(), _massIon, minMass));
    
    // mutableDecayer()->PostInit(dynamic_cast<ReactionInfo*>(&_reactionInfo));

    auto decayer= dynamic_cast<ScatteredElectron_xy* >(mutableDecayer());
    decayer->SetModel(Model());
    std::cout<<"Electron Scattering InitGen 3"<<std::endl;

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
      std::cout<<"Electron Scattering InitGen 4"<<std::endl;

      if(_ePmax!=0||_Ymin!=0) { //convert to y limit
	//Find lowest allowed y
	double y =0;
	if(_ePmax!=0) y = (_nuclRestElec.E()-escat::E_el(_ePmax))/_nuclRestElec.E();
	if(_Ymin!=0){
	  if(y<_Ymin){
	    y=_Ymin;
	    // W^2 - M^2 + Q2 = 2M(Eg) = 2M*ebeam*y
	    auto tempminMass=TMath::Sqrt(2*_massIon*_nuclRestElec.E()*y + _massIon*_massIon);
	    if(tempminMass<minMass){ //ymin is below threshold, set to 0 to use threshold
	      _Ymin=0;
	      y=0;
	    }
	    else{
	      minMass=tempminMass;
	    }
	    
	    _ePmax=_nuclRestElec.E() - _nuclRestElec.E()*y;
	  }

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
    //if(_gStarN!=nullptr){
    Q2WModel()->GetGammaN().SetMinMass(_Wmin);
    //  if(auto Q2WModel=dynamic_cast<FormationQ2W*>(Model())){
    Q2WModel()->setThreshold(_Wmin);
      //   }
      //      generator().SetModelForMassPhaseSpace(_gStarN->Model());
      //}

    //try here so chance to redefine minimum masses etc
    ProductionProcess::PostInit(dynamic_cast<ReactionInfo*>(&_reactionInfo));
 
 
  }
  //////////////////////////////////////////////////////////////////////////
  ///Use Frixione + sigma(W) to integrate cross section over x , y and t
  double ElectronScattering::IntegrateCrossSectionFast(TwoBodyProduction* model2body){
    gBenchmark->Start("IntegrateCrossSectionFast");
    auto collision=MakeCollision();

    //auto Q2WModel =dynamic_cast<FormationQ2W*>(Model());
    Q2WModel()->ZeroPhoton();//need Q2=0 for P1CM
    
    // //   auto gStarModel =dynamic_cast<DecayModelst*>(_gStarN->Model());
    //auto model2body =dynamic_cast<ProductionModel*>(Q2WModel()->GetGammaN().Model());
    auto threshold=model2body->GetMeson()->PdgMass()+model2body->GetBaryon()->PdgMass();
    
    TH1D hWdist("sdisthigh","sdisthigh",100,threshold,collision.M());
    std::cout<<"ElectronScattering::IntegrateCrossSectionFast() "<<threshold<<" "<<collision.M()<<std::endl;
    hWdist=model2body->CrossSectionW( hWdist);
    
    double integrated_xsection = 0; // get sigma_ep from integral over W: f(W)*sigma_gp(W)

    
    for(int i=0; i<hWdist.GetNbinsX(); i++) {
      double W = hWdist.GetXaxis()->GetBinCenter(i+1);
      double WbinWidthScale = hWdist.GetBinWidth(i+1);
      double W_xsection = hWdist.GetBinContent(i+1);
      
       //change to ion mass 8.3.2023
      double y = (W*W-_massIon)/_nuclRestElec.E()/2/_massIon;
      double W_fluxWeight = escat::Frixione(_nuclRestElec.E(),y,_massIon) * W /_nuclRestElec.E() /_massIon;
      // std::cout<<"ElectronScattering::IntegrateCrossSectionFast() W = "<< W<<" photoXS = "<<W_xsection<<" photoFlux = "<<W_fluxWeight<<" xs "<<W_xsection * W_fluxWeight* WbinWidthScale<<std::endl;
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
      _beamElec.SetP4(_electronptr->GetInteracting4Vector());
    }
    if(_targetptr!=nullptr){
      _targetptr->GenerateComponents();
      _beamNucl.SetP4(_targetptr->GetInteracting4Vector());
    }
    //First, Eventually want to sample from beam divergence distributions
    LorentzVector collision = _beamElec.P4() + _beamNucl.P4();
    //Boost into ion rest frame
    auto prBoost=_beamNucl.P4().BoostToCM();
    collision=boost(collision,prBoost);
    SetBoostToLab(-prBoost);
    _nuclRestNucl=LorentzVector(0,0,0,_beamNucl.Mass());
    _nuclRestElec= boost(_beamElec.P4(),prBoost);
    //set decay parent for e -> e'g*
    SetXYZT(collision.X(),collision.Y(),collision.Z(),collision.T());
    // std::cout<<" ElectronScattering::MakeCollision()  "<<(_targetptr->GetInteracting4Vector())<<" "<<_electronptr<<_beamElec.P4()<<_beamNucl.P4()<<_nuclRestNucl<<_nuclRestElec<<collision<<std::endl;
    return collision;
  }
 //////////////////////////////////////////////////////////////////////////
  ///Use RooFit integrator to integrate cross section over x , y and t
  double ElectronScattering::IntegrateCrossSection(TwoBodyProduction* model2body){
    
    auto collision=MakeCollision();

    auto photonFlux= dynamic_cast<ScatteredElectron_xy* >(mutableDecayer());
    //photonFlux->Dist().SetWThresholdVal(model2body->GetMeson()->PdgMass()+model2body->GetBaryon()->PdgMass());
    //photonFlux->Dist().SetWThreshold(model2body->GetMeson()->PdgMass()+model2body->GetBaryon()->PdgMass());

    
    auto xvar = RooRealVar(Form("xIntegral%lf_%lf",photonFlux->Dist().GetMinLnX(),photonFlux->Dist().GetMaxLnX()),"xIntegral",(photonFlux->Dist().GetMinLnX()),(photonFlux->Dist().GetMinLnX()),(photonFlux->Dist().GetMaxLnX()),"");
    auto yvar = RooRealVar(Form("yIntegral%lf_%lf",photonFlux->Dist().GetMinLnY(),photonFlux->Dist().GetMaxLnY()),"yIntegral",(photonFlux->Dist().GetMinLnY()),(photonFlux->Dist().GetMinLnY()),(photonFlux->Dist().GetMaxLnY()),"");

    xvar.Print("v");yvar.Print("v");

    
    
    //DEBUG
    // auto model2body =dynamic_cast<DecayModelst*>(_gStarN->Model());
    //auto model2body =dynamic_cast<TwoBodyProduction*>(Q2WModel()->GetGammaN().Model());
    auto fastIntegral=IntegrateCrossSectionFast(model2body);
    std::cout<<"       check fast cross section "<<fastIntegral<<std::endl;

  
    
    auto Eel=_nuclRestElec.E();
    
    double_t threshW= model2body->GetMeson()->PdgMass()+model2body->GetBaryon()->PdgMass();
    Double_t maxVal=0;
    auto fXYcosth = [this,&photonFlux,&model2body,&Eel,&threshW,&maxVal](const double *x)
      {
	if(x[0]==0) return 0.; //x
	if(x[1]==0) return 0.; //y
	auto val = photonFlux->Dist().Eval(x);
	if(TMath::IsNaN(val)) return 0.;
	if(val==0) return 0.;
	//calculate scatered electron at x and y 
	photonFlux->GenerateGivenXandY(P4(),Model()->Products(),TMath::Exp(x[0]),TMath::Exp(x[1]));
	
 	
	//calculate virtual photon
	Q2WModel()->Intensity();
	//	std::cout<<"fXYcosth "<<Q2WModel()->getW()<<" "<<Q2WModel()->getQ2()<<" "<<val<<" xs "<<x[0]<<" "<<x[1]<<" "<<x[2]<<std::endl;
	if(Q2WModel()->getW()<threshW) return 0.0;
	//get value of dsigma(s)/dcosth cross section at x,y,costh
	//Double_t dsigma_dcosth=model2body->dsigma_costh(x[2]);
	Double_t dsigma_dcosth=model2body->dsigma_dcosth(model2body->get_W_FromParent(),x[2]);
	//std::cout<<"fXYcosth "<<Q2WModel()->getW()<<" "<<Q2WModel()->getQ2()<<" "<<val<<" "<<dsigma_dcosth<<std::endl;
	val*=dsigma_dcosth;
	//std::cout<<"fXYcosth "<<Q2WModel()->getW()<<" "<<model2body->get_W_FromParent()<<" "<<Q2WModel()->getQ2()<<" "<<val<<" "<<dsigma_dcosth<<" "<<x[2]<<" "<<model2body->dsigma_dcosth(model2body->get_W_FromParent(),1)<<std::endl;
	//additional (not real photo) Q2dependence of cross section
	if(TMath::IsNaN(val)) return 0.;
	if(val<0) return 0.;
	val*=Q2WModel()->Q2H1Rho();

	//	std::cout<< "check W "<<Q2WModel()->getW()<<" Q2 "<<Q2WModel()->getQ2()<<" pdg1 "<<Model()->Products()[0]->Pdg()<<" pdg2 "<<Model()->Products()[1]->Pdg()<<" "<<model2body->get_W_FromParent()<<" "<<threshW<<" "<<x[2]<<" "<<val<<std::endl;

	return val;
      };
    //std::vector<double> xvals = {-15.6763, -0.0841589, 0.49065};
  //auto res =  fXYcosth(xvals.data());
    //std::cout<<"res = "<<res<<std::endl;
    // exit(0);
   
    auto wrapPdf=ROOT::Math::Functor( fXYcosth , 3);

    //Append integral number to name to prevent RooFit cahce if not wanted
    //Note call SetCacheIntegrals() to use cahced values
    double RFintegral=0.0;

    //Split the integration up for more accurate low t integration
    //Numbers can be quite different if try and itnegrate
    //full range only 
    std::vector<double> cosThMin={-1,0.98,0.995};
    std::vector<double> cosThMax={0.98,0.995,1};
    //std::vector<double> cosThMin={0.9};
    //std::vector<double> cosThMax={1};
    gBenchmark->Start("RooFitIntegral");
    for(ushort iint=0;iint<1;iint++){
      TString pdfname(Form("ElScatterIntegral%d",NintegralsElectronScattering));
      auto cthvar = RooRealVar("CosThIntegral","CosThIntegral",0.8,cosThMin[iint],cosThMax[iint],"");
      auto pdf = RooFunctorPdfBinding(pdfname, "ElScatterIntegral", wrapPdf, RooArgList(xvar,yvar,cthvar));
      if(_cacheIntegrals==0) NintegralsElectronScattering++;//work around RooFit agressive caching!
      // pdf->Print();
      
      auto roovars= RooArgSet(xvar,yvar,cthvar);
      
 
      RFintegral+=pdf.getNorm(roovars);
    }
    
    gBenchmark->Stop("RooFitIntegral");
    gBenchmark->Print("RooFitIntegral");
     
    std::cout<<" ElectronScattering::IntegrateCrossSection()  "<<RFintegral<<" nb "<<std::endl<<" giving a photon flux weighted average photoproduction cross section of "<<RFintegral/photonFlux->Dist().Integral()<<" nb"<<std::endl;
    std::cout<<" W range "<<_Wmin<<" - "<< collision.M() <<" =  "<< ( collision.M()- _Wmin)<<std::endl;
    //xvar.Print();
    //yvar.Print();
    //cthvar.Print();
    photonFlux->Dist().SetWThresholdVal(Q2WModel()->getThreshold());
   
    return RFintegral;
  }
/////////////////////////////////////////////////////////////////////////
  DecayStatus  ElectronScattering::GenerateProducts(const ProductionProcess* production){

    auto collision=MakeCollision();

    //Assign the Production vertex
    InitVertex();

    //choose formation channel
    //this will be based on integrated xsect
    //for each 2-body final state
    Product().ChooseDecay();
    
    //proceed through decay chain
    while(DecayingParticle::GenerateProducts(this)!=DecayStatus::Decayed){
      _nsamples++;
      collision=MakeCollision();
      // std::cout<<"ElectronScattering::GenerateProducts() next event "<<_nsamples<<std::endl;
    }//DecayModelQ2W
    
     
    //Boost all stable particles back to lab
    auto prBoost=-_beamNucl.P4().BoostToCM();
    //    Manager::Instance().Particles().BoostStable(-prBoost);
    //Manager::Instance().Particles().BoostToFrame(-prBoost,collision);

    particle_ptrs final_state;
    EventParticles(final_state);//collect all final state particles
   // std::cout<<"ElectronScattering final particles "<<final_state.size()<<std::endl;
   // for(auto& p:final_state){std::cout<<" "<<p->Pdg();}
   // std::cout<<std::endl;
   kine::BoostParticles(prBoost,final_state); //boost back to lab
   SetFinalParticles(final_state); //assign to process
    
    return DecayStatus::Decayed;
  }

}

