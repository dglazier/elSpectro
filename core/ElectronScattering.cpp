#include "ElectronScattering.h"
#include "FunctionsForGenvector.h"
#include "Manager.h"
#include "Interface.h" //for generator
#include "ScatteredElectron_xy.h"
#include <TDatabasePDG.h>
#include <TH1F.h>
#include <TFile.h>


namespace elSpectro{

  // ElectronScattering::ElectronScattering(double ep,double ionp,
  // 					 double eangle,double ionangle,
  // 					 DecayModel* model, int ionpdg):
  //   _pElectron{ep},
  //   _pIon{ionp},
  //   _angleElectron{eangle},
  //   _angleIon{ionangle},
  //   _pdgIon{ionpdg},
  //   ProductionProcess{model}
  // {
  //   SetBeamCondtion();
  // }
  // //////////////////////////////////////////////////////////////////
  // ElectronScattering::ElectronScattering(double ep,double ionp,
  // 					 DecayModel* model, int ionpdg):
  //   _pElectron{ep},
  //   _pIon{ionp},
  //   _angleElectron{TMath::Pi()},
  //   _angleIon{0},
  //   _pdgIon{ionpdg},
  //   ProductionProcess{model}
  // {
  //   SetBeamCondtion();
  // }
  /////////////////////////////////////////////////////////////////////
  ElectronScattering::ElectronScattering(double ep,double ionp, DecayVectors* decayer,  DecayModel* model, int ionpdg):
    _pElectron{ep},
    _pIon{ionp},
    _angleElectron{TMath::Pi()},
    _angleIon{0},
    _pdgIon{ionpdg},
    _beamElec{11},
    _beamNucl{ionpdg},
    ProductionProcess{0,decayer,model}
  {
      
      SetBeamCondtion();
  }
  /////////////////////////////////////////////////////////////////////
  ElectronScattering::ElectronScattering(double ep,double ionp,
		     double eangle,double ionangle,
		     DecayVectors* decayer,  DecayModel* model, int ionpdg):
    _pElectron{ep},
    _pIon{ionp},
    _angleElectron{eangle},
    _angleIon{ionangle},
    _pdgIon{ionpdg},
    _beamElec{11},
    _beamNucl{ionpdg},
    ProductionProcess{0,decayer,model}
  {
  
      SetBeamCondtion();
  }
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
      std::cerr<<"ElectronScattering::PostInit need a Q2W model with just a gamma*N decay product"<<std::endl;
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
    std::cout<<"E scatter go t a threshold "<<minMass<<std::endl;

    
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
      std::cout<<"ELECTRON "<<_Q2min<<" "<<_Q2max<<std::endl;
      if(_Q2min!=0)  decayer->Dist().SetQ2min(_Q2min);
      if(_Q2max!=0)  decayer->Dist().SetQ2max(_Q2max);
      if(_Xmin!=0)  decayer->Dist().SetXmin(_Xmin);
      if(_Xmax!=0)  decayer->Dist().SetXmax(_Xmax);
      if(_eThmin!=0)  decayer->Dist().SetThmin(_eThmin);
      if(_eThmax!=0)  decayer->Dist().SetThmax(_eThmax);

      //Do momemntum =>y, W limits last
      auto Wmin=  minMass;
      /* if(_ePmax!=0) {
	
	auto minMass2=minMass*minMass;
	//W^2 - M^2 + Q2 = 2M(Eg) = 2M(ebeam-escat) = 2M*ebeam*y
	if( (_nuclRestElec.E() < escat::E_el(_ePmax)) ){
	  std::cerr<<"ElectronScattering::InitGen() Error, requested Maximum electron momentum (in proton rest frame) higher than beam energy "<< escat::E_el(_ePmax)<<" > "<<_nuclRestElec.E()<<std::endl;
	  exit(0);
	}
	auto W2min= 2*_massIon * (_nuclRestElec.E() - escat::E_el(_ePmax) )
	  + _massIon*_massIon + 0;//take Q2=0 for W2 threshold, any Q2 limit will be handled in DistVirtPhotFlux
	//take largest from requested or threshold
	if( W2min>minMass2){
	  Wmin=TMath::Sqrt(W2min);
	  minMass=Wmin;
	}

      }
      */
      if(_ePmax!=0) { //convert to y limit
	auto y = (_nuclRestElec.E()-escat::E_el(_ePmax))/_nuclRestElec.E();
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
	    Wmin=TMath::Sqrt(W2min);
	    minMass=Wmin;
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
      Wmin=minMass;
    }
    
  
    if(_gStarN!=nullptr){
      _gStarN->SetMinMass(minMass);
      std::cout<<"CHECK "<<dynamic_cast<DecayModelQ2W*>(Model())<<std::endl;
      if(auto Q2WModel=dynamic_cast<DecayModelQ2W*>(Model())){
	Q2WModel->setThreshold(minMass);
      }
      generator().SetModelForMassPhaseSpace(_gStarN->Model());
    }

    
    DecayingParticle::PostInit(dynamic_cast<ReactionInfo*>(&_reactionInfo));
    //now scattered electron should be set
    //_scattered= _reactionInfo._scattered;

    //integration

    /*
    long _NIntegrateW=100;
    long _NIntegratet=10000;
    long Npass=0;
    double integral=0;

    
    static_cast<DistVirtPhotFlux_xy*>(&decayer->Dist())->ForIntegrate(true);
    // _gStarN->mutableDecayer()->ForIntegrate(true);
    for(long iW=0;iW<_NIntegrateW;++iW){
      LorentzVector collision = _beamElec.P4() + _beamNucl.P4();

      //Boost into ion rest frame
      auto prBoost=_beamNucl.P4().BoostToCM();
      collision=boost(collision,prBoost);
      _nuclRestNucl=LorentzVector(0,0,0,_beamNucl.Mass());
      _nuclRestElec= boost(_beamElec.P4(),prBoost);
      //set decay parent for e -> e'g*
      SetXYZT(collision.X(),collision.Y(),collision.Z(),collision.T());
 
      Decay(); //e'
      auto pass=Model()->Intensity();
      if(pass==0) continue;
      
      double intet=0;
      Npass++;
      //TH1D histt("t","t",100,-30,0);
     for(long it=0;it<_NIntegratet;++it){
	_gStarN->Decay(); //random t/cosTh
	_gStarN->Model()->Intensity();
	//histt.Fill(static_cast<DecayModelst*>(_gStarN->Model())->get_t());
	
	intet+=dsigma();
      }
      //integral+=intet;
      TH1D histCross("XSection","XSection",1,_gStarN->Mass()-0.001,_gStarN->Mass()+0.001);
      //double Wmid=gRandom->Uniform(5.5,6.5);
      //TH1D histCross("XSection","XSection",1,Wmid-0.001,Wmid+0.001);
      static_cast<DecayModelst*>(_gStarN->Model())->HistIntegratedXSection(histCross);
      std::cout<<"Integral for t "<<_gStarN->Mass()<<" "<<intet/_NIntegratet<<" integration "<<histCross.GetBinContent(1)<<" at "<<_reactionInfo._photon->M2()<<"phase Q2 "<<histCross.GetBinContent(1)*static_cast<DecayModelQ2W*>(Model())->PhaseSpaceFactorToQ2eq0(_gStarN->Mass(),_beamNucl.Mass())<<std::endl <<std::endl <<std::endl ;
      
      integral+=intet/_NIntegratet;
      // TFile * crossFile=new TFile("XSection.root","recreate");
      //histt.Write();
      //delete crossFile;
      //exit(0); 
    }


    //samint=myFns.Integral(0,100);sum=0;for(int i=0;i<1E6;i++){auto val=myFns.GetRandom();sum+=myFn2.Eval(val)/myFns.Eval(val)*samint;}
    
    static_cast<DistVirtPhotFlux_xy*>(&decayer->Dist())->ForIntegrate(false);
    double Wmax = (_nuclRestNucl + _nuclRestElec).M();
    double Wmin=minMass;
    std::cout<<"Integral "<< Wmin<<" - "<<Wmax<<" mean "<<integral/Npass<<"  int  "<<integral/_NIntegrateW*(Wmax-Wmin) << " "<<integral/_NIntegrateW/_NIntegratet*(Wmax-_nuclRestNucl.M()) <<std::endl <<std::endl <<std::endl;

    /*TH1D histCross("XSection","XSection",10,Wmin,Wmax);
    std::cout<< "hist integral for "<<_gStarN->Model()->GetName()<<std::endl;
    static_cast<DecayModelst*>(_gStarN->Model())->HistIntegratedXSection(histCross);
    TFile * crossFile=new TFile("XSection.root","recreate");
    histCross.Write();
    delete crossFile;
   
    exit(0);
    */
  }
/////////////////////////////////////////////////////////////////////////
DecayStatus  ElectronScattering::GenerateProducts(){
    //First, Eventually want to sample from beam divergence distributions
    LorentzVector collision = _beamElec.P4() + _beamNucl.P4();

    //Boost into ion rest frame
    auto prBoost=_beamNucl.P4().BoostToCM();
    collision=boost(collision,prBoost);
    _nuclRestNucl=LorentzVector(0,0,0,_beamNucl.Mass());
    _nuclRestElec= boost(_beamElec.P4(),prBoost);
    
    //set decay parent for e -> e'g*
    SetXYZT(collision.X(),collision.Y(),collision.Z(),collision.T());

    //proceed through decay chain
    //first get masses for all products
    // DetermineAllMasses(); //this will define W threshold...
    //make sure masses are below maximum kinematically possible
    // if(_gStarN->MinimumMassPossible() > W2Max() ){
      // DetermineAllMasses();
    //}
      
    //find a W candidate
    while(DecayingParticle::GenerateProducts()!=DecayStatus::Decayed){
      //std::cout<<"not decayed "<<_nsamples<<std::endl;
      _nsamples++;
    }//DecayModelQ2W
    
    
    //Boost all stable particles back to lab
    Manager::Instance().Particles().BoostStable(-prBoost);

    return DecayStatus::Decayed;
  }

}
