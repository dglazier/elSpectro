#include "DecayModelQ2W.h"
#include "DecayingParticle.h"
#include "FunctionsForElectronScattering.h"
#include <TDatabasePDG.h>


namespace elSpectro{
  //////////////////////////////////////////////////////
  ////Constructor for e- scattering kinematics only
  DecayModelQ2W::DecayModelQ2W( double thresh) :
    _threshold{thresh},
    DecayModel{{},{-2211,11}}
  {
    _name={"DecayModelQ2W"};

    Init();
  }
  ///////////////////////////////////////////////////////
  ///constructor includes subseqent decay of Ngamma* system
  /* DecayModelQ2W::DecayModelQ2W( double thresh, DecayModel* gNmodel) :
    _threshold{thresh},
    DecayModel{{ new DecayingParticle{-2211,primary} },{11}}
  {
    _name={"DecayModelQ2W_with_primary_decay"};
    Init();
    }*/
  ///////////////////////////////////////////////////////
  ///constructor includes subseqent decay of Ngamma* system
  DecayModelQ2W::DecayModelQ2W( double thresh,
				DecayModel* gNmodel,DecayVectors* gNdecayer) :
    _threshold{thresh},
    DecayModel{{ new DecayingParticle{-2211,gNmodel,gNdecayer} },{11}}
  {
    _name={"DecayModelQ2W_with_primary_decay_and_decayer"};
    Init();
  }
  ////////////////////////////////////////////////////////
  ///complete constructor
  void DecayModelQ2W::Init(){
    if(Products()[0]->Pdg()==11){
      _gstarNuc=Products()[1]; //-2211
      _electron=Products()[0]; //11 
    }
    else{
      _gstarNuc=Products()[0]; //-2211
      _electron=Products()[1]; //11
    }
 
    _electron->Print();
    _gstarNuc->Print();
    
    auto gNprods=dynamic_cast<DecayingParticle*>(_gstarNuc)->Model()->Products();
    if( TString("Baryon")==TDatabasePDG::Instance()
	->GetParticle(gNprods[0]->Pdg())->ParticleClass() ){
      //Make sure meson is product 0 and baryon product 1
      dynamic_cast<DecayingParticle*>(_gstarNuc)->Model()->SwapProducts(0,1);
    }
 
    if(_threshold<MinimumMassPossible() )_threshold=MinimumMassPossible(); 
  }

  ////////////////////////////////////////////////////////
  void DecayModelQ2W::PostInit(ReactionInfo* info){
      _prodInfo = dynamic_cast<ReactionElectroProd*> (info);
      //std::cout<<"DecayModelQ2W::PostInit "<<_electron->P4ptr()<<" "<<_gstarNuc->P4ptr()<<std::endl;
      if( _prodInfo==nullptr) std::cerr<<"DecayModelQ2W PostInit not an ElectronScattering reaction"<<std::endl;
      _prodInfo->_scattered=_electron->P4ptr();
      _prodInfo->_photoN=_gstarNuc->P4ptr();
      _prodInfo->_photon=&_gamma;
      _prodInfo->_photonPol=&_photonPol;
      
      auto gNprods=dynamic_cast<DecayingParticle*>(_gstarNuc)->Model()->Products();
      // std::cout<<"DecayModelQ2W::PostInit "<<_prodInfo<<" "<<gNprods[0]->Pdg()<<" "<<gNprods[1]->Pdg()<<std::endl;
      
      _prodInfo->_baryon=gNprods[1]->P4ptr();
      _prodInfo->_meson=gNprods[0]->P4ptr();
      
      
      std::cout<<"DecayModelQ2W::PostInit "<<std::endl;
      DecayModel::PostInit(_prodInfo);
      
 
      FindExcitationSpectra();
     

  }
  
  ////////////////////////////////////////////////////////
  double  DecayModelQ2W::Intensity() const{
    // std::cout<<"DecayModelQ2W::Intensity "<<MinimumMassPossible()<<" "<<ParentVector().M()<<" "<<getW()<<std::endl;
    /*if(CheckThreshold()==false){
      return 0.;
      }*/

    double W = getW();
    //if(W <GetGammaN()->MinimumMassPossible()) return 0;
    if(TMath::IsNaN(W)) return 0.0;
    if(W  < _threshold ) return 0.;
 
  
    //calculate virtual photon
    const auto& p4beam=*(_prodInfo->_ebeam);
    const auto& p4tar=*(_prodInfo->_target);
    const auto& p4scat=_electron->P4();

    _gamma = p4beam-p4scat;//can now use getQ2
    
    //calculate photon polarisation
    auto epsilon = escat::virtualPhotonPolarisation(p4beam,p4tar,p4scat);
    auto delta = 2*escat::M2_el()/getQ2()*(1-epsilon);
    
    _photonPol.SetEpsilon(epsilon);
    _photonPol.SetDelta(delta);

    //Get envelope weight from integrated cross section
    double weight=_Wrealphoto_Dist->GetWeightFor( W  );

    if(getQ2() > 2*p4tar.M()*_gamma.E()){
      std::cout<<"Q2 above max how ? "<<getQ2()<<" 2Mmu "<<2*p4tar.M()*_gamma.E() <<" W "<<W<<std::endl;
      exit(0);
    }
    
 
    //  std::cout<<" Q2 DEPENDENECE "<<PhaseSpaceFactorToQ2eq0(W,p4tar.M() )<<"      "<<getQ2()<<" 2Mmu "<<2*p4tar.M()*_gamma.E() <<" W "<<W<<"             PDKs     "<< kine::PDK(W, -getQ2(),p4tar.M())<<" "<< kine::PDK(W, getQ2(),p4tar.M())<<" "<< kine::PDK(W, 0 ,p4tar.M())<<" "<< p1cm.P()<<std::endl;


    _prodInfo->_sWeight=weight; //might be used in s and t

    //now add Q2 depedence to get weighted here
 
    //Q2 dependence of phase space needed to effectively multiply st max value
    //Here the Q2 dependence would multiply weight value and max ,thus cancelling
    //i.e do not do weight*=PhaseSpaceFactorToQ2eq0(W,p4tar.M() );
    //_prodInfo->_sWeight*=PhaseSpaceFactorToQ2eq0(W,p4tar.M() );
    
    
    if(weight>1){
    auto cmBoost=_gstarNuc->P4().BoostToCM();
    auto p1cm=boost(_gamma,cmBoost);
      std::cout<<" Q2 DEPENDENECE "<<PhaseSpaceFactorToQ2eq0(W,p4tar.M() )<<"      "<<getQ2()<<" 2Mmu "<<2*p4tar.M()*_gamma.E() <<" W "<<W<<"             PDKs     "<< kine::PDK(W, -getQ2(),p4tar.M())<<" "<< kine::PDK(W, getQ2(),p4tar.M())<<" "<< kine::PDK(W, 0 ,p4tar.M())<<" "<< p1cm.P()<<std::endl;
    }
    
    //Q2 dependence of cross section
    weight*=Q2H1Rho();
  
    return weight;
  
  }

  void DecayModelQ2W::FindExcitationSpectra(){

    //Make excitation spectra envelope which should be greater than the cross section
    //for all values of meson mass, to do this take maximum from values at
    //meson threshold and pdg mass values
    //note the point is phase space and therefore xsection changes with mass
    double maxW = ( *(_prodInfo->_target) + *(_prodInfo->_ebeam) ).M();

    std::cout<<"DecayModelQ2W::PostInit generating total cross section, may take some time... "<<std::endl;
    auto gNprods=dynamic_cast<DecayingParticle*>(_gstarNuc)->Model()->Products();
      
    //auto baryon = gNprods[1];
    auto meson=gNprods[0];
 
    DecayModelst* mesonBaryon = nullptr;
    TH1D histlow("Wdistlow","Wdistlow",100,_threshold,maxW);
    TH1D histpeak("Wdisthigh","Wdisthigh",100,_threshold,maxW);
    double minMesonMass=-1;
    if( ( mesonBaryon=dynamic_cast<DecayModelst*>(GetGammaN()->Model())) != nullptr){
      //check for low mass meson limits
      if(dynamic_cast<DecayingParticle*>(meson)){ //meson
	dynamic_cast<DecayingParticle*>(meson)->TakeMinimumMass();//to get threshold behaviour
	minMesonMass=meson->Mass();
	//mesonBaryon->HistIntegratedXSection(histlow);
	mesonBaryon->HistMaxXSection(histlow);
	//back to PDg mass if exists
	if(meson->PdgMass()>minMesonMass)
	  dynamic_cast<DecayingParticle*>(meson)->TakePdgMass();

      }
      //now for PDG mass
      //only needs to be done if meson does not decay
      //or pdg mass is different from minMesonMass (possible if !=0)
      if(meson->PdgMass()!=minMesonMass)
 	mesonBaryon->HistMaxXSection(histpeak);
      //mesonBaryon->HistIntegratedXSection(histpeak);
   
     auto hist = HistFromLargestBinContents(histlow,histpeak);
     
     _Wrealphoto_Dist.reset( new DistTH1(hist) );
    }
    else{
      std::cerr<<"DecayModelQ2W::FindExcitationSpectra()Need a DecayModelst"<<std::endl;
      exit(0);
    }
  }
    
  
   TH1D HistFromLargestBinContents(const TH1D& h1,const TH1D& h2){
      auto hist= TH1D{h1};
      auto maxVal= h1.GetMaximum();
      for(int ibin=1;ibin<=hist.GetNbinsX();ibin++){
	auto val = h1.GetBinContent(ibin)>h2.GetBinContent(ibin) ? h1.GetBinContent(ibin):h2.GetBinContent(ibin);
	hist.SetBinContent(ibin,val + 0.05*maxVal);
	std::cout<<"HistFromLargestBins "<<hist.GetBinCenter(ibin)<<" "<<hist.GetBinContent(ibin)<<" "<< maxVal<<std::endl;
	
      }
      return hist;
   }
  
}
