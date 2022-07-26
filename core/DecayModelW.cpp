#include "DecayModelW.h"
#include "DecayingParticle.h"
#include "FunctionsForElectronScattering.h"
#include "DistConst.h"
#include <TDatabasePDG.h>
#include "TFile.h"


namespace elSpectro{
  //////////////////////////////////////////////////////
  ////Constructor for e- scattering kinematics only
  DecayModelW::DecayModelW( double thresh) :
    _threshold{thresh},
    DecayModel{{},{-2211}}
  {
    _name={"DecayModelW"};

    Init();
  }
  ///////////////////////////////////////////////////////
  ///constructor includes subseqent decay of Ngamma* system
  DecayModelW::DecayModelW( double thresh,
				DecayModel* gNmodel,DecayVectors* gNdecayer) :
    _threshold{thresh},
    DecayModel{{ new DecayingParticle{-2211,gNmodel,gNdecayer} },{}}
  {
    _name={"DecayModelW_with_primary_decay_and_decayer"};
    Init();
  }
  ////////////////////////////////////////////////////////
  ///complete constructor
  void DecayModelW::Init(){
    _gstarNuc=Products()[0]; //-2211
   
    auto gNprods=dynamic_cast<DecayingParticle*>(_gstarNuc)->Model()->Products();
    if( TString("Baryon")==TDatabasePDG::Instance()
	->GetParticle(gNprods[0]->Pdg())->ParticleClass() ){
      //Make sure meson is product 0 and baryon product 1
      dynamic_cast<DecayingParticle*>(_gstarNuc)->Model()->SwapProducts(0,1);
    }
 
    if(_threshold<MinimumMassPossible() )_threshold=MinimumMassPossible(); 
  }

  ////////////////////////////////////////////////////////
  void DecayModelW::PostInit(ReactionInfo* info){
      _prodInfo = dynamic_cast<ReactionPhotoProd*> (info);
      _photon =  _prodInfo->_photon;
      std::cout<<"DecayModelW::PostInit "<<" "<<_gstarNuc->P4ptr()<<std::endl;
      if( _prodInfo==nullptr) std::cerr<<"DecayModelW PostInit not a Photoproduction reaction"<<std::endl;
      _prodInfo->_photoN=_gstarNuc->P4ptr();
      _prodInfo->_photonPol=&_photonPol;
      if(_prodInfo->_Wmax==0){
	//For backward compatability, should be done with
	//colliding particles now and set in ProductionProcess::PostInit().
	_prodInfo->_Wmax=( *(_prodInfo->_target) + *(_prodInfo->_ebeam) ).M();
      }
      auto gNprods=dynamic_cast<DecayingParticle*>(_gstarNuc)->Model()->Products();
      _prodInfo->_baryon=gNprods[1]->P4ptr();
      _prodInfo->_meson=gNprods[0]->P4ptr();
      
      
      std::cout<<"DecayModelW::PostInit with W threshold = "<<getThreshold()<<std::endl;
      DecayModel::PostInit(_prodInfo);
      
      if( dynamic_cast<DecayModelst*>(GetGammaN()->Model()) == nullptr){
	_Wrealphoto_Dist.reset( new DistConst(1) );
      }
      else{	
	FindExcitationSpectra();
      }
  }
  
  ////////////////////////////////////////////////////////
  double  DecayModelW::Intensity() const{
    // std::cout<<"DecayModelW::Intensity "<<MinimumMassPossible()<<" "<<ParentVector().M()<<" "<<getW()<<" > "<<_threshold<<" "<<GetGammaN()->P4().E() <<std::endl;
    /// std::cout<<"DecayModelW::Intensity "<<ParentVector()<<" gN "<<GetGammaN()->P4()<<std::endl;
    /*if(CheckThreshold()==false){
      return 0.;
      }*/

    double W = getW();
    if(TMath::IsNaN(W)) return 0.0;
    if(W  < _threshold ) return 0.;
 
  
    //Get envelope weight from integrated cross section
    double weight=_Wrealphoto_Dist->GetWeightFor( W  );
         

    _prodInfo->_sWeight=weight; //might be used in s and t
   
  
     return weight;
    
  }
  void DecayModelW::FindExcitationSpectra(){

    //Make excitation spectra envelope which should be greater than the
    //max cross section for that W value
    //for all values of meson mass, to do this take maximum from values at
    //meson threshold and pdg mass values
    //note the point is phase space and therefore xsection changes with mass
    //double maxW = ( *(_prodInfo->_target) + *(_prodInfo->_ebeam) ).M();
    double maxW = _prodInfo->_Wmax;
   
    std::cout<<"DecayModelW::PostInit generating max cross section @W, may take some time... "<<std::endl;
    auto gNprods=dynamic_cast<DecayingParticle*>(_gstarNuc)->Model()->Products();
      
    //auto baryon = gNprods[1];
    auto meson=gNprods[0];
 
    DecayModelst* mesonBaryon = nullptr;
    TH1D histlow("Wdistlow","Wdistlow",400,_threshold,maxW);
    TH1D histpeak("Wdisthigh","Wdisthigh",400,_threshold,maxW);
    double minMesonMass=-1;
    if( ( mesonBaryon=dynamic_cast<DecayModelst*>(GetGammaN()->Model())) != nullptr){
      //check for low mass meson limits
      if(dynamic_cast<DecayingParticle*>(meson)){ //meson
	dynamic_cast<DecayingParticle*>(meson)->TakeMinimumMass();//to get threshold behaviour
	minMesonMass=meson->Mass();
 	mesonBaryon->HistMaxXSection(histlow);
	//back to PDg mass if exists
	if(meson->PdgMass()>minMesonMass)
	  dynamic_cast<DecayingParticle*>(meson)->TakePdgMass();

      }
      //now for PDG mass
      //only needs to be done if meson does not decay
      //or pdg mass is different from minMesonMass (possible if !=0)
      if(meson->PdgMass()!=minMesonMass){
	mesonBaryon->HistMaxXSection(histpeak);
      }
      auto hist = HistFromLargestBinContents(histpeak,histlow);
      std::cout<<"DecayModelW::FindExcitationSpectra()  result   "<<hist.GetMaximum()<<" "<<hist.GetBinCenter(hist.GetMaximumBin())<<" "<<hist.GetNbinsX()<<std::endl;
      hist.SetName("Wdist");

        
     _Wrealphoto_Dist.reset( new DistTH1(hist) );
    }
    else{
      std::cerr<<"DecayModelW::FindExcitationSpectra()Need a DecayModelst"<<std::endl;
      exit(0);
    }
  }
    
  
  TH1D HistFromLargestBinContents(const TH1D& h1,const TH1D& h2){
      auto hist= TH1D{h1};
      auto maxVal= h1.GetMaximum();
      double max_so_far=0.;
      for(int ibin=1;ibin<=hist.GetNbinsX();ibin++){
	auto val = h1.GetBinContent(ibin)>h2.GetBinContent(ibin) ? h1.GetBinContent(ibin):h2.GetBinContent(ibin);
	//hist.SetBinContent(ibin,val + 0.05*maxVal);

	//depending on differential cross section and particle mass
	//phase space factors can increase cross section at high
	//particle mass and W
	//We want envelope to contain this so once we get to the max
	//just stay there. Flux is low at high W so no big effect
	//on efficiency from this
	
	if(val<max_so_far){
	  hist.SetBinContent(ibin,max_so_far );
	}
	else{
	  max_so_far = val + 0.05*maxVal;
	  hist.SetBinContent(ibin,val + 0.05*maxVal);
	}
      }
      return hist;
   }
  
}
