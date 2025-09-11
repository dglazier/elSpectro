#include "TwoBodyProduction.h"
#include "JpacTwoBody.h"
#include "FunctionsForGenvector.h"
#include "TwoBodyEnvelope.h"
#include "TwoBodytEnvelope.h"
#include "VectorUtils.h"
#include <TDatabasePDG.h>
#include <TBenchmark.h>
#include <Math/GSLIntegrator.h>
#include <Math/IntegrationTypes.h>
#include <Math/Functor.h>
#include <Math/Minimizer.h>
#include <Math/Factory.h>

namespace elSpectro{

  ///////////////////////////////////////////////////////
  ///constructor includes subseqent decay of Ngamma* system
  TwoBodyProduction::TwoBodyProduction(const decaying_objs& decs, const particle_objs& stables) :
    ProductionModel{ decs, stables }
  {
    _name={"TwoBodyProduction"};

    if(Products().size()!=2){
      Fatal("TwoBodyProduction","Can only have two decay particles");
    }
    
    MakeMesonFirst();
    
     //set intial values
    _p4baryon=GetBaryon()->P4();
    _p4meson=GetMeson()->P4();
      
    
  }
  /////////////////////////////////////////////////////////////////
  void TwoBodyProduction::PostInit(ReactionInfo* info){
    DecayModel::PostInit(info);
 
     _prodInfo= dynamic_cast<ReactionElectroProd*> (info); //I need Reaction info
     if(_prodInfo==nullptr){
       _isElProd=kFALSE;
       _prodInfo= dynamic_cast<ReactionPhotoProd*> (info);
     }
     
     
     //  Build_WCosTh_Envelope();
     MakeMesonFirst();

  
  }
  void TwoBodyProduction::MakeMesonFirst(){
    
    //need to find meson and baryon
    if(TDatabasePDG::Instance()->GetParticle(Products()[0]->Pdg())->ParticleClass()==TString("Baryon") && TDatabasePDG::Instance()->GetParticle(Products()[1]->Pdg())->ParticleClass()==TString("Meson")){
      //We need the meson first for TwoBody decay vectors t-distribution
      std::cout<<"TwoBodyProduction swap " <<Products()[0]->Pdg()<<" "<<Products()[0]<<" "<<Products()[1]->Pdg()<<" "<<Products()[1]<<" 223 products  "<<dynamic_cast<DecayingParticle*>(Products()[1])->Model()->Product(0)->Pdg()<<" "<<dynamic_cast<DecayingParticle*>(Products()[1])->Model()->Product(0)<<" "<<dynamic_cast<DecayingParticle*>(Products()[1])->Model()->Product(1)<<" "<<dynamic_cast<DecayingParticle*>(Products()[1])->Model()->Product(1)->Pdg()<<std::endl;
      SwapProducts(0,1);
      SetBaryonIdx(1);
      SetMesonIdx(0);
      std::cout<<"TwoBodyProduction swap " <<Products()[0]->Pdg()<<" "<<Products()[0]<<" "<<Products()[1]->Pdg()<<" "<<Products()[1]<<" products  "<<dynamic_cast<DecayingParticle*>(Products()[0])->Model()->Product(0)->Pdg()<<" "<<dynamic_cast<DecayingParticle*>(Products()[0])->Model()->Product(0)<<" "<<dynamic_cast<DecayingParticle*>(Products()[0])->Model()->Product(1)<<" "<<dynamic_cast<DecayingParticle*>(Products()[0])->Model()->Product(1)->Pdg()<<std::endl;
    }
    else  if(TDatabasePDG::Instance()->GetParticle(Products()[1]->Pdg())->ParticleClass()==TString("Baryon") && TDatabasePDG::Instance()->GetParticle(Products()[0]->Pdg())->ParticleClass()==TString("Meson")){
      SetBaryonIdx(1);
      SetMesonIdx(0);
    }
    else {
      std::cerr<<"TwoBodyProduction Need a MEson and Baryon as defined by root TDatabasePDG "<<Products()[0]->Pdg()<<" "<<TDatabasePDG::Instance()->GetParticle(Products()[0]->Pdg())->ParticleClass()<<Products()[1]->Pdg()<<" "<<TDatabasePDG::Instance()->GetParticle(Products()[1]->Pdg())->ParticleClass()<<std::endl;
      exit(0);
    }
 
  }
  
  double TwoBodyProduction::Intensity() const
  {
    _isSampling=true;
    //return 1.;
    _p4parent=Parent()->P4();
    _W = _p4parent.M();

    //take copies if all currently known particles
    _p4photon = _prodInfo->_photon;
    _p4target = _prodInfo->_target;
    //    _ebeam = _prodInfo->_ebeam;
    //    _photonPol = _prodInfo->_photonPol;
    //From my own decayer
    _p4baryon=GetBaryon()->P4();
    _p4meson=GetMeson()->P4();

    //std::cout<<"Start TwoBodyProduction::Intensity w = "<<get_W()<<" masses "<<_p4meson.M()<<" "<<_p4baryon.M()<<std::endl;

    //make sure mass is above threshold, as meson and baryon masses may have changed
    if( _W < (_p4meson.M()+_p4baryon.M()) ) return 0;
 
    //calculate thetaCM etc. for this event
    CalcKine();

    //    if(_isElProd==kTRUE) ElectroProduction();//photon polarisation etc.
    double weight = DiffXS();
    //double weight =  MatrixElementsSquared_T();
     ++_Ntries;
     // _totalXS +=  DiffXS();

    //  if(weight>get_max()){
    //    //  set_max(weight);
    //   ++_NhighWeight;
    //   std::cout<<"TwoBodyProduction::Intensity() high weight % "<<_NhighWeight/_Ntries<<" "<<weight/get_max() <<std::endl;
    // }

     auto wmaxDist = _distHighXS.GetValueFor(_W);
     double wmax = 0.;//DiffXS_at_tmin();
     wmax = wmax < wmaxDist ? wmaxDist : wmax;
     wmax = wmax < (wmaxDist=_distHighXS.GetValueInterpolated(_W)) ? wmaxDist : wmax ;
     // wmax = wmax < (wmaxDist=_distHighXS.GetValueForBinAbove(_W)) ? wmaxDist : wmax ;
     // wmax = wmax < (wmaxDist=_distHighXS.GetValueForBinBelow(_W)) ? wmaxDist : wmax ;

     
     //     wmax = wmax < wmaxDist=_distHighXS.GetValueFor(_W) ? wmaxDist :wmax ;
     //std::cout<< _distHighXS.GetValueFor(_W)<<" "<<_distHighXS.GetValueInterpolated(_W)<<" "<<_distHighXS.GetValueForBinAbove(_W)<<" "<<_distHighXS.GetValueForBinBelow(_W)<<" "<<_distHighXS.GetTH1().FindFixBin(_W)<<" out of "<<_distHighXS.GetTH1().GetNbinsX() <<std::endl;
     if(wmax==0)wmax=1.;
     wmax*=Q2PhaseSpaceCorrect();
     wmax*=MassPhaseSpaceCorrect();
     wmax*=2;//fudge factor
    //
     //std::cout<<" TwoBodyProduction::Intensity w = "<<get_W()<<" weight "<<weight<<" "<<wmax<<" new weight"<< weight/wmax<<" cos  "<<get_cosThCM()<<" t "<<get_t()<<" or "<<"kin_tFromWCosTh(get_W(),get_cosThCM()) "<<" min "<<kin_tFromWCosTh(get_W(),1)<<" Q2 "<<-_p4photon.M2() <<std::endl;
      weight/=wmax;
      //weight*=0.5; ///the reduces the efficiecny of the sampling, but also decreases the probability of current weight being grerater than the sampling weight. This may happen due to different values of t for a given cosTheta due to the mass not being PDG
      
      if(weight>1.){
	auto cmBoost=_p4parent.BoostToCM();
	auto p1cm=boost(_p4photon,cmBoost);
 	auto pdgt = p1cm.M2() + 0.77526000*0.77526000 - 2 * (p1cm.E()* kinCM_MesonE(get_W())-p1cm.P()* kinCM_MesonP(get_W())*get_cosThCM());

	auto tQ20= kine::tFromcosthW(get_cosThCM(),get_W(),0.0,_p4target.M(),_p4meson.M(),_p4baryon.M());
	auto tQ20PDG= kine::tFromcosthW(get_cosThCM(),get_W(),0.0,_p4target.M(),GetMeson()->PdgMass(),GetBaryon()->PdgMass());
	
	auto pdgEquivCosTh=kine::costhFromt(tQ20PDG, get_W(),0.0,_p4target.M(),_p4meson.M(),GetBaryon()->PdgMass());
	std::cout<<kine::tFromcosthW(pdgEquivCosTh,get_W(),0.0,_p4target.M(),GetMeson()->PdgMass(),GetBaryon()->PdgMass())<<std::endl;
	//     auto tempt = kin_tFromWCosThatQ20(_W,samplesCth[i]);
	//     // decMeson_ptr->SetP4M(1);
	//     ROOT::Math::PxPyPzM4D<double> tempp4{_p4meson.X(),_p4meson.Y(),_p4meson.Z(),1};
	//     _p4meson.SetXYZT(tempp4.X(),tempp4.Y(),tempp4.Z(),tempp4.T());
	//     std::cout<<_W<<" "<<hist.GetBinContent(ih)<<" "<<i<< " "<<samplesCth[i]<<" xs "<<samplesVal[i]<<" "<<tempt<<" with max mass "<<F(samplesCth[i])<<std::endl;
	//     // decMeson_ptr->SetP4M(tempmass);
	//     _p4meson.SetXYZT(tempmass.X(),tempmass.Y(),tempmass.Z(),tempmass.T());

	//	std::cout<<" TwoBodyProduction::Intensity w = "<<get_W()<<" weight "<<weight<<" max "<<wmax<<"xs at min "<< DiffXS_at_tmin()<< " new xs"<< weight*wmax<< " alternative "<< weight*wmax/DiffXS_at_tmin()<<" cos  "<<get_cosThCM()<<" t "<<get_t()<<" or "<<kin_tFromWCosTh(get_W(),get_cosThCM())<<" or pdg  "<< pdgt<<" or Q20 "<<tQ20<<" or Q20PDG "<<tQ20PDG<<" min "<<kin_tFromWCosTh(get_W(),1)<<" Q2 "<<-_p4photon.M2()<<" t from particles "<< (_p4meson-_p4photon).M2()<<" target "<<_p4target.M()<<" baryon "<<_p4baryon.M()<<" baryon "<<_p4meson.M()<<" surronding max "<<_distHighXS.GetValueFor(_W)<<" "<<_distHighXS.GetValueForBinAbove(_W)<<" "<<_distHighXS.GetValueForBinBelow(_W)<<" pdg mass t xs "<<pdgEquivCosTh<<" "<<_distEnvelope.GetValueFor(_W,pdgEquivCosTh+1)<<" "<<std::endl;
	//	exit(0);
      }
     //     weight/=get_max(); //normalise range 0-1
     //   return _distMassFrac->GetWeightFor(_W);
     //return 1;
    return weight;
     
  }

  

  ///////////////////////////////////////////////////////////////////////
  //creates cross section as function of W
  //and makes envelope from integration sample points.
  
  TH1D  TwoBodyProduction::CrossSectionW(TH1D hist) {
    gBenchmark->Start("CrossSectionW");
    
    std::cout<<"TwoBodyProduction::CrossSectionW "<< std::endl ;
    auto M1 = 0;//assume real photon for calculation
    _p4target= _prodInfo->_target;
    auto M2 = _p4target.M();
    auto M3 = GetMeson()->Mass(); //should be pdg value here
    auto M4 = GetBaryon()->Mass();
    //auto Wmin = M3+M4;
    auto Wmin = Parent()->MinimumMassPossible();
    auto Wmax = hist.GetXaxis()->GetXmax();
    //integrate over costh
    // std::cout<<"TwoBodyProduction::CrossSectionW "<< std::endl ;
 
    //histogram fraction of mass distribution allowed as function of W
    //used to suppress sub-threshold cross sections
   auto* decMeson_ptr = dynamic_cast<DecayingParticle*>(GetMutableMeson());
   auto* decBaryon_ptr = dynamic_cast<DecayingParticle*>(GetMutableBaryon());
   std::cout<<"TwoBodyProduction::CrossSectionW "<< GetMutableMeson()->Pdg()<<" "<<decMeson_ptr<<" "<<GetMutableBaryon()->Pdg()<<" "<<decBaryon_ptr<<std::endl ;
   //turn off below threshold
   decMeson_ptr = nullptr;
   decBaryon_ptr = nullptr;
   
   //copy Particle version for using here
   DecayingParticle decMeson;
   DecayingParticle decBaryon;
   double mesMinMass=0;
   double barMinMass=0;
   Double_t mesonFull=0;
   Double_t baryonFull=0;
   if(decMeson_ptr!=nullptr){
     decMeson= *decMeson_ptr;
     mesMinMass=decMeson_ptr->MinimumMassPossible();
     mesonFull = decMeson.MassDistribution()->Integrate1DX();
    }
   else{
     mesMinMass=GetMeson()->MinimumMassPossible();
   }
   if(decBaryon_ptr!=nullptr){
     decBaryon= *decBaryon_ptr;
     barMinMass=decBaryon_ptr->MinimumMassPossible();
     baryonFull = decBaryon.MassDistribution()->Integrate1DX();
    }
   else{
     barMinMass=GetBaryon()->MinimumMassPossible();    
   }

   Bool_t fullReached=kFALSE; //past kinematics with full mass distribution
   
   Double_t randWmin=0;
   Double_t randWmax=0;
   
   std::vector<double> samplesVal;
   std::vector<double> samplesCth;
   
   
   auto F = [this,&samplesVal,&samplesCth](double costh)
    {
      //nned to calculate t from costh and W
      set_t(kin_tFromWCosThatQ20(get_W(),costh));
      //set_t(kine::tFromcosthW(costh,get_W(),0.0,_p4target.M(),_p4meson.M(),_p4baryon.M()));
      //now can get cross section, no further kinematics caclualted
      auto result = dsigma_dcosth();
      //samplesCth.push_back(costh);
      samplesCth.push_back(get_t());
      samplesVal.push_back(result);
        
      //samplesVal.push_back(MatrixElementsSquared_T());
      //      std::cout<<"F "<<get_W()<<" "<<costh<<" t= "<<get_t()<<" "<<result<<std::endl;
      return  result;
    };
    
    ROOT::Math::GSLIntegrator ig(ROOT::Math::IntegrationOneDim::kADAPTIVE,
				 ROOT::Math::Integration::kGAUSS61);
    ROOT::Math::Functor1D wF(F);
    ig.SetFunction(wF);

    //initialise with xbin at histogram low edge
    //and y given by the costheta limits -1,1
    std::vector<double> wbins={hist.GetXaxis()->GetBinLowEdge(1)};
    std::vector<std::vector<double> > cthbins={std::vector<double>{-1,1}};
    _W = wbins.front();
    std::vector<std::vector<double> > values={std::vector<double>{F(-1),F(1)}};

 
    for(int ih=1;ih<=hist.GetNbinsX();ih++){
      //for(int ih=1;ih<=5;ih++){
      samplesCth.clear();
      samplesVal.clear();
      //std::cout<<"W " <<ih<<" "<<hist.GetXaxis()->GetBinCenter(ih)<<std::endl;
      set_W(hist.GetXaxis()->GetBinCenter(ih));
      wbins.push_back(hist.GetXaxis()->GetBinCenter(ih));
      //if(_W>1.7) exit(0);
      // if( _W < Wmin )
      // 	hist.SetBinContent(ih, 0);
      // else{
	//to account for meson mass distribution :
	//   evaluate at  mean mass allowed up to current W (allows sub-threshold)
	//   multiply by fraction of mass distribution to currentW (suppress sub-threshold )
	// i.e we assume the 2-body cross section integrates over full mass distribution. If the mass distribution is kinematically inhibited, then we reduce the 2-body cross section accordingly
	Double_t massFrac = 1.0;
	if(decMeson_ptr!=nullptr) decMeson_ptr->TakePdgMass();
	if(decBaryon_ptr!=nullptr) decBaryon_ptr->TakePdgMass();
	if(fullReached==kFALSE){
	    
	  if(decMeson_ptr!=nullptr){
	    auto mesonMean = decMeson.MeanMass(_W,barMinMass);
	    decMeson.SetP4M(mesonMean);
	    _p4meson=decMeson.P4();
	    // decMeson_ptr->SetP4M(mesonMean);
	    std::cout<<"TwoBodyProduction::CrossSectionW " <<_W<<" range "<<mesMinMass<<" "<<get_W()-barMinMass<<" barmin "<<barMinMass<<" mean "<<mesonMean<<" "<<decMeson.Mass()<<std::endl;
	  }
	  if(decBaryon_ptr!=nullptr){
	    auto baryonMean = decBaryon.MeanMass(_W,mesMinMass);
	    decBaryon.SetP4M(baryonMean);
	    _p4baryon=decBaryon.P4();
	    //decBaryon_ptr->SetP4M(baryonMean);
	  }
	    

	  //When integrating it is not clear what to use for other mass
	  //if other does not decay, it must be PDG mass
	  //if other does decay, then not well defined, currently just use minimum mass. PDG might be better
	  if(decMeson_ptr!=nullptr)
	    massFrac*= decMeson.IntegratedMass(_W,barMinMass)/mesonFull;
	  if(decBaryon_ptr!=nullptr)
	    massFrac*= decBaryon.IntegratedMass(_W,mesMinMass)/mesonFull;
	} //finshed calculating correction
	
	  //reached full mass distribution, don't need to correct anymore.
	if(massFrac==1.0) fullReached=kTRUE;
	
	F(1);F(-1);//make sure we sample costheta limits
       	hist.SetBinContent(ih, ig.Integral(-1,1)*massFrac*massFrac);
	//	hist.SetBinContent(ih, ig.Integral(-1,1));
	
	// for(auto& val:samplesVal){
	//   val+=hist.GetBinContent(ih)*0.5;
	// }
	
	cthbins.push_back(samplesCth);
	values.push_back(samplesVal);
	  
        // std::cout<<"TwoBodyProduction::CrossSectionW " <<ih<<" "<<_W <<" "<<massFrac<<" mesonFull "<<mesonFull<<" "<<hist.GetBinContent(ih)<<" nsamples "<<samplesCth.size()<<std::endl;
	// if(ih==2){
	//   auto nsamples  = samplesCth.size();
	//   for(auto i =0;i<nsamples;++i){

	//     auto tQ20= kine::tFromcosthW(samplesCth[i],_W,0.0,_p4target.M(),1,_p4baryon.M());
	//     auto tQ20PDG= kine::tFromcosthW(samplesCth[i],_W,0.0,_p4target.M(),GetMeson()->PdgMass(),_p4baryon.M());

	//     auto tempmass = _p4meson;
	//     auto tempt = kin_tFromWCosThatQ20(_W,samplesCth[i]);
	//     decMeson_ptr->SetP4M(decMeson.P4().M());
	//     ROOT::Math::PxPyPzM4D<double> tempp4{_p4meson.X(),_p4meson.Y(),_p4meson.Z(),1};
	//     _p4meson.SetXYZT(tempp4.X(),tempp4.Y(),tempp4.Z(),tempp4.T());
	//     std::cout<<_W<<" "<<samplesVal[i]<<" "<<i<< " "<<samplesCth[i]<<" xs "<<samplesVal[i]<<" t "<<tempt<<" more t "<<tQ20<<" "<<tQ20PDG<<std::endl;
	//     // decMeson_ptr->SetP4M(0.967365);
	//     // _W = 2.71125;
	//     // auto tdata= kine::tFromcosthW(0.968506,_W,0.0,_p4target.M(),0.967365,_p4baryon.M());
	//     // std::cout<<"gen event "<<_W<<" "<<tdata<<" "<<F(0.968506)<<" target "<<_p4target.M()<<" baryon "<<_p4baryon.M()<<" "<<std::endl;
	//     // // _p4meson.SetXYZT(tempmass.X(),tempmass.Y(),tempmass.Z(),tempmass.T());
	//      decMeson_ptr->SetP4M(0.77526);
	//   }
	//   exit(0);
	// }
	
	// if(ih>2)exit(0);
	//      }
      //	hist.SetBinContent(ih, 1 );
      if( TMath::IsNaN(hist.GetBinContent(ih)) )hist.SetBinContent(ih, 0 );

      // std::cout<<"TwoBodyProduction::CrossSectionW " <<ih<<" "<<_W <<" "<<hist.GetBinContent(ih)<<std::endl;

    }
    std::cout<<"TwoBodyProduction::CrossSectionW total cross-section" <<" "<<hist.Integral("w")<<" threshold "<<hist.GetXaxis()->GetBinCenter(1)<<std::endl;
    std::cout<<"DEBUG TwoBodyProduction::CrossSectionW" <<" "<<cthbins.size()<<" "<<wbins.size()<<std::endl;
    
    if(decMeson_ptr!=nullptr) decMeson_ptr->TakePdgMass();
    if(decBaryon_ptr!=nullptr) decBaryon_ptr->TakePdgMass();

    //create DistYGivenX
    _distEnvelope = DistYGivenX{wbins,cthbins,values};						      
    // _distHighXS.reset( new DistTH1(dist.GetMaxValVersusX()) );
    _distHighXS=DistTH1{_distEnvelope.GetMaxValVersusX()};	      
    Parent()->SetDecayer(CloneDecayer(TwoBodytEnvelope(_distEnvelope)));
    gBenchmark->Stop("CrossSectionW");
    gBenchmark->Print("CrossSectionW");

    return hist;
  }
 
 //////////////////////////////////////////////////////////
 void TwoBodyProduction::Use_WCosTh_Envelope(bool use){
   _createMyEnvelope=use;
 }

  /*
  //////////////////////////////////////////////////////////
  void TwoBodyProduction::Build_WCosTh_Envelope()  {
    gBenchmark->Start("envelope");//timer
    std::cout<<"TwoBodyProduction::Build_WCosTh_Envelope"<<std::endl;
    auto meson = Products()[idx::Meson()];
    auto baryon = Products()[idx::Meson()];
    auto M1 = 0;//assum real photon for calculation
    auto M2 = _p4target.M();
    auto M3 = meson->Mass(); //should be pdg value here
    //auto M3 = _meson->MaximumMassPossible();// _meson->Mass(); //should be pdg value here
    auto M4 = baryon->Mass();

    auto minW = Parent()->MinimumMassPossible();

    auto maxW = _prodInfo->_Wmax;

    //Define cosTheta binning
    //Need more resolution close to cosTheta=1
    Int_t Ncth=20; //safe at 100, but try 50 for speed
    std::vector<double > cthBins;
   
   //  Float_t convertToLog=18./Ncth;//18 as = range -10 to 10 (-2) due to the +1 in TMath::Log10. The makes the argument have a minimum value of 1 and log10(10)=0
   // for(int ib=-Ncth/2;ib<=Ncth/2;++ib){
   //    auto sign = ib==0 ? 0 : (ib/TMath::Abs(ib));
   //    cthBins.push_back(sign*TMath::Log10(TMath::Abs(float(ib)*convertToLog)+1));
   //   }
    
    cthBins.push_back(TMath::Cos(TMath::Pi()/100)); //add in small width bin at low t. Improves performance when large peak in smallest bin.
    cthBins.push_back(TMath::Cos(TMath::Pi()-TMath::Pi()/100/4));
    for(int ib=0;ib<=Ncth;ib++){
      double theta = TMath::Pi()*ib/Ncth;
      double cth = TMath::Cos(theta);
      cthBins.push_back(cth);
    }
    std::sort(cthBins.begin(),cthBins.end());
    

    //W bins want focussed on threshold
    std::vector<double > WBins;
    int NW=20;//safe at 60 but much faster with 20!
    double WRange = maxW-minW;
    double deltaW = WRange/NW;
    for(int iW=0;iW<NW+1;++iW){
      int Nsteps = NW-iW;
      if(Nsteps==0)Nsteps==1;
      if(iW==0) Nsteps = 100; //extra at threshold
      for(int iWi=0;iWi<Nsteps;++iWi){
	WBins.push_back(minW + iW*deltaW+static_cast<double>(iWi*deltaW)/Nsteps );
     }
    }
    WBins.push_back(maxW);
    //increasing order
    std::sort(WBins.begin(),WBins.end());
    //create envelope histogram
    TH2D hist("WCosTh","WCosTh",WBins.size()-1,WBins.data(),cthBins.size()-1,cthBins.data());

    //consider cross section at 3 points.
    //minimum, PDG and maximum to make sure we get relaiable maximum
    //for full meson mass range
    auto M3s=std::vector<double>{meson->Mass()}; //safe with all 3, but testing with just 1
    
    //  if(_meson->Mass()!=_meson->MinimumMassPossible()) M3s.push_back(_meson->MinimumMassPossible());
    // if(_meson->Mass()!=_meson->MaximumMassPossible()) M3s.push_back(_meson->MaximumMassPossible());

    //Loop over meson masses
    for(int Mpos=0;Mpos<M3s.size();++Mpos){ //take max of upper,mid and lower bin values
      M3=M3s[Mpos];

      //jpac 2-body matrix elements squared probably dont depend on mass
      auto pdgM3 = _meson->Mass();

      //histogram with x-axis = W ; y-axis = cos(theta)
      for(int Wpos=-1;Wpos<1;++Wpos){ //take max of upper,mid and lower bin values
	for(int ih=1;ih<=hist.GetNbinsX();ih++){
	  _W=hist.GetXaxis()->GetBinCenter(ih) + Wpos*hist.GetXaxis()->GetBinWidth(ih)/2;
	  // if(_W>105) exit(0);
  //create a correction to account for differences in phase space
	  //due to mass
	  auto PScorr=kine::PDK(_W,_meson->MinimumMassPossible(),_p4baryon->P4().M())/kine::PDK(_W,_meson->P4().M(),_baryon->P4().M());
	  if(TMath::IsNaN(PScorr)) PScorr=1;
	  
	  _meson->SetXYZT(0,0,0,M3);

	  if(ih==1&&Wpos==-1){ //move slightly above threshold
	    //threshold will not have high cross section
	    _W+=hist.GetXaxis()->GetBinWidth(ih)/10;
	  }
	  _s=_W*_W;      

	  //loop over cos theta
	  for(int ic=1;ic<=hist.GetNbinsY();ic++){
	    //function for returning cross section
	    auto evalXS = [&ic,&M1,&M2,&M3,&M4,&hist,this](const double binFact, double& val){
	      _cosThCM=hist.GetYaxis()->GetBinCenter(ic)
		+binFact*hist.GetYaxis()->GetBinWidth(ic);
	      if( _cosThCM>1) _cosThCM=1.0;
	      if( _cosThCM<-1) _cosThCM=-1.0;
	      _t = kine::tFromcosthW(_cosThCM, _W, M1, M2, M3, M4);
	      val=MatrixElementsSquared_T();
		//	      val=dsigma_dcosth();
	      if(TMath::Abs(val)==TMath::Infinity()) val=0;
	      return;
	    };
	  
	    //0 below threshold
	    if( _W < minW ){
	      hist.SetBinContent(ih,ic, 0);
	    }
	    
	    //zoom in near t=0
	    //exponential can increase very fast and dip just at end
	    //so need to look with extra resolution
	    if(ic==hist.GetNbinsY()){
	      auto cthmax=hist.GetYaxis()->GetBinCenter(ic)
		+0.5*hist.GetYaxis()->GetBinWidth(ic);
	      auto cthmin=hist.GetYaxis()->GetBinCenter(ic)
		-0.5*hist.GetYaxis()->GetBinWidth(ic);
	      //std::cout<<_W<<" cos limits "<<cthmax<<" "<<cthmin<<" "<<hist.GetYaxis()->GetBinCenter(ic)<<std::endl;
	      
	      auto tmax = kine::tFromcosthW(cthmax, _W, M1, M2, M3, M4);
	      auto tmin = kine::tFromcosthW(cthmin, _W, M1, M2, M3, M4);
	      Double_t trange= tmax-tmin;
	      Int_t nsteps = 10;
	      Double_t maxvalatt =0;
	      //step in t closest to edge looking for highest value
	      for(UInt_t istep = 0;istep<nsteps;istep++){
		_t = tmax-istep*trange/nsteps/20;
		_cosThCM=kine::costhFromt(_t, _W,M1,M2,M3,M4);
		//Double_t val=dsigma_dcosth();
		Double_t val=MatrixElementsSquared_T();

		//std::cout<<val<< " "<<_W<<" "<<" costh "<<_cosThCM<<" t "<<_t<<std::endl;
	      if(TMath::Abs(val)==TMath::Infinity()) val=0;
	      if(val>maxvalatt) maxvalatt=val;
	      }
	      //use the highest value
	      if(hist.GetBinContent(ih,ic)<maxvalatt){
		hist.SetBinContent(ih,ic, maxvalatt);
	      }
	    
	    }//if last bin
	    
	  
	  //look for the highest value across the bin
	    std::vector<double> vals(3);
	    
	    evalXS(-0.5,vals.at(0));
	    evalXS(0,vals.at(1));
	    evalXS(0.5,vals.at(2));
	    
	    std::sort(vals.begin(),vals.end(),std::greater<double>());
	    if(TMath::IsNaN(vals[0])) {
	      if(hist.GetBinContent(ih,ic)==0) hist.SetBinContent(ih,ic, 0);
	    }
	    else{
	      if(hist.GetBinContent(ih,ic)<vals[0])
		hist.SetBinContent(ih,ic, vals[0]*PScorr);//applying PS fudge
	    }

	    //    std::cout<<"TwoBodyProduction result for bin W = "<<_W<<" "<<_cosThCM<<" "<<hist.GetBinContent(ih,ic)<<std::endl;
	}
	
      }
    }//Wpos
    }//Mpos
 
  
 
    for(int ih=1;ih<=hist.GetNbinsX();ih++){
     auto maxVal= 0;
     for(int ic=1;ic<=hist.GetNbinsY();ic++){
       if(hist.GetBinContent(ih,ic)>maxVal)
	 maxVal=hist.GetBinContent(ih,ic);
     }
     for(int ic=1;ic<=hist.GetNbinsY();ic++){
       //  std::cout<<" twobody envelope "<<hist.GetXaxis()->GetBinCenter(ih)<<" "<<hist.GetYaxis()->GetBinCenter(ic)<<" "<<hist.GetBinContent(ih,ic)<<std::endl;
       //hist.SetBinContent(ih,ic,hist.GetBinContent(ih,ic)+maxVal*0.02 );
       hist.SetBinContent(ih,ic,hist.GetBinContent(ih,ic));
       // hist.SetBinContent(ih,ic,1);
     }
   }
     //adjust max after addition
   auto  maxhist = hist.GetMaximum();
   _max=maxhist*1.1;

   _histHighXS = TH1D{"high2BodyXSvW","high2BodyXSvW",hist.GetXaxis()->GetNbins(),hist.GetXaxis()->GetXbins()->GetArray()};
   
   for(int ih=1;ih<=hist.GetNbinsX();ih++){
     Double_t maxXs=0.0;
     Double_t integral =0.0;
     
     for(int ic=1;ic<=hist.GetNbinsY();ic++){
       auto currXs = hist.GetBinContent(ih,ic);
       integral+=currXs*hist.GetXaxis()->GetBinWidth(ic);
       if(currXs>maxXs)maxXs=currXs;
     }
     _histHighXS.SetBinContent(ih,maxXs);//*1.05);
     //_histHighXS.SetBinContent(ih,hist.GetXaxis()->GetBinCenter(ih));
     std::cout<<ih<<" max "<<maxXs<<" at "<<hist.GetXaxis()->GetBinCenter(ih)<<" integral "<<integral<<" "<<_histHighXS.GetBinContent(ih)<<" to "<<hist.GetXaxis()->GetBinCenter(ih)+hist.GetXaxis()->GetBinWidth(ih)/2<<std::endl;
   }
   _distHighXS.reset( new DistTH1(_histHighXS) );
   //helps to appply some smoothing to histogram!
   // hist.Smooth(1);
   //apply envelope if requested
   // if(_createMyEnvelope)Parent()->SetDecayer(new TwoBodyEnvelope(DistTH2Slice{hist}));
   gBenchmark->Stop("envelope");//timer
   gBenchmark->Print("envelope");//timer
 
   std::cout<<"TwoBodyProduction::Build_WCosTh_Envelope max xsec = "<<_max<<std::endl;

  }
  */
  void TwoBodyProduction::Print() const{
    DecayModel::Print();
    std::cout<<"TwoBodyProduction Average cross section sampled "<<_totalXS/_Ntries<<std::endl;
 
  }

}
