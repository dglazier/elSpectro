#include "TwoBodyProduction.h"
#include "FunctionsForGenvector.h"
#include <TDatabasePDG.h>
#include <Math/GSLIntegrator.h>
#include <Math/IntegrationTypes.h>
#include <Math/Functor.h>
#include <Math/Minimizer.h>
#include <Math/Factory.h>

namespace elSpectro{

  ///////////////////////////////////////////////////////
  ///constructor includes subseqent decay of Ngamma* system
  TwoBodyProduction::TwoBodyProduction(particle_ptrs parts, const std::vector<int> pdgs) :
    DecayModel{ parts, pdgs }
  {
    _name={"TwoBodyProduction"};

    if(Products().size()!=2){
      Fatal("TwoBodyProduction","Can only have two decay particles");
    }
    //need to find meson and baryon
    if(TDatabasePDG::Instance()->GetParticle(Products()[0]->Pdg())->ParticleClass()==TString("Baryon") ){
      _baryon=Products()[0];
      _meson=Products()[1]; 
    }
    else {
      _baryon=Products()[1];
      _meson=Products()[0];
    }
    
  }
  /////////////////////////////////////////////////////////////////
  void TwoBodyProduction::PostInit(ReactionInfo* info){
     std::cout<<"TwoBodyProduction::PostInit "<<" "<<GetName()<<std::endl;
     DecayModel::PostInit(info);
 
     _prodInfo= dynamic_cast<ReactionElectroProd*> (info); //I need Reaction info
     if(_prodInfo==nullptr){
       _isElProd=kFALSE;
       _prodInfo= dynamic_cast<ReactionPhotoProd*> (info);
     }
     
     _photon = _prodInfo->_photon;
     _target = _prodInfo->_target;
     _ebeam = _prodInfo->_ebeam;
     _photonPol = _prodInfo->_photonPol;
     


     Build_WCosTh_Envelope();
         
  
  }
  
  double TwoBodyProduction::Intensity() const
  {
    _isSampling=true;
    
    _W = Parent()->P4().M();
     //make sure mass is above threshold, as meson and baryon masses may have changed
    if( _W < (_meson->P4().M()+_baryon->P4().M()) ) return 0;
 
    //calculate thetaCM etc. for this event
    CalcKine();

    //    if(_isElProd==kTRUE) ElectroProduction();//photon polarisation etc.
    double weight = DiffXS();
     ++_Ntries;

     if(weight>get_max()){
      //std::cout<<"TwoBodyProduction::Intensity() dsig="<<weight<<" max "<<_max<<" W= "<<_W<<" t= "<<_t<<"cosTh = "<<_cosThCM<<std::endl;
      //  set_max(weight);
      ++_NhighWeight;
      std::cout<<"TwoBodyProduction::Intensity() high weight % "<<_NhighWeight/_Ntries<<" "<<weight/get_max() <<std::endl;
    }

     //apply suppression for regions with low mass phase space
     //This will be close to/sub threshold
     if(_distMassFrac.get()) weight*=_distMassFrac->GetWeightFor(_W);
     
     weight/=get_max(); //normalise range 0-1
     
     
     return weight;
     
  }

  
  //////////////////////////////////////////////////////////////////////////
   double TwoBodyProduction::FindMaxOfIntensity(){
    
    auto M1 = 0;//assum real photon for max calculation
    auto M2 = _target->M();
    auto M3 = _meson->Mass(); //should be pdg value here
    auto M4 = _baryon->Mass();
    auto Wmin = Parent()->MinimumMassPossible();
   
    auto Wmax = _prodInfo->_Wmax;
    if(Wmax == 0){
      std::cerr<<"TwoBodyProduction::FindMaxOfIntensity(), Wmax=0, should be set in ReactionInfo"<<std::endl;
      exit(0);
    }
    double Wrange=Wmax-Wmin;
 
    auto Fmax = [&M1,&M2,&M3,&M4,&Wmin,&Wmax,this](const double *x){
		  //x ={W,cosTheta}
		  _s = x[0]*x[0];
		  _W=x[0];
		  _cosThCM=x[1];
		  if( _W < Wmin ) return 0.;
		  if( _W < M3+M4 ) return 0.;
		  if( _W > Wmax ) return 0.;

		  //calculate t from nominal values
		  auto currt=kine::tFromcosthW(_cosThCM,_W,M1,M2,M3,M4);
		  auto myt0=kine::t0(_W,M1,M2,M3,M4);
		  auto mytmax=kine::tmax(_W,M1,M2,M3,M4);
		  if(currt>myt0) return 0.;
		  if(currt<mytmax) return 0.;
		  if( TMath::IsNaN(x[1]) ) return 0.;

		  _t=currt;
		  double val = DiffXS();
		  if( TMath::IsNaN(val) ) return 0.;
		  return -(val); //using a minimiser!
		};
    

    //First try with nominal (probably PDG meson and baryon mass)
    ROOT::Math::Minimizer* minimum =
      ROOT::Math::Factory::CreateMinimizer("Minuit2", "");
    //      ROOT::Math::Factory::CreateMinimizer("Genetic", "");
    
    
    // set tolerance , etc...
    minimum->SetMaxFunctionCalls(1000000); // for Minuit/Minuit2
    minimum->SetMaxIterations(1000);  // for GSL
    minimum->SetTolerance(0.0001);
    minimum->SetPrintLevel(0);
    
    // create function wrapper for minimizer
    // a IMultiGenFunction type
    ROOT::Math::Functor wrapf(Fmax,2);
    
    //variable = W, variable 1 = cosTheta
    double step[2] = {Wrange/100,2./100};
    // starting point at mid point in range
    double variable[2] = {Wmin+Wrange/2,0};
      
    minimum->SetFunction(wrapf);
    
    // Set the free variables to be minimized !
    minimum->SetVariable(0,"W",variable[0], step[0]);
    minimum->SetVariable(1,"cosTh",variable[1], step[1]);
 
    // do the minimization
    minimum->Minimize();
    const double *xs = minimum->X();
    
    auto minVal = minimum->MinValue();
    auto minW= xs[0];
    auto mint= kine::tFromcosthW(xs[1],minW,M1,M2,M3,M4);//xs[1];
    auto minCosTh=  xs[1];
    
    std::cout << "Maximum : Probabiltiy Dist at ( W=" << minW << " , t = "  << mint<< " , cosTh = "  << minCosTh << "): "<< -minimum->MinValue()  << " note t0 "<<kine::t0(minW,M1,M2,M3,M4)<< " and M3 "<<M3<< std::endl;
    
    //check for low mass meson limits
    if(dynamic_cast<DecayingParticle*>(_meson)){ //meson
      dynamic_cast<DecayingParticle*>(_meson)->TakeMinimumMass();//to get threshold behaviour
      M3=_meson->Mass();
	
      // do the minimization at min mass in case higher max
      minimum->Minimize();
      const double *xs = minimum->X();
	
      auto minminVal = minimum->MinValue();
	

      if(minminVal<minVal){
	std::cout << "Minimum Mass Maximum : Probabiltiy Dist at ( W=" << minW << " , t = "  << mint << "): "<< -minimum->MinValue()<< " note t0 "<<kine::t0(minW,M1,M2,M3,M4)  << std::endl;
	std::cout<<"minmin "<<-minminVal<<" < "<<-minVal<<std::endl;
	minVal=minminVal;
	minW= xs[0];
	mint= xs[1];
      }
      //back to PDg mass if exists
      if(_meson->PdgMass()>M3)
	M3=_meson->PdgMass();
	
      dynamic_cast<DecayingParticle*>(_meson)->TakePdgMass();

    }
    //make +ve again
    return -minVal;
  }
  ///////////////////////////////////////////////////////////////////////
  void TwoBodyProduction::HistIntegratedXSection(TH1D& hist){

 
    auto M1 = 0;//assume real photon for calculation
    auto M2 = _target->M();
    auto M3 = _meson->Mass(); //should be pdg value here
    auto M4 = _baryon->Mass();
    auto Wmin = M3+M4;
 
    //integrate over costh
    auto F = [this,M1,M2,M3,M4](double costh)
      {
	//if(_W<M3+M4) return 0.;
	
	//	_s=_W*_W;
	//	_t = kine::tFromcosthW(costh, _W, M1, M2, M3, M4);
  	//_cosThCM=costh;
  	return  dsigma_dcosth(get_W(),costh);
     };
    
      ROOT::Math::GSLIntegrator ig(ROOT::Math::IntegrationOneDim::kADAPTIVE,
  				   ROOT::Math::Integration::kGAUSS61);
      ROOT::Math::Functor1D wF(F);
      ig.SetFunction(wF);

      
  
      for(int ih=1;ih<=hist.GetNbinsX();ih++){
  	set_W(hist.GetXaxis()->GetBinCenter(ih));
  	if( _W < Wmin )
  	  hist.SetBinContent(ih, 0);
  	else
  	  hist.SetBinContent(ih, ig.Integral(-1,1) );
	
    }
      std::cout<<std::endl;
      //done
  }
  
  //////////////////////////////////////////////////////////
  void TwoBodyProduction::HistMaxXSection(TH1D& hist){
    hist = _histHighXS;
    return;
    
    auto M1 = 0;//assum real photon for calculation
    auto M2 = _target->M();
    auto M3 = _meson->Mass(); //should be pdg value here
    auto M4 = _baryon->Mass();
    //auto Wmin = M3+M4;
    auto Wmin = Parent()->MinimumMassPossible();
 
    auto F = [this,Wmin,M1,M2,M3,M4](double costh)
      {
	//_W=W;
	if(_W<Wmin)return 0.;
	_s=_W*_W;      
	_t = kine::tFromcosthW(costh, _W, M1, M2, M3, M4);
	_cosThCM=costh;
	return  DiffXS();
      };
    
   
      
 
      for(int ih=1;ih<=hist.GetNbinsX();ih++){
	_W=hist.GetXaxis()->GetBinCenter(ih);
	if( _W < Wmin )
	  hist.SetBinContent(ih, 0);
	else{
	  double max_at_W=0;
	  int Ntpoints=100;
	  double cth_step=(2)/Ntpoints;
	  double cth_val=-1;
	  for(int itt=0;itt<Ntpoints;itt++){
	    _W=hist.GetXaxis()->GetBinCenter(ih);
	    
	    double val_at_cth = F(cth_val);
	    if(val_at_cth>max_at_W)
	      max_at_W=val_at_cth;

	    
	    _W=_W+hist.GetXaxis()->GetBinWidth(ih)/2; //take right limit
	    val_at_cth = F(cth_val);
	    if(val_at_cth>max_at_W)
	      max_at_W=val_at_cth;
	    
	    _W=_W-hist.GetXaxis()->GetBinWidth(ih)/2; //take left limit 
	    val_at_cth = F(cth_val);
	    if(val_at_cth>max_at_W)
	      max_at_W=val_at_cth;
	    
	    //move on
	    cth_val+=cth_step;
	  }
	  hist.SetBinContent(ih, max_at_W );

	}
      }
    
  }

 //////////////////////////////////////////////////////////
 void TwoBodyProduction::Use_WCosTh_Envelope(){
   _createMyEnvelope=kTRUE;
 }
  //////////////////////////////////////////////////////////
  void TwoBodyProduction::Build_WCosTh_Envelope()  {
    
    std::cout<<"TwoBodyProduction::Build_WCosTh_Envelope"<<std::endl;
    auto M1 = 0;//assum real photon for calculation
    auto M2 = _target->M();
    auto M3 = _meson->Mass(); //should be pdg value here
    //auto M3 = _meson->MaximumMassPossible();// _meson->Mass(); //should be pdg value here
    auto M4 = _baryon->Mass();

    auto minW = Parent()->MinimumMassPossible();

    auto maxW = ( *(_target) + *(_ebeam) ).M();

    //Define cosTheta binning
    //Need more resolution close to cosTheta=1
    Int_t Ncth=100;
    Float_t convertToLog=18./Ncth;//18 as = range -10 to 10 (-2) due to the +1 in TMath::Log10. The makes the argument have a minimum value of 1 and log10(10)=0
    std::vector<double > cthBins;
    for(int ib=-Ncth/2;ib<=Ncth/2;++ib){
      auto sign = ib==0 ? 0 : (ib/TMath::Abs(ib));
      cthBins.push_back(sign*TMath::Log10(TMath::Abs(float(ib)*convertToLog)+1));
     }

    //W bins want focussed on threshold
    std::vector<double > WBins;
    int NW=60;
    double WRange = maxW-minW;
    double deltaW = WRange/NW;
    for(int iW=0;iW<NW;++iW){
      int Nsteps = NW-iW;
      for(int iWi=0;iWi<Nsteps;++iWi){
	WBins.push_back(minW + iW*deltaW+iWi*deltaW/Nsteps );
     }
    }

    //create envelope histogram
    TH2D hist("WCosTh","WCosTh",WBins.size()-1,WBins.data(),Ncth,cthBins.data());

    //consider cross section at 3 points.
    //minimum, PDG and maximum to make sure we get relaiable maximum
    //for full meson mass range
    auto M3s=std::vector<double>{_meson->Mass()};
    if(_meson->Mass()!=_meson->MinimumMassPossible()) M3s.push_back(_meson->MinimumMassPossible());
    if(_meson->Mass()!=_meson->MaximumMassPossible()) M3s.push_back(_meson->MaximumMassPossible());

    //Loop over meson masses
    for(int Mpos=0;Mpos<M3s.size();++Mpos){ //take max of upper,mid and lower bin values
      M3=M3s[Mpos];

      //jpac 2-body matrix elements squared probably dont depend on mass
      auto pdgM3 = _meson->Mass();

      //histogram with x-axis = W ; y-axis = cos(theta)
      for(int Wpos=-1;Wpos<1;++Wpos){ //take max of upper,mid and lower bin values
	for(int ih=1;ih<=hist.GetNbinsX();ih++){
	  _W=hist.GetXaxis()->GetBinCenter(ih) + Wpos*hist.GetXaxis()->GetBinWidth(ih)/2;
	  //create a correction to account for differences in phase space
	  //due to mass
	  auto PScorr=kine::PDK(_W,_meson->MinimumMassPossible(),_baryon->P4().M())/kine::PDK(_W,_meson->P4().M(),_baryon->P4().M());
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
	      _t = kine::tFromcosthW(_cosThCM, _W, M1, M2, M3, M4);
	      val=dsigma_dcosth();
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
	      
	      auto tmax = kine::tFromcosthW(cthmax, _W, M1, M2, M3, M4);
	      auto tmin = kine::tFromcosthW(cthmin, _W, M1, M2, M3, M4);
	      Double_t trange= tmax-tmin;
	      Int_t nsteps = 10;
	      Double_t maxvalatt =0;
	      //step in t closest to edge looking for highest value
	      for(UInt_t istep = 0;istep<nsteps;istep++){
		_t = tmax-istep*trange/nsteps/20;
		_cosThCM=kine::costhFromt(_t, _W,M1,M2,M3,M4);
		Double_t val=dsigma_dcosth();
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

	    
	}
	
      }
    }//Wpos
    }//Mpos
 
    //histogram fraction of mass distribution allowed as function of W
    //used to suppress sub-threshold cross sections
    DecayingParticle* decPart = dynamic_cast<DecayingParticle*>(_meson);
    if(decPart!=nullptr){
      TH1D histMassFrac("WMassFrac","WMassFrac",WBins.size()-1,WBins.data());
      for(int ih=1;ih<=histMassFrac.GetNbinsX();ih++){
	auto frac = decPart->IntegratedMass(histMassFrac.GetXaxis()->GetBinCenter(ih),_baryon->P4().M());
	histMassFrac.SetBinContent(ih,frac);
      }
      _distMassFrac.reset(new DistTH1{histMassFrac});
    }

   
    for(int ih=1;ih<=hist.GetNbinsX();ih++){
     auto maxVal= 0;
     for(int ic=1;ic<=hist.GetNbinsY();ic++){
       if(hist.GetBinContent(ih,ic)>maxVal)
	 maxVal=hist.GetBinContent(ih,ic);
     }
     for(int ic=1;ic<=hist.GetNbinsY();ic++){
       hist.SetBinContent(ih,ic,hist.GetBinContent(ih,ic)+maxVal*0.02 );
     }
   }
   //adjust max after addition
   auto  maxhist = hist.GetMaximum();
   _max=maxhist*1.1;

   //_histHighXS = TH1D{"high2BodyXSvW","high2BodyXSvW",hist.GetXaxis()->GetNbins(),hist.GetXaxis()->GetXbins()->GetArray()};
   
    // for(int ih=1;ih<=hist.GetNbinsX();ih++){
    //   Double_t maxXs=0.0;
    //   Double_t integral =0.0;
    //   for(int ic=1;ic<=hist.GetNbinsY();ic++){
    // 	auto currXs = hist.GetBinContent(ih,ic);
    // 	integral+=currXs*hist.GetXaxis()->GetBinWidth(ic);
    // 	if(currXs>maxXs)maxXs=currXs;
    //   }
    //   std::cout<<"max "<<maxXs<<" at "<<hist.GetXaxis()->GetBinCenter(ih)<<" integral "<<integral<<std::endl;
    //   _histHighXS.SetBinContent(ih,maxXs);
    // }
    // _distHighXS.reset( new DistTH1(_histHighXS) );

   //helps to appply some smoothing to histogram!
   hist.Smooth(2);
   //apply envelope if requested
   if(_createMyEnvelope)Parent()->SetDecayer(new TwoBodyEnvelope(DistTH2Slice{hist}));

   
  }


}
