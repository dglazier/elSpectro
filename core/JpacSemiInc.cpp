#include "JpacSemiInc.h"
#include <TDatabasePDG.h>
#include <Math/GSLIntegrator.h>
#include <Math/IntegrationTypes.h>
#include <Math/Functor.h>
#include <Math/Minimizer.h>
#include <Math/Factory.h>

namespace elSpectro{
  ///////////////////////////////////////////////////////
  ///constructor includes subseqent decay of Ngamma* system
  JpacSemiInc::JpacSemiInc( jpacInc_ptr inc ,
			    particle_ptrs parts, const std::vector<int> pdgs) :
    _inc{inc},
    TwoBodyProduction{ parts, pdgs }
  {
    _name={"JpacSemiInc"};
    _minMassB= GetBaryon()->MinimumMassPossible();
  }

  //////////////////////////////////////////////////////
  void JpacSemiInc::PostInit(ReactionInfo* info) {
    _prodInfo= dynamic_cast<ReactionElectroProd*> (info); //I need Reaction info
     if(_createMassEnvelope==true) CreateMassEnvelope();
    TwoBodyProduction::PostInit(info);
 
    
  }

  //////////////////////////////////////////////////////
  double JpacSemiInc::dsigma_dcosth(double W,double costheta) const {
    set_cosThCM(costheta);
    set_W(W);
    return dsigma_dcosth();
  }
  //////////////////////////////////////////////////////
  double JpacSemiInc::dsigma_dcosth() const {
    //integrate over r
    
    //below threshold
    if( get_W() < (GetMeson()->P4().M()+GetBaryon()->P4().M()) ) return 0;

    //  std::cout<<"JpacSemiInc::dsigma_dcosth "<<std::endl;
    // How we integrate depends on what variables are being used
    double result = 0.;
    ROOT::Math::GSLIntegrator ig(ROOT::Math::IntegrationOneDim::kADAPTIVE, ROOT::Math::Integration::kGAUSS15);

    // integrate over usual t bounds
    auto dSigma = [&](double r)
    {
      _r=r;
      //returns d2sigma/drdcth
      //_s and _cosThetaCM should be set
      // std::cout<<"JpacSemiInc::dsigma_dcosth dsigma( "<<get_s()<<" "<<get_r()<<" "<<get_cosThCM()<<std::endl;
      return DiffXSFromRCosThatQ20();
    };
    
    ROOT::Math::Functor1D wF(dSigma);
    ig.SetFunction(wF);

    result = ig.Integral(0,1); //r defined in range 0-1

    //std::cout<<"JpacSemiInc::dsigma_dcosth = "<<result<< " at "<<get_cosThCM()<<std::endl;
    if( TMath::IsNaN(result) ) return 0.;
    if( TMath::Infinity()==result ) return 0.;

    return result;
  }
    //////////////////////////////////////////////////////////
  void JpacSemiInc::HistMaxXSection(TH1D& hist){
    
    
    // auto M1 = 0;//assum real photon for calculation
    // auto M2 = GetTarget()->M();
    // auto M3 = GetMeson->Mass(); //should be pdg value here
    //auto M4 = _baryon->Mass();
    //auto Wmin = M3+M4;
    auto Wmin = Parent()->MinimumMassPossible();
    std::cout<<"JpacSemiInc::HistMaxXSection "<<std::endl;
      
    auto F = [this,Wmin]()
      {
	if(get_W()<Wmin)return 0.;
	//if(TMath::Sqrt(get_MB2())>2) return 0.;
	set_s(get_W()*get_W());      
	auto val = DiffXSFromRCosThatQ20();
	if( TMath::IsNaN(val) ) return 0.;
	if( TMath::Infinity()==val ) return 0.;
	//std::cout<<"JpacSemiInc::HistMaxXSection f( "<<get_s()<<" "<<get_r()<<" "<<get_cosThCM()<<" MB2 "<<get_MB2()<<" MB "<<TMath::Sqrt(get_MB2())<< " val "<<val<<std::endl;
	return val;
      };

    //store global maximum and position
    double max_global=0.;
    double ictmax=0;
    double irmax=0;
    double iWmax=0;

    double max_at_W=0;
    int Ntpoints=10;
    double cth_step=(0.01)/Ntpoints;
    int Nrpoints=100;
    double r_step=(1.)/Nrpoints;

	//Loop over given histogram bins
    //and find max xsection for that W
    for(int ih=1;ih<=hist.GetNbinsX();ih++){
      //for(int ih=1;ih<10;ih++){
      set_W(hist.GetXaxis()->GetBinCenter(ih));
      max_at_W=0;
      if( get_W() < Wmin )
	hist.SetBinContent(ih, 0);
      else{

	auto stepAction =[this,&max_at_W,&ictmax,&irmax,&iWmax,&max_global,&F](int ict,int ir){
			   double val_here = F();
			   //std::cout<<" xs "<<val_here<<std::endl;
			   if(val_here>max_at_W){
			     max_at_W=val_here;
			     if(val_here>max_global){
			       max_global = val_here;
			       ictmax = get_cosThCM();
			       irmax=get_r();
			       iWmax=get_W();
			     }
			     
			   }

			 };

	//set_cosThCM(-1.);//reset cosTheta
	set_cosThCM(0.99);//reset cosTheta
	//loop over cosTh and r and evaluate cross section to find max.

	//Cross section highest close to threshold in Mbaryon (i.e. Delta), make only consider upto max Mbaryon
	//std::cout<<"r_step = "<<r_step<<" "<<std::endl;
	auto rmin = RFromM2(1.4); //Make this a parameter for changing
	r_step=(1-rmin)/Nrpoints; //max 1, where mass is smallest

	//std::cout<<"r_step = "<<r_step<<" "<<rmin<<std::endl;
	for(int itt=0;itt<Ntpoints;itt++){
	  _r = rmin;//reset r
	  
	  for(int ir=0;ir<Nrpoints;ir++){
	    
	    set_W(hist.GetXaxis()->GetBinCenter(ih));

	    stepAction(itt,ir);
	    set_W(get_W()+hist.GetXaxis()->GetBinWidth(ih)/2); //take upper limit
	    stepAction(itt,ir);
	    set_W(get_W()-hist.GetXaxis()->GetBinWidth(ih)); //take lower limit 
	    stepAction(itt,ir);

	    _r+=r_step;
	  }
	  set_cosThCM(get_cosThCM()+cth_step);
	}
	
	//found the max for this W, store it in histogram
	hist.SetBinContent(ih, max_at_W );
	std::cout<<"JpacSemiInc::HistMaxXSection W bin "<< ih <<" out of "<<hist.GetNbinsX()<<" with max "<<max_at_W<<std::endl;
 
      }//end above threshold action
      
    }//end W loop

    //could now scan around binned max to find global max...
    std::cout<<"JpacSemiInc::HistMaxXSection max value of cross section = "<<max_global<<" at cos(theta),r,W = "<<ictmax<<","<<irmax<<","<<iWmax<<std::endl;


    set_max( max_global);
    //auto newmax = FindMaxOfIntensity();
    //if(newmax>get_max())set_max(newmax);

  }
  //////////////////////////////////////////////////////////////////////////

  double JpacSemiInc::FindMaxOfIntensity(){
    std::cout<<"JpacSemiInc::FindMaxOfIntensity()************************************************"<<std::endl;
    auto Wmin = Parent()->MinimumMassPossible();
   
    auto Wmax = ProductionInfo()->_Wmax;
    if(Wmax == 0){
      std::cerr<<" JpacSemiInc:::FindMaxOfIntensity(), Wmax=0, should be set in ReactionInfo"<<std::endl;
      exit(0);
    }
    double Wrange=Wmax-Wmin;
    
    auto F = [this,Wmin](const double *x)
      {
	set_W( x[0]);
	set_cosThCM(x[1]);
	set_r(x[2]);
	
	if(get_W()<Wmin)return 0.;

	auto val = DiffXSFromRCosThatQ20();
	if( TMath::IsNaN(val) ) return 0.;
	if( TMath::Infinity()==val ) return 0.;
	//	std::cout<<"JpacSemiInc::FindMaxOfIntensity() f( "<<get_s()<<" "<<get_r()<<" "<<get_cosThCM()<<" MB2 "<<get_MB2()<<" MB "<<TMath::Sqrt(get_MB2())<< " val "<<val<<std::endl;
	return -val;
      };

    //First try with nominal (probably PDG meson and baryon mass)
    ROOT::Math::Minimizer* minimum =
      ROOT::Math::Factory::CreateMinimizer("Genetic", "");
    
    
    // set tolerance , etc...
    minimum->SetMaxFunctionCalls(1000000); // for Minuit/Minuit2
    minimum->SetMaxIterations(1000);  // for GSL
    minimum->SetTolerance(0.0001);
    minimum->SetPrintLevel(0);
    
    // create function wrapper for minimizer
    // a IMultiGenFunction type
    ROOT::Math::Functor wrapf(F,2);
    
    //variable = W, variable 1 = cosTheta, variable 2 = r
    double step[3] = {Wrange/100,2./100,1./100};
    // starting point at mid point in range, needs to be mid point for genetic
    double variable[3] = {Wmin+Wrange/2,0,0.5}; //low mass (Delta) => high r
      
    minimum->SetFunction(wrapf);
    
    // Set the free variables to be minimized !
    minimum->SetVariable(0,"W",variable[0], step[0]);
    minimum->SetVariable(1,"cosTh",variable[1], step[1]);
    minimum->SetVariable(2,"r",variable[2], step[2]);
 
    // do the minimization
    minimum->Minimize();
    const double *xs = minimum->X();
    
    auto minVal = minimum->MinValue();
    auto minW= xs[0];
    auto minCosTh=  xs[1];
    auto minr=  xs[2];
    auto mint= TFromWRCosThatQ20(minW,minr,minCosTh);;

    std::cout << "Maximum : Probabiltiy Dist at ( W=" << minW << " , t = "  << mint<< " , cosTh = "  << minCosTh<< " , r = "  << minr << "): "<< -minimum->MinValue()  << " note t0 "<< std::endl;

    return get_max();
   }
  
  ////////////////////////////////////////////////////////////////////////////////////////////////////
  void JpacSemiInc::CreateMassEnvelope(){
    auto oldDist = GetBaryon()->MassDistribution();
    auto hist = TH1D("baryonMass","baryonMass",200,oldDist->GetMinX(),oldDist->GetMaxX());

    auto xsFunc = [&](){
			set_W(GetBaryon()->Mass()+GetMeson()->Mass());
			auto xs = _inc->dsigma_dM2(get_s(),get_MB2() );
			set_W(ProductionInfo()->_Wmax/2);
 			xs += _inc->dsigma_dM2(get_s(),get_MB2() );
			set_W(ProductionInfo()->_Wmax);
			xs += _inc->dsigma_dM2(get_s(),get_MB2() );
			return xs;
			
    };
    auto setMB = [this](double m){set_MB2(m*m);};
      
    for(int ibin=0; ibin<=hist.GetNbinsX();++ibin){
      //evaluate differential cross section at maximum W
      set_W(ProductionInfo()->_Wmax);
      set_W(GetBaryon()->Mass()+GetMeson()->Mass()+2);
      setMB(hist.GetXaxis()->GetBinCenter(ibin));
      
      //std::cout<<"JpacSemiInc::CreateMassEnvelope() bin = "<<ibin<<" M="<<get_MB2()<< " xs = "<<get_W()<<std::endl;
      auto xsm = xsFunc() ;
      setMB(hist.GetXaxis()->GetBinCenter(ibin)- hist.GetXaxis()->GetBinWidth(ibin)/2);
      auto xsl = xsFunc() ;
      setMB(hist.GetXaxis()->GetBinCenter(ibin)+ hist.GetXaxis()->GetBinWidth(ibin)/2);
      auto xsh = xsFunc() ;

      auto xs = xsm>xsl ? xsm:xsl;
      xs = xsh>xs ? xsl:xs;
      std::cout<<"JpacSemiInc::CreateMassEnvelope() bin = "<<ibin<<" M="<<get_MB2()<< " xs = "<<xs<<std::endl;
      hist.SetBinContent(ibin, xs );
    }
    _massEnvelope.reset(new DistTH1{hist});
    std::cout<<" JpacSemiInc::CreateMassEnvelope()  min mass possible "<<_massEnvelope->GetMinX()<<std::endl;
   const_cast<Particle*>(GetBaryon())->SetMassDist((_massEnvelope).get());
    const_cast<Particle*>(GetBaryon())->SetPdgMass(_massEnvelope->GetMinX());
   // exit(0);

  }
}
