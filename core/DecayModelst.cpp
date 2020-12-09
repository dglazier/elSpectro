#include "DecayModelst.h"
#include "FunctionsForGenvector.h"
#include <TDatabasePDG.h>
#include <Math/GSLIntegrator.h>
#include <Math/IntegrationTypes.h>
#include <Math/Functor.h>
#include <Math/Minimizer.h>
#include <Math/Factory.h>
#include <TRandom3.h>

namespace elSpectro{
  ///////////////////////////////////////////////////////
  ///constructor includes subseqent decay of Ngamma* system
  DecayModelst::DecayModelst(particle_ptrs parts, const std::vector<int> pdgs) :
    DecayModel{ parts, pdgs }
  {
    _name={"DecayModelst"};

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
  void DecayModelst::PostInit(ReactionInfo* info){

    if( dynamic_cast<DecayingParticle*>(_meson) ){
      if( dynamic_cast<DecayingParticle*>(_meson)->Model()->CanUseSDME() ){
	_sdmeMeson = _meson->InitSDME(1,4);
	//could have electroproduced baryon spin 3/2
	//_sdmeBaryon = _baryon->InitSDME(3,9);
      }
    }
    
    DecayModel::PostInit(info);
    
    _prodInfo= dynamic_cast<ReactionElectroProd*> (info); //I need Reaction info

    _photon = _prodInfo->_photon;
    _target = _prodInfo->_target;
    _ebeam = _prodInfo->_ebeam;
    _photonPol = _prodInfo->_photonPol;
    
     double maxW = ( *(_prodInfo->_target) + *(_prodInfo->_ebeam) ).M();

     _max = FindMaxOfIntensity(); //add 20% for Q2,meson mass effects etc.

     std::cout<<"DecayModelst::PostInit max value "<<_max<<" "<<_meson<<" "<<_meson->Pdg()<<" "<<_sdmeMeson<<std::endl;
  }
  
  //////////////////////////////////////////////////////////////////
  double DecayModelst::Intensity() const
  {
    /*A      B        A/2    B/2        A*2/3   B   
      1      2         1/2    1          1/3    1   30each    5     10       
      1/2    1   -->   1/4    1/2  --->  1/6    1/2  ---->    5/2   5
      0      2         0      1                 1   row       0     10

     */
    _W = Parent()->P4().M();
    _s=_W*_W;
    _t = (_meson->P4()-*_photon).M2();//_amp->kinematics->t_man(s,cmMeson.Theta());
       //now we can define production/polarisation plane
    MomentumVector decayAngles;
    kine::electroCMDecay(&Parent()->P4(),_ebeam,_photon,_meson->P4ptr(),&decayAngles);
    _photonPol->SetPhi(decayAngles.Phi());

    //   std::cout<<_W <<" and "<<_t<<" "<<_meson->P4().E()<<" "<<_photon->E()<<" cos costheta "<<TMath::Cos(decayAngles.Theta())<<std::endl;

   //weight==differential cross section //See Seyboth and Wolf eqn (70)
   // double trangeRatio=4* kine::PDK(_W,0,_target->M())  * kine::PDK(_W,_meson->P4().M(),_baryon->P4().M() );// / (  kine::PDK(_Wmax,0,_target->M() )* kine::PDK(_Wmax,_meson->P4().M(),_baryon->P4().M() ) );
    _dt=4* TMath::Sqrt(PgammaCMsq())  * kine::PDK(_W,_meson->P4().M(),_baryon->P4().M() );
    //double otherdt=4*kine::PDK(_W,_meson->P4().M(),_baryon->P4().M() )   * kine::PDK(_W,0,_baryon->P4().M() );
    // std::cout<<" dt "<<  _dt<<" "<<otherdt<<std::endl;
    
    double weight = DifferentialXSect() * _dt ;

    // _dt=otherdt;
    //_dt/=(2*(TMath::Sqrt(PgammaCMsq())  * kine::PDK(_W,_meson->P4().M(),_baryon->P4().M() )));
    CalcMesonSDMEs();
    CalcBaryonSDMEs();
 
    weight/=_max; //normalise range 0-1
    weight/=TMath::Sqrt(PgammaCMsq()/kine::PDK2(_W,0,_baryon->Mass())); //correct max for finite Q2 phase space
    

    if(weight>1){
      //don't change weight, likely due to large Q2 value....
      
      //auto oldmax=_max;
      //_max=weight*oldmax;
      std::cout<<"DecayModelst::Intensity weight too high but won't change maxprobable low meson mass and W from  "<<_max<<" to "<<weight*_max<<" meson "<<_meson->Mass()<<" W "<<_W<<std::endl;
      //std::cout<<"DX "<<DifferentialXSect()<<" "<< _dt<<" pgam "<<PgammaCMsq()<<" M^2 "<<MatrixElementsSquared_T()<<" Q2 factor "<<PgammaCMsq()/kine::PDK(_W,_meson->Mass(),_baryon->Mass())<<std::endl;
    }
    
    //Correct for W weighting which has already been applied
    weight/=_prodInfo->_sWeight;
    
    if(weight>1){
      std::cout<<" s weight "<<_prodInfo->_sWeight<<" Q2 "<<-_photon->M2()<<" 2Mmu "<<2*_target->M()*_photon->E() <<" W "<<_W<<" t "<<_t<<" new weight "<<weight*_prodInfo->_sWeight<<" meson "<<_meson->Mass()<<std::endl;
      std::cout<<"DecayModelst::Intensity sWeight corrected weight too large "<<weight <<" "<<_prodInfo->_sWeight<<"  max "<<_max<<" val "<< weight*_prodInfo->_sWeight*_max<<std::endl;
      std::cout<<"DX "<<DifferentialXSect()<<" "<< _dt<<" pgam "<<TMath::Sqrt(PgammaCMsq())<<" M^2 "<<MatrixElementsSquared_T()<<" Q2 factor "<<TMath::Sqrt(PgammaCMsq())/kine::PDK(_W,0,_baryon->Mass())<<std::endl;
    //std::cout<<"flux "<<kine::FluxPhaseSpaceFactor(*_photon,*_target)<<" "<<4*kine::PDK(W,_meson->Mass(),_baryon->Mass())*W<<" "<< kine::PhaseSpaceFactorDt(W,P1,_meson->Mass(),_baryon->Mass())<<std::endl;				    
      //std::cout<<"Pi check "<<P1 <<" versus "<< kine::PDK(W,_photon->M(),_target->M())<<std::endl;
    }
     
 
    return weight;
    
  }
  
  double DecayModelst::FindMaxOfIntensity(){
    
    auto M1 = 0;//assum real photon for max calculation
    auto M2 = _target->M();
    auto M3 = _meson->Mass(); //should be pdg value here
    auto M4 = _baryon->Mass();
    // auto Wmin = M3+M4;
    auto Wmin = Parent()->MinimumMassPossible();
  
    _Wmax = ( *(_prodInfo->_target) + *(_prodInfo->_ebeam) ).M();

    std::cout<<" DecayModelst::FindMaxOfIntensity()  "<<Wmin<<" "<<_Wmax<<std::endl;

    auto Fmax = [&M1,&M2,&M3,&M4,&Wmin,this](const double *x)
	{
	  _s = x[0]*x[0];
	  _W=x[0];
	  if( _W < Wmin ) return 0.;
	  if( _W < M3+M4 ) return 0.;
	  if( x[1] > kine::t0(_W,M1,M2,M3,M4) ) return 0.;
	  if( x[1]< kine::tmax(_W,M1,M2,M3,M4) ) return 0.;
	  _t=x[1];
	  double val = DifferentialXSect()* (kine::t0(_W,M1,M2,M3,M4)-kine::tmax(_W,M1,M2,M3,M4));
	  
	  if( TMath::IsNaN(val) ) return 0.;
	  return -val; //using a minimiser!
	};
      
      //First perform grid search for intital values
      double Wrange=_Wmax-Wmin;
      double tmax=kine::tmax(_Wmax,M1,M2,M3,M4);
       
      double gridMin=0;
      double gridW=0;
      double gridt=0;
      double WtVals[2];
      int Npoints=50;
      for(int iW=1;iW<Npoints;iW++){
	WtVals[0]=Wmin+iW*Wrange/Npoints;
	double tming=kine::t0(WtVals[0],M1,M2,M3,M4);
	double tmaxg=kine::tmax(WtVals[0],M1,M2,M3,M4);
	double trange= tming-tmaxg;

	for(int it=0;it<Npoints;it++){
	  WtVals[1]=tming-it*trange/Npoints;
	  auto val = Fmax(WtVals);
	  // if(it==0)std::cout<<WtVals[0]<<" "<<WtVals[1]<<" val "<<val<<std::endl;
	  if(val<gridMin) {
	    gridMin=val;
	    gridW=WtVals[0];
	    gridt=WtVals[1];
	  }
	}
      }

      std::cout<<"FindMaxOfProbabilityDistribution grid search max= "<<-gridMin<<" at W = "<<gridW<<" and t = "<<gridt<<std::endl;
      ROOT::Math::Minimizer* minimum =
	ROOT::Math::Factory::CreateMinimizer("Minuit2", "");

      if(minimum==nullptr) //Minuit2 not always installed!
	minimum = ROOT::Math::Factory::CreateMinimizer("Minuit", "");
      
      // set tolerance , etc...
      minimum->SetMaxFunctionCalls(1000000); // for Minuit/Minuit2
      minimum->SetMaxIterations(10000);  // for GSL
      minimum->SetTolerance(0.0001);
      minimum->SetPrintLevel(1);
      
      // create function wrapper for minimizer
      // a IMultiGenFunction type
      ROOT::Math::Functor wrapf(Fmax,2);

      //variable = W, variable 1 = t
      double step[2] = {0.1,0.1};
      // starting point
      
      double variable[2] = { gridW,gridt};
      
      minimum->SetFunction(wrapf);
      
      // Set the free variables to be minimized !
      minimum->SetVariable(0,"W",variable[0], step[0]);
      minimum->SetVariable(1,"t",variable[1], step[1]);
      minimum->SetVariableLimits(0,Wmin,_Wmax);
      minimum->SetVariableLimits(1,tmax,0);
      
      // do the minimization
      minimum->Minimize();
      const double *xs = minimum->X();
 
      auto minVal = minimum->MinValue();
      auto minW= xs[0];
      auto mint= xs[1];
   
      std::cout << "Maximum : Probabiltiy Dist at ( W=" << minW << " , t = "  << mint << "): "<< -minimum->MinValue()  << std::endl;

      //check for low mass meson limits
      if(dynamic_cast<DecayingParticle*>(_meson)){ //meson
	dynamic_cast<DecayingParticle*>(_meson)->TakeMinimumMass();//to get threshold behaviour
	M3=_meson->Mass();
	
	// do the minimization at min mass in case higher max
	minimum->Minimize();
	const double *xs = minimum->X();
	
	auto minminVal = minimum->MinValue();
	minW= xs[0];
	mint= xs[1];
	
	std::cout << "Maximum : Probabiltiy Dist at ( W=" << minW << " , t = "  << mint << "): "<< -minimum->MinValue()  << std::endl;

	if(minminVal<minVal)minVal=minminVal;
	//back to PDg mass if exists
	if(_meson->PdgMass()>M3)
	  dynamic_cast<DecayingParticle*>(_meson)->TakePdgMass();

      }

      
      return -minVal ;
  }

  void DecayModelst::HistIntegratedXSection(TH1D& hist){

 
    auto F = [this](double t)
      {
	//_W=W;
	_s=_W*_W;
	_t=t;
	return DifferentialXSect();
      };
    
      ROOT::Math::GSLIntegrator ig(ROOT::Math::IntegrationOneDim::kADAPTIVE,
				   ROOT::Math::Integration::kGAUSS61);
      ROOT::Math::Functor1D wF(F);
      ig.SetFunction(wF);

      
      auto M1 = 0;//assum real photon for calculation
      auto M2 = _target->M();
      auto M3 = _meson->Mass(); //should be pdg value here
      auto M4 = _baryon->Mass();
      auto Wmin = M3+M4;

      for(int ih=1;ih<=hist.GetNbinsX();ih++){
	_W=hist.GetXaxis()->GetBinCenter(ih);
	//	_W=_W+hist.GetXaxis()->GetBinWidth(ih); //take right limit so do not miss threshold
	if( _W < Wmin )
	  hist.SetBinContent(ih, 0);
	if( TMath::IsNaN(kine::tmax(_W,M1,M2,M3,M4)) )
	  hist.SetBinContent(ih, 0);
	if( TMath::IsNaN(kine::t0(_W,M1,M2,M3,M4)) )
	  hist.SetBinContent(ih, 0);
	else
	  hist.SetBinContent(ih, ig.Integral(kine::tmax(_W,M1,M2,M3,M4),kine::t0(_W,M1,M2,M3,M4)) );
	
	std::cout<<ih<<" "<<"Going to integrate from t "<<kine::tmax(_W,M1,M2,M3,M4)<<" "<<kine::t0(_W,M1,M2,M3,M4)<<" at W "<<" and "<<_W<<" "<<Wmin<<" with result "<<hist.GetBinContent(ih)<<"t rage "<< kine::tmax(_W,M1,M2,M3,M4) - kine::t0(_W,M1,M2,M3,M4)<<std::endl;
	//or try histogram method
	TH1F histi("histi","integral",100, kine::tmax(_W,M1,M2,M3,M4),kine::t0(_W,M1,M2,M3,M4));
	for(int it=1;it<=histi.GetNbinsX();it++){histi.SetBinContent(it,F(histi.GetBinCenter(it)));}
	std::cout<<"alternative "<<histi.Integral("width")<<std::endl;
	if(ih%10==0)std::cout<<(hist.GetNbinsX() - ih)/10<<" "<<std::endl;
      }
      std::cout<<std::endl;
      //done
  }
  
  void DecayModelst::HistMaxXSection(TH1D& hist){

 
    auto M1 = 0;//assum real photon for calculation
    auto M2 = _target->M();
    auto M3 = _meson->Mass(); //should be pdg value here
    auto M4 = _baryon->Mass();
    //auto Wmin = M3+M4;
    auto Wmin = Parent()->MinimumMassPossible();
 
    auto F = [this,&Wmin](double t)
      {
	//_W=W;
	if(_W<Wmin)return 0.;
	_s=_W*_W;
	_t=t;
	return DifferentialXSect();
      };
    
   
      
 
      for(int ih=1;ih<=hist.GetNbinsX();ih++){
	_W=hist.GetXaxis()->GetBinCenter(ih);
	if( _W < Wmin )
	  hist.SetBinContent(ih, 0);
	if( TMath::IsNaN(kine::tmax(_W,M1,M2,M3,M4)) )
	  hist.SetBinContent(ih, 0);
	if( TMath::IsNaN(kine::t0(_W,M1,M2,M3,M4)) )
	  hist.SetBinContent(ih, 0);
	else{
	  double max_at_W=0;
	  double tmax=kine::tmax(_W,M1,M2,M3,M4);
	  double tmin=kine::t0(_W,M1,M2,M3,M4);
	  int Ntpoints=100;
	  double tstep=(tmax-tmin)/Ntpoints;
	  double tval=tmin;
	  for(int itt=0;itt<Ntpoints;itt++){
	    _W=hist.GetXaxis()->GetBinCenter(ih);

	    double val_at_t = F(tval)*(tmin-tmax);
	    if(val_at_t>max_at_W)
	      max_at_W=val_at_t;

	    
	    _W=_W+hist.GetXaxis()->GetBinWidth(ih)/2; //take right limit
	    val_at_t = F(tval)*(tmin-tmax);
	    if(val_at_t>max_at_W)
	      max_at_W=val_at_t;
	    
	    _W=_W-hist.GetXaxis()->GetBinWidth(ih); //take left limit 
	    val_at_t = F(tval)*(tmin-tmax);
	    if(val_at_t>max_at_W)
	      max_at_W=val_at_t;
	    
	    //move on
	    tval+=tstep;
	  }
	  hist.SetBinContent(ih, max_at_W );

	}
      }
      //	if(ih%10==0)std::cout<<(hist.GetNbinsX() - ih)/10<<" "<<std::endl;
  
      std::cout<<std::endl;
      //done
  }
}
/* perhaps this can go in script for fixed values of sdmes
    //Meson spin density marix elements, note this is photoproduced
    if(_sdmeMeson){
      _sdmeMeson->SetElement(0,0,0,);
      _sdmeMeson->SetElement(0,1,0,);
      _sdmeMeson->SetElement(0,1,-1,);
      _sdmeMeson->SetElement(1,1,1,);
      _sdmeMeson->SetElement(1,0,0,);
      _sdmeMeson->SetElement(1,1,0,);
      _sdmeMeson->SetElement(1,1,-1,);
      _sdmeMeson->SetElement(2,1,0,);
      _sdmeMeson->SetElement(2,1,-1,);
      //_sdmeMeson->SetElement(3,1,0,(_amp->SDME(3, 1, 0, s, t)));
      // _sdmeMeson->SetElement(3,1,-1,(_amp->SDME(3, 1, -1, s, t)));
    }
    */
