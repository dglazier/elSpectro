#include "DecayModelstM2.h"
#include "SDMEDecay.h"
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
  DecayModelstM2::DecayModelstM2(particle_ptrs parts, const std::vector<int> pdgs) :
    DecayModelst{ parts, pdgs }
  {
    _name={"DecayModelstM2"};

    if(Products().size()!=2){
      Fatal("DecayModelstM2","Can only have two decay particles");
    }
    // //need to find meson and baryon
    // if(TDatabasePDG::Instance()->GetParticle(Products()[0]->Pdg())->ParticleClass()==TString("Baryon") ){
    //   _baryon=Products()[0];
    //   _meson=Products()[1]; 
    // }
    // else {
    //   _baryon=Products()[1];
    //   _meson=Products()[0];
    // }
   _minMB= _baryon->MinimumMassPossible();
  }
  /////////////////////////////////////////////////////////////////
  void DecayModelstM2::PostInit(ReactionInfo* info){
   
     
    DecayModelst::PostInit(info);
 
    _isElProd=kTRUE;
    _prodInfo= dynamic_cast<ReactionElectroProd*> (info); //I need Reaction info
    if(_prodInfo==nullptr){
      _isElProd=kFALSE;
      _prodInfo= dynamic_cast<ReactionPhotoProd*> (info);
    }
    _photon = _prodInfo->_photon;
    _target = _prodInfo->_target;
    _ebeam = _prodInfo->_ebeam;
    _photonPol = _prodInfo->_photonPol;
    
     double maxW = ( *(_prodInfo->_target) + *(_prodInfo->_ebeam) ).M();

     _max = FindMaxOfIntensity()*1.08; //add 5% for Q2,meson mass effects etc.

     std::cout<<"DecayModelstM2::PostInit max value "<<_max<<" "<<_meson<<" "<<_meson->Pdg()<<" "<<_sdmeMeson<<std::endl;
  }
  
  //////////////////////////////////////////////////////////////////
  double DecayModelstM2::Intensity() const
  {
    /*A      B        A/2    B/2        A*2/3   B   
      1      2         1/2    1          1/3    1   30each    5     10       
      1/2    1   -->   1/4    1/2  --->  1/6    1/2  ---->    5/2   5
      0      2         0      1                 1   row       0     10

     */
    _W = Parent()->P4().M();
    _s=_W*_W;
    _t = (_meson->P4()-*_photon).M2();//_amp->kinematics->t_man(s,cmMeson.Theta());
    auto M2= get_M2();
    
    _dt=0;
    _dsigma=0;

    //check above threshold for meson and baryon masses
    if( _W < (_meson->P4().M()+_baryon->P4().M()) ) return 0;
    //std::cout<<"DecayModelstM2 "<<Parent()->Pdg()<<" "<<_meson->P4().M()<<" "<<_baryon->P4().M()<<std::endl;

    // if(_isElProd==kTRUE){
    //   //now we can define production/polarisation plane
    //   MomentumVector decayAngles;
    //   kine::electroCMDecay(&Parent()->P4(),_ebeam,_photon,_meson->P4ptr(),&decayAngles);
    //   _photonPol->SetPhi(decayAngles.Phi());
    // }
    // else{//photoproduction
    //   _photonPol->SetPhi(_meson->P4().Phi());
    // }

    // //now kinemaics
    // _dt=4* TMath::Sqrt(PgammaCMsq())  * kine::PDK(_W,_meson->P4().M(),_baryon->P4().M() );

    //do we need to multiply by dM2 here too ?
    double weight = DifferentialXSect();// * _dt ;//must multiply by t-range for correct sampling
  

    weight/=_max; //normalise range 0-1
    if(_isElProd==kTRUE)
      weight/= TMath::Sqrt(PgammaCMsq()/kine::PDK2(_W,0,_target->M())); //correct max for finite Q2 phase space
 
    if(weight>1){
      //don't change weight, likely due to large Q2 value....
      std::cout<<"DecayModelstM2::Intensity weight too high but won't change maxprobable low meson mass and W from  "<<_max<<" to "<<weight*_max<<" meson "<<_meson->Mass()<<" W "<<_W<<std::endl;
      }
    
    //Correct for W weighting which has already been applied
    weight/=_prodInfo->_sWeight;
    //std::cout<<" s weight "<<_prodInfo->_sWeight<<" weight "<<weight<<" "<<_W<<std::endl;
   
    return weight;
    
  }
  
  double DecayModelstM2::FindMaxOfIntensity(){
    
    auto M1 = 0;//assum real photon for max calculation
    auto M2 = _target->M();
    auto M3 = _meson->Mass(); //should be pdg value here
    auto M4 = _baryon->Mass();
    auto Wmin = Parent()->MinimumMassPossible();
   
    //    _Wmax = ( *(_prodInfo->_target) + *(_prodInfo->_ebeam) ).M();
    _Wmax = _prodInfo->_Wmax;
    
    auto Fmax = [&M1,&M2,&M3,&M4,&Wmin,this](const double *x)
      {
	_W = x[0];
	_s = _W*_W;
	auto X = x[1];
	auto Y = x[2];
	
	_M2 =  M2fromXY(X,Y);
	M4=TMath::Sqrt(_M2);
	
	if( _W < Wmin ) return 0.;
	if( _W < M3+M4 ) return 0.;
	if( _W > _Wmax ) return 0.;

	auto currt=kine::TFromXY(X,Y);
	auto myt0=kine::t0(_W,M1,M2,M3,M4);
	auto mytmax=kine::tmax(_W,M1,M2,M3,M4);
	
	if(currt>myt0) return 0.;
	if(currt<mytmax) return 0.;
	if( TMath::IsNaN(x[1]) ) return 0.;

	_t=currt;
	//auto dt=4* TMath::Sqrt(PgammaCMsq())  * kine::PDK(_W,M3,M4 );

	double val = DifferentialXSect();//*dt;
	if( TMath::IsNaN(val) ) return 0.;
	return -(val); //using a minimiser!
      };
   
      //First perform grid search for intital values
      double Wrange=_Wmax-Wmin;
      double tmax=kine::tmax(_Wmax,M1,M2,M3,M4);
       
      double gridMin=0;
      double gridW=0;
      double gridt=0;
      double WtVals[2];
   
      ROOT::Math::Minimizer* minimum =
	ROOT::Math::Factory::CreateMinimizer("Genetic", "");
      //	ROOT::Math::Factory::CreateMinimizer("Minuit2", "Combined");

      if(minimum==nullptr) //Minuit2 not always installed!
	minimum = ROOT::Math::Factory::CreateMinimizer("Minuit", "");
      
      // set tolerance , etc...
      minimum->SetMaxFunctionCalls(1000000); // for Minuit/Minuit2
      minimum->SetMaxIterations(1000);  // for GSL
      minimum->SetTolerance(0.0001);
      minimum->SetPrintLevel(0);
      
      // create function wrapper for minimizer
      // a IMultiGenFunction type
      ROOT::Math::Functor wrapf(Fmax,2);

      //variable = W, variable 1 = X, variable 2 = Y

      double step[2] = {Wrange/100,1./100,1.100};
      // starting point
      
      // double variable[2] = { gridW,gridt};
      double variable[2] = {Wmin+Wrange/2,0.5,0.5};
      
      minimum->SetFunction(wrapf);
      
      // Set the free variables to be minimized !
      minimum->SetVariable(0,"W",variable[0], step[0]);
      minimum->SetVariable(1,"X",variable[1], step[1]);
      minimum->SetVariable(2,"Y",variable[2], step[2]);
      // minimum->SetVariableLimits(0,Wmin,_Wmax);
      // minimum->SetVariableLimits(1,tmax,0);

      // do the minimization
      minimum->Minimize();
      const double *xs = minimum->X();
 
      auto minVal = minimum->MinValue();
      auto minW= xs[0];
      auto mint= kine::TFromXY(xs[1],xs[2]);
      auto minM2= kine::M2FromXY(xs[1],xs[2]);
   
      std::cout << "Maximum : Probabiltiy Dist at ( W=" << minW << " , t = "  <<<< " , M2 = "  << minM2<< " , x = "  <<xs[1]<< " , y = "  <<xs[2] << "): "<< -minimum->MinValue()  << std::endl;

      //check for low mass meson limits
      if(dynamic_cast<DecayingParticle*>(_meson)){ //meson
	dynamic_cast<DecayingParticle*>(_meson)->TakeMinimumMass();//to get threshold behaviour
	M3=_meson->Mass();
	
	// do the minimization at min mass in case higher max
	minimum->Minimize();
	const double *xs = minimum->X();
	
	auto minminVal = minimum->MinValue();
	

	if(minminVal<minVal){
	  auto minW= xs[0];
	  auto mint= kine::TFromXY(xs[1],xs[2]);
	  auto minM2= kine::M2FromXY(xs[1],xs[2]);
	  minVal=minminVal;
	  std::cout<<"minmin "<<-minminVal<<" < "<<-minVal<<std::endl;

	  std::cout << "Minumum Maximum : Probabiltiy Dist at ( W=" << minW << " , t = "  <<<< " , M2 = "  << minM2<< " , x = "  <<xs[1]<< " , y = "  <<xs[2] << "): "<< -minimum->MinValue()  << std::endl;
	}
	//back to PDg mass if exists
	if(_meson->PdgMass()>M3){
	  M3=_meson->PdgMass();
	}  
	dynamic_cast<DecayingParticle*>(_meson)->TakePdgMass();
	
	
      }

      if(gridMin<minVal){
	Warning("DecayModelstM2::FindMaxOfIntensity()","grid search value already bigger than minimised, so will revert to that max value +5 percent");
	std::cout<<"gridMin "<<gridMin<<" "<<minVal<<std::endl;
	minVal=gridMin*1.05;
      }

      return -minVal;
  }
 

  void DecayModelstM2::HistMaxXSection(TH1D& hist){

 
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
	else if( TMath::IsNaN(kine::tmax(_W,M1,M2,M3,M4)) )
	  hist.SetBinContent(ih, 0);
	else if( TMath::IsNaN(kine::t0(_W,M1,M2,M3,M4)) )
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
