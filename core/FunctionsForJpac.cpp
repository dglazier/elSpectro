#include "FunctionsForJpac.h"

#include <Math/GSLIntegrator.h>
#include <Math/IntegrationTypes.h>
#include <Math/Functor.h>
#include <Math/Minimizer.h>
#include <Math/Factory.h>
#include <TRandom3.h>

namespace elSpectro {

  namespace jpacFun{

    void HistProbabilityDistribution_s(jpacPhoto::amplitude* amp, TH1D&  hist){
      double s=0;
      auto F = [amp,&s](double t)
	{
	  return amp->differential_xsection(s, t);
	};

      ROOT::Math::GSLIntegrator ig(ROOT::Math::IntegrationOneDim::kADAPTIVE,
				   ROOT::Math::Integration::kGAUSS61);
      ROOT::Math::Functor1D wF(F);
      ig.SetFunction(wF);

  
      for(int ih=1;ih<=hist.GetNbinsX();ih++){
	double Wval=hist.GetXaxis()->GetBinCenter(ih);
	s=Wval*Wval;
	if( s < amp->kinematics->sth )
	  hist.SetBinContent(ih, 0);
	else
	  hist.SetBinContent(ih,ig.Integral(amp->kinematics->t_man(s,TMath::Pi()),amp->kinematics->t_man(s,0)));
	if(ih%10==0)std::cout<<(hist.GetNbinsX() - ih)/10<<" ";
      }
      std::cout<<std::endl;
      //done
    }

    //////////////////////////////////////////////////////////////////
    double FindMaxOfProbabilityDistribution(jpacPhoto::amplitude* amp,double Wmax){
      auto Fmax = [amp](const double *x)
	{
	  double s = x[0]*x[0];
	  if(s < amp->kinematics->sth) return 0.;
	  if(x[1]> amp->kinematics->t_man(s,0)) return 0.;
	  if(x[1]< amp->kinematics->t_man(s,TMath::Pi())) return 0.;
	  
	  double val=amp->differential_xsection(s , x[1]);
	  
	  if( TMath::IsNaN(val) ) return 0.;
	  return -val; //using a minimiser!
	};
      
      //First perform grid search for intital values
      double Wmin=TMath::Sqrt(amp->kinematics->sth);
      double Wrange=Wmax-Wmin;
      double tmax=amp->kinematics->t_man(Wmax*Wmax,TMath::Pi());
      
      double gridMin=0;
      double gridW=0;
      double gridt=0;
      double WtVals[2];
      int Npoints=50;
      for(int iW=1;iW<Npoints;iW++){
	WtVals[0]=Wmin+iW*Wrange/Npoints;
	double tmin=amp->kinematics->t_man(WtVals[0]*WtVals[0],0);
	double tmaxg=amp->kinematics->t_man(WtVals[0]*WtVals[0],TMath::Pi());
	double trange= tmin-tmaxg;

	for(int it=0;it<Npoints;it++){
	  WtVals[1]=tmin-it*trange/Npoints;
	  auto val = Fmax(WtVals);
	  // std::cout<<WtVals[0]<<" "<<WtVals[1]<<" "<<val<<std::endl;
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
      minimum->SetVariableLimits(0,Wmin,Wmax);
      minimum->SetVariableLimits(1,tmax,0);
      
      // do the minimization
      minimum->Minimize();
     const double *xs = minimum->X();
 
      auto minVal = minimum->MinValue();
      auto minW= xs[0];
      auto mint= xs[1];
      /*
      for(int i=0;i<4;i++){
	auto valW=gRandom->Uniform(TMath::Sqrt(amp->kinematics->sth),Wmax);
	auto valt=gRandom->Uniform(tmax,0);
	std::cout<<"Starting values "<<valW<<" "<<valt<<std::endl;
	minimum->SetVariableValue(0,valW);
	minimum->SetVariableValue(1,valt);
	// do the minimization
	minimum->Minimize();
	if(minVal>minimum->MinValue()){
	  minVal=minimum->MinValue();
	  minW= xs[0];
	  mint= xs[1];
	}
      }
      */
      std::cout << "Maximum : Probabiltiy Dist at ( W=" << minW << " , t = "  << mint << "): "<< -minimum->MinValue()  << std::endl;
      return -minVal ;
    }
    /////////////////////////////////////////////////
  }
}
