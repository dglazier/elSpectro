#include "DecayModelDnpee.h"
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
  ///fix decay products
  DecayModelDnpee::DecayModelDnpee() :
    //_proton(new Particle{2212}),_neutron(new Particle{2112}),_ele(new Particle{11}),_pos(new Particle{-11}),
    // DecayModel( {_proton,_neutron,_ele,_pos},{} )
    DecayModel( {},{2212,2112,11,-11} )
  {
    _name={"DecayModelDnpee"};

    if(Products().size()!=4){
      Fatal("DecayModelDnpee","Can only have four decay particles");
    }

    auto products = StableProducts();
    for(auto* prod:products){
      if(prod->Pdg()==2212) _proton=prod;
      if(prod->Pdg()==2112) _neutron=prod;
      if(prod->Pdg()==11) _ele=prod;
      if(prod->Pdg()==-11) _pos=prod;
    }
  }
  
   void DecayModelDnpee::SetParent(DecayingParticle* pa){
    DecayModel::SetParent(pa);
    //now can make cascade of decays
    _phaseSpace.SetParentAndProducts(Parent(),{_proton,_neutron,_ele,_pos},{});
    //  std::cout<<"DecayModelDnpee::SetParent  ResetProducts "<<_phaseSpace.Products().size()<<" "<<_phaseSpace.Product(0)->Pdg()<<" "<<_phaseSpace.Product(0)->MinimumMassPossible()<<" "<<_phaseSpace.Product(1)->Pdg()<<std::endl;
    ResetProducts(_phaseSpace.Products());
  }
  /////////////////////////////////////////////////////////////////
  void DecayModelDnpee::PostInit(ReactionInfo* info){
    
    DecayModel::PostInit(info);
    _phaseSpace.PostInit(info);
   
    _isElProd=kTRUE;
    _prodInfo= dynamic_cast<ReactionElectroProd*> (info); //I need Reaction info
    if(_prodInfo==nullptr){
      _isElProd=kFALSE;
      _prodInfo= dynamic_cast<ReactionPhotoProd*> (info);
    }
    _photon = _prodInfo->_photon;
    _target = _prodInfo->_target;
    _ebeam = _prodInfo->_ebeam;
    
    // double maxW = ( *(_prodInfo->_target) + *(_prodInfo->_ebeam) ).M();
    _Wmax = _prodInfo->_Wmax;
 
     _max = FindMaxOfIntensity()*1.08; //add 5% for Q2,meson mass effects etc.

     std::cout<<"DecayModelDnpee::PostInit max value "<<_max<<std::endl;
  }

  double DecayModelDnpee::Intensity() const
  {
    //First you must calculate all the variables that the cross section depends on
    //add them as data members and use them in DifferentialXSect()
    _W = Parent()->P4().M();
    _s=_W*_W;
 
    // _dvolume = _dt; //NEED TO DEFINE PHASE SPACE RANGE AT CURRENT VALUE OF S. DEPENDS ON VARIABLES etc
    _dvolume=(_W-_Wmin)*(_W-_Wmin);//THSIS IS NOT RIGHT, JUST A PROXY!!!

    double weight = DifferentialXSect() * _dvolume ;//must multiply by ranges for correct sampling
    //std::cout<<"********************DecayModelDnpee::Intensity() "<<_W<<" weight = "<<weight<<" max = "<<_max<<std::endl;
    weight/=_max; //normalise range 0-1

    if(_isElProd==kTRUE)
      weight/= TMath::Sqrt(PgammaCMsq()/kine::PDK2(_W,0,_target->M())); //correct max for finite Q2 phase space

    if(weight>1){
      //don't change weight, likely due to large Q2 value....
      std::cout<<"DecayModelDnpee::Intensity weight too high but won't change maxprobable low meson mass and W from  "<<_max<<" to "<<weight*_max<<" W "<<_W<<std::endl;
    }

   //Correct for W weighting which has already been applied
    weight/=_prodInfo->_sWeight;
    // std::cout<<"********************DecayModelDnpee::Intensity() "<<_W<<" "<<weight<<" "<<_prodInfo->_sWeight<<std::endl;
    return weight;
    

  }

    double DecayModelDnpee::FindMaxOfIntensity(){

      _Wmax = _prodInfo->_Wmax;
      _Wmin = Parent()->MinimumMassPossible();
      double Wrange=_Wmax-_Wmin;
      std::cout<<"Search "<<_Wmin<<" to "<<_Wmax<<std::endl;
      auto Fmax = [this](const double *x)
		  {
		    _s = x[0]*x[0];
		    _W=x[0];
		    
  //_dvolume = _dt; //NEED TO DEFINE PHASE SPACE RANGE AT CURRENT VALUE OF S. DEPENDS ON VARIABLES etc
  		    _dvolume=(_W-_Wmin)*(_W-_Wmin);//THSIS IS NOT RIGHT, JUST A PROXY!!!
		    double val = DifferentialXSect()*_dvolume;
		    // std::cout<<"DecayModelDnpee::FindMaxOfIntensity() "<<val<<" volume "<<_dvolume<<" at W = "<<_W<<std::endl;
		    if( TMath::IsNaN(val) ) return 0.;
		    return -(val); //using a minimiser!
		  };
      
      
      
      ROOT::Math::Minimizer* minimum =
	ROOT::Math::Factory::CreateMinimizer("Genetic", "");
  
      // set tolerance , etc...
      minimum->SetMaxFunctionCalls(1000000); // for Minuit/Minuit2
      minimum->SetMaxIterations(1000);  // for GSL
      minimum->SetTolerance(0.0001);
      minimum->SetPrintLevel(0);
      
      // create function wrapper for minimizer
      // a IMultiGenFunction type
      ROOT::Math::Functor wrapf(Fmax,1); //currently only dependent on W
      double step[1] = {Wrange/100};
      double variable[1] = {_Wmin+Wrange/2};
      
      minimum->SetFunction(wrapf);
      
      // Set the free variables to be minimized !
      minimum->SetVariable(0,"W",variable[0], step[0]);

      // do the minimization
      minimum->Minimize();
      const double *xs = minimum->X();
 
      auto minVal = -minimum->MinValue();
      auto minW= xs[0];
      std::cout << "Maximum : Probabiltiy Dist at ( W=" << minW <<  "): "<< minVal  <<  std::endl;
      return minVal;
    }
}
