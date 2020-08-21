#include "JpacModelst.h"
#include "DecayingParticle.h"
#include "FunctionsForJpac.h"
#include "FunctionsForGenvector.h"
#include <TDatabasePDG.h>

namespace elSpectro{
  ///////////////////////////////////////////////////////
  ///constructor includes subseqent decay of Ngamma* system
  JpacModelst::JpacModelst( jpacPhoto::amplitude* amp ,
			      particle_ptrs parts, const std::vector<int> pdgs) :
    _amp{amp},
    DecayModel{ parts, pdgs }
  {
    _name={"JpacModelst"};

    //need to find meson and baryon
    if(TDatabasePDG::Instance()->GetParticle(Products()[0]->Pdg())->ParticleClass()==TString("Baryon") ){
      _baryon=Products()[0];
      _meson=Products()[1]; 
    }
    else {
      _baryon=Products()[1];
      _meson=Products()[0];
    }
    
    //find max value of amplitude as function of s and t
    _max = 2E5;//_jpac->MaxValue();
    // _max = 1E8;//_jpac->MaxValue();
    std::cout<<"JpacModelst::JpacModelst "<<_amp<<std::endl;
  }
  /////////////////////////////////////////////////////////////////
  void JpacModelst::PostInit(ReactionInfo* info){

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

     _max = jpacFun::FindMaxOfProbabilityDistribution(_amp,maxW);

     std::cout<<"JpacModelQ2W::PostInit max value "<<_max<<" "<<_meson<<" "<<_meson->Pdg()<<" "<<_sdmeMeson<<std::endl;
  }
  //////////////////////////////////////////////////////////////////
  double JpacModelst::Intensity() const
  {
  
   
    auto parent = ParentVector();
    double s = parent.M();
    s*=s;

  
    //For t, must rotate into the parent g*N z-axis frame
    /*
    auto cmBoost=parent.BoostToCM();
    ROOT::Math::RotationY rotateToZaxis;
    rotateToZaxis.SetAngle(parent.Theta());

    auto cmMeson = _meson->P4();
    cmMeson=ROOT::Math::VectorUtil::boost(cmMeson,cmBoost);
    cmMeson=rotateToZaxis*cmMeson;
     */
    double t = (_meson->P4()-*_photon).M2();//_amp->kinematics->t_man(s,cmMeson.Theta());

    MomentumVector decayAngles;
    kine::electroCMDecay(&parent,_ebeam,_photon,_meson->P4ptr(),&decayAngles);
    _photonPol->SetPhi(decayAngles.Phi());
    
     //jpac photo amp depends on s, t, and meson mass
    _amp->kinematics->set_vectormass( _meson->Mass() );
    double weight = _amp->differential_xsection(s,t)/_max;

    
    //Meson spin density marix elements, note this is photoproduced
    if(_sdmeMeson){
      _sdmeMeson->SetElement(0,0,0,(_amp->SDME(0, 0, 0, s, t)));
      _sdmeMeson->SetElement(0,1,0,(_amp->SDME(0, 1, 0, s, t)));
      _sdmeMeson->SetElement(0,1,-1,(_amp->SDME(0, 1, -1, s, t)));
      _sdmeMeson->SetElement(1,1,1,(_amp->SDME(1, 1, 1, s, t)));
      _sdmeMeson->SetElement(1,0,0,(_amp->SDME(1, 0, 0, s, t)));
      _sdmeMeson->SetElement(1,1,0,(_amp->SDME(1, 1, 0, s, t)));
      _sdmeMeson->SetElement(1,1,-1,(_amp->SDME(1, 1, -1, s, t)));
      _sdmeMeson->SetElement(2,1,0,(_amp->SDME(2, 1, 0, s, t)));
      _sdmeMeson->SetElement(2,1,-1,(_amp->SDME(2, 1, -1, s, t)));
      //_sdmeMeson->SetElement(3,1,0,(_amp->SDME(3, 1, 0, s, t)));
      // _sdmeMeson->SetElement(3,1,-1,(_amp->SDME(3, 1, -1, s, t)));
    }

    if(weight>1){
      auto oldmax=_max;
      _max=weight*oldmax;
      std::cout<<"JpacModelst::Intensity changing max from  "<<oldmax<<" to "<<_max<<std::endl;
    }
    
    //Correct for Q2&W weighting which has already been applied
    weight/=_prodInfo->_sWeight;
    
    if(weight>1){
      std::cout<<"JpacModelst::Intensity sWeight corrected weight too large "<<weight <<" "<<_prodInfo->_sWeight<<std::endl;
    }
    return weight;
  }
  
}
