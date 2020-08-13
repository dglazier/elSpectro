#include "JpacModelst.h"
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
    DecayModel::PostInit(info);
    
    _prodInfo= dynamic_cast<ReactionElectroProd*> (info); //I need Reaction info
 
     double maxW = ( *(_prodInfo->_target) + *(_prodInfo->_ebeam) ).M();

     _max = jpacFun::FindMaxOfProbabilityDistribution(_amp,maxW);
    std::cout<<"JpacModelQ2W::PostInit max value "<<_max<<std::endl;
  }
  //////////////////////////////////////////////////////////////////
  double JpacModelst::Intensity() const
  {
  
   
    auto parent = ParentVector();
    double s = parent.M();
    s*=s;

  
    //must rotate into the g*N z-axis frame
    auto cmBoost=parent.BoostToCM();
    auto rfm = _meson->P4();
    ROOT::Math::RotationY rotateToZaxis;
    rotateToZaxis.SetAngle(parent.Theta());

    rfm=ROOT::Math::VectorUtil::boost(rfm,cmBoost);
    rfm=rotateToZaxis*rfm;
    double t = _amp->kinematics->t_man(s,rfm.Theta());

    //jpac photo amp depends on s, t, and meson mass
    _amp->kinematics->set_vectormass( _meson->Mass() );
    double weight = _amp->differential_xsection(s,t)/_max;
   
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
