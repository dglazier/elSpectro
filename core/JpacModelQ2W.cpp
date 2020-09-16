#include "JpacModelQ2W.h"
#include "TwoBody_stu.h"
#include "FunctionsForJpac.h"
#include "Interface.h" //for generator

namespace elSpectro{
  ///////////////////////////////////////////////////////
  ///constructor includes subseqent decay of g*N system
  ///this class controls decay of e+p ->e' + g*N
  JpacModelQ2W::JpacModelQ2W( jpacPhoto::amplitude* amp ,
			      particle_ptrs parts, const std::vector<int> pdgs):
    _amp{amp},
    DecayModelQ2W{0, new JpacModelst{amp, parts, pdgs},
	new TwoBody_stu{0,1,0.5,0,0} }
    // ->   DecayModel{{ new DecayingParticle{-2211,Decay_st} },{11}}
  {
    _name={"JpacModelQ2W"};

   
  }
  ///////////////////////////////////////////////////////
  ///constructor includes subseqent decay of Ngamma* system with s and t channel
  JpacModelQ2W::JpacModelQ2W( jpacPhoto::amplitude* amp ,
			      particle_ptrs parts, const std::vector<int> pdgs,
			      double s_strength,
			      double t_strength,double t_slope):
    _amp{amp},
    DecayModelQ2W{0, new JpacModelst{amp, parts, pdgs},
	new TwoBody_stu{s_strength,t_strength, t_slope,0,0} }
    // ->   DecayModel{{ new DecayingParticle{-2211,Decay_st} },{11}}
  {
    _name={"JpacModelQ2W_given_s_and_t"};
  }

  ////////////////////////////////////////////////////////
  void JpacModelQ2W::PostInit(ReactionInfo* info){

    DecayModelQ2W::PostInit(info);

    //declare the model for controlling phase space
    //this is the decay of the overall hadronic system
    // move to ElectronScatter  generator().SetModelForMassPhaseSpace(GetGammaN()->Model());
      generator().SetModelForMassPhaseSpace(GetGammaN()->Model());
    
    //JPAC model is for photoproduction
    //Leave Q2 dependence for now
    //W dependence requires integration over W
    //Do we need this or just evaluate at each s and t via JpacModelst
    auto prodInfo=ProdInfo();
    double maxW = ( *(prodInfo->_target) + *(prodInfo->_ebeam) ).M();

    std::cout<<"JpacModelQ2W::PostInit generating total cross section, may take some time... "<<std::endl;
    TH1D histlow("Wdistlow","Wdistlow",50,0,maxW);
    _amp->kinematics->set_vectormass(GetDecayMeson()->MinimumMassPossible());//to get threshold behaviour
    jpacFun::HistProbabilityDistribution_s(_amp,histlow);

    TH1D histpeak("Wdist","Wdist",50,0,maxW);
    _amp->kinematics->set_vectormass(GetDecayMeson()->PdgMass());//to get threshold behaviour
    jpacFun::HistProbabilityDistribution_s(_amp,histpeak);
    
    auto hist = jpacFun::HistFromLargestBins(histlow,histpeak);
    
    //    auto maxVal= hist.GetMaximum();
   
    //for(int ibin=1;ibin<=hist.GetNbinsX();ibin++)
    // hist.SetBinContent(ibin,hist.GetBinContent(ibin)+0.01*maxVal);
    
    _W_Dist.reset( new DistTH1(hist) );
    
  
  }
  double JpacModelQ2W::Intensity() const{

    if(CheckThreshold()==false){
      return 0.;
    }
    if(GetGammaN()->P4().M()<GetGammaN()->MinimumMassPossible()) return 0;

    auto phasespace_weight=DecayModelQ2W::Intensity(); //does virtual photon calculations
   
  
    if(_W_Dist.get()==nullptr) return 1.;

    auto prodInfo=ProdInfo(); 
 
    double weight=_W_Dist->GetWeightFor( prodInfo->_photoN->M() );
    
    prodInfo->_sWeight=weight; //might be used in s and t
    
    return weight;
  }
  
  
}
