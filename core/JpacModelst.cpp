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
    DecayModelst{ parts, pdgs }
  {
    _name={"JpacModelst"};

 
    std::cout<<"JpacModelst::JpacModelst "<<_amp<<std::endl;
  }
  /////////////////////////////////////////////////////////////////
  /*void JpacModelst::PostInit(ReactionInfo* info){

    DecayModelst::PostInit(info);
    
    //double maxW = ( *(ProductionInfo()->_target) + *(ProductionInfo()->_ebeam) ).M();

    // _max = jpacFun::FindMaxOfProbabilityDistribution(_amp,maxW)*1.2; //add 10% for Q2 effects etc.

    std::cout<<"JpacModelst::PostInit max value "<<" "<<_meson<<" "<<_meson->Pdg()<<" "<<_sdmeMeson<<std::endl;
    }*/
  //////////////////////////////////////////////////////////////////
  void JpacModelst::CalcMesonSDMEs() const {
    //Meson spin density marix elements, note this is photoproduced
    auto *sdme=GetMesonSDMEs();
    
    if(sdme){ //note this is vector formalism
      if(sdme->Spin()==1){
	sdme->SetElement(0,0,0,(_amp->SDME(0, 0, 0, get_s(), get_t())));
	sdme->SetElement(0,1,0,(_amp->SDME(0, 1, 0, get_s(), get_t())));
	sdme->SetElement(0,1,-1,(_amp->SDME(0, 1, -1, get_s(), get_t())));
	sdme->SetElement(1,1,1,(_amp->SDME(1, 1, 1, get_s(), get_t())));
	sdme->SetElement(1,0,0,(_amp->SDME(1, 0, 0, get_s(), get_t())));
	sdme->SetElement(1,1,0,(_amp->SDME(1, 1, 0, get_s(), get_t())));
	sdme->SetElement(1,1,-1,(_amp->SDME(1, 1, -1, get_s(), get_t())));
	sdme->SetElement(2,1,0,(_amp->SDME(2, 1, 0, get_s(), get_t())));
	sdme->SetElement(2,1,-1,(_amp->SDME(2, 1, -1, get_s(), get_t())));
      }
      else if(sdme->Spin()==2){
	sdme->SetElement(0,0,0,(_amp->SDME(0, 0, 0, get_s(), get_t())));
	sdme->SetElement(0,1,0,(_amp->SDME(0, 1, 0, get_s(), get_t())));
	sdme->SetElement(0,1,-1,(_amp->SDME(0, 1, -1, get_s(), get_t())));
	sdme->SetElement(0,1,1,(_amp->SDME(0, 1, 1, get_s(), get_t())));
	sdme->SetElement(0,2,-2,(_amp->SDME(0, 2, -2, get_s(), get_t())));
	sdme->SetElement(0,2,0,(_amp->SDME(0, 2, 0, get_s(), get_t())));
	sdme->SetElement(0,2,1,(_amp->SDME(0, 2, 1, get_s(), get_t())));
	sdme->SetElement(0,2,2,(_amp->SDME(0, 2, 2, get_s(), get_t())));
	
	sdme->SetElement(1,0,0,(_amp->SDME(0, 0, 0, get_s(), get_t())));
	sdme->SetElement(1,1,0,(_amp->SDME(0, 1, 0, get_s(), get_t())));
	sdme->SetElement(1,1,-1,(_amp->SDME(0, 1, -1, get_s(), get_t())));
	sdme->SetElement(1,1,1,(_amp->SDME(0, 1, 1, get_s(), get_t())));
	sdme->SetElement(1,2,-1,(_amp->SDME(0, 2, -1, get_s(), get_t())));
	sdme->SetElement(1,2,-2,(_amp->SDME(0, 2, -2, get_s(), get_t())));
	sdme->SetElement(1,2,0,(_amp->SDME(0, 2, 0, get_s(), get_t())));
	sdme->SetElement(1,2,1,(_amp->SDME(0, 2, 1, get_s(), get_t())));
	sdme->SetElement(1,2,2,(_amp->SDME(0, 2, 2, get_s(), get_t())));

	sdme->SetElement(2,1,0,(_amp->SDME(0, 1, 0, get_s(), get_t())));
	sdme->SetElement(2,1,-1,(_amp->SDME(0, 1, -1, get_s(), get_t())));
	sdme->SetElement(2,2,-1,(_amp->SDME(0, 2, -1, get_s(), get_t())));
	sdme->SetElement(2,2,-2,(_amp->SDME(0, 2, -2, get_s(), get_t())));
	sdme->SetElement(2,2,0,(_amp->SDME(0, 2, 0, get_s(), get_t())));
	sdme->SetElement(2,2,1,(_amp->SDME(0, 2, 1, get_s(), get_t())));



      }
    }


  }
  //////////////////////////////////////////////////////////////////
  void JpacModelst::CalcBaryonSDMEs() const {
    auto  *sdme=GetBaryonSDMEs();
    
    if(sdme){

    }
  }
    /*double JpacModelst::Intensity() const
  {
  
   
    double s=Parent->P4().M2();
  
    double t = (_meson->P4()-*_photon).M2();//_amp->kinematics->t_man(s,cmMeson.Theta());

    _amp->kinematics->set_vectormass( _meson->Mass() );


    double weight = _amp->probability_distribution(s,t)*PhaseSpaceFactor();
    
      

    
    // double weight2 = _amp->differential_xsection(s,t);//* (2.56819E-6)*4*64*TMath::Pi() * s * (TMath::Power(std::real(_amp->kinematics->initial->momentum(s)), 2.));
    //double weight2 = _amp->differential_xsection(s,t);
    // std::cout<<"Compare "<<weight2<<" "<<weight<<" "<<weight2/weight <<std::endl;

    weight/=_max; //normalise range 0-1
    
    //  std::cout<<"flux "<<kine::FluxPhaseSpaceFactor(*_photon,*get_t()arget)<<" "<<kine::PDK(W,_meson->Mass(),_baryon->Mass())*W<<std::endl;
    //std::cout<<"dt "<<kine::FluxPhaseSpaceFactor(*_photon,*get_t()arget)*kine::PDK(W,_meson->Mass(),_baryon->Mass())/W<<" "<< kine::PhaseSpaceFactorDt(W,_photon->M(),get_t()arget->M(),_meson->Mass(),_baryon->Mass())<<std::endl;
    
    //Meson spin density marix elements, note this is photoproduced
    if(get_s()dmeMeson){
      get_s()dmeMeson->SetElement(0,0,0,(_amp->SDME(0, 0, 0, s, t)));
      get_s()dmeMeson->SetElement(0,1,0,(_amp->SDME(0, 1, 0, s, t)));
      get_s()dmeMeson->SetElement(0,1,-1,(_amp->SDME(0, 1, -1, s, t)));
      get_s()dmeMeson->SetElement(1,1,1,(_amp->SDME(1, 1, 1, s, t)));
      get_s()dmeMeson->SetElement(1,0,0,(_amp->SDME(1, 0, 0, s, t)));
      get_s()dmeMeson->SetElement(1,1,0,(_amp->SDME(1, 1, 0, s, t)));
      get_s()dmeMeson->SetElement(1,1,-1,(_amp->SDME(1, 1, -1, s, t)));
      get_s()dmeMeson->SetElement(2,1,0,(_amp->SDME(2, 1, 0, s, t)));
      get_s()dmeMeson->SetElement(2,1,-1,(_amp->SDME(2, 1, -1, s, t)));
      //get_s()dmeMeson->SetElement(3,1,0,(_amp->SDME(3, 1, 0, s, t)));
      // get_s()dmeMeson->SetElement(3,1,-1,(_amp->SDME(3, 1, -1, s, t)));
    }

    if(weight>1){
      //don't change weight, likely due to large Q2 value....
      
      //auto oldmax=_max;
      //_max=weight*oldmax;
      std::cout<<"JpacModelst::Intensity weight too high but won't change max (prob Q2 too high) from  "<<_max<<" to "<<weight*_max<<std::endl;
    }
    
    //Correct for Q2&W weighting which has already been applied
    weight/=_prodInfo->get_s()Weight;
    
    if(weight>1){
      std::cout<<" s weight "<<_prodInfo->get_s()Weight<<" Q2 "<<-_photon->M2()<<" 2Mmu "<<2*get_t()arget->M()*_photon->E() <<" W "<<W<<" t "<<t<<" new weight "<<weight*_prodInfo->get_s()Weight<<" meson "<<_meson->Mass()<<std::endl;
      std::cout<<"JpacModelst::Intensity sWeight corrected weight too large "<<weight <<" "<<_prodInfo->get_s()Weight<<" diff Xsect "<<_amp->differential_xsection(s,t)<<std::endl;
																				    std::cout<<"flux "<<kine::FluxPhaseSpaceFactor(*_photon,*get_t()arget)<<" "<<4*kine::PDK(W,_meson->Mass(),_baryon->Mass())*W<<" "<< kine::PhaseSpaceFactorDt(W,P1,_meson->Mass(),_baryon->Mass())<<std::endl;				    
      std::cout<<"Pi check "<<P1 <<" versus "<< kine::PDK(W,_photon->M(),get_t()arget->M())<<std::endl;
    }
     
 
    return weight;
  }
  */
}
