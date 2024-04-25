//////////////////////////////////////////////////////////////
///
///Class:		JpacSemiInc
///Description:
///             Control behaviour of Particle decay to Particle products
///             Defined by
///             1) preconfigured jpacPhoto amplitude
///             2) it decay as a function of s and t JpacDecayst
///
///            Note derived classes should include a constructor to initialise
///            JpacSemiInc( particle_ptrs , const std::vector<int> pdgs );
/// Reference https://journals.aps.org/prd/abstract/10.1103/PhysRevD.106.094009
#pragma once

#include "TwoBodyProduction.h"
#include "FunctionsForElectronScattering.h"
#include "FunctionsForKinematics.h"
#include "inclusive/inclusive_production.hpp"
#include "DistTH1.h"

namespace elSpectro{

  using jpacInc_ptr = jpacPhoto::inclusive_production*;

  class JpacSemiInc : public TwoBodyProduction {

  public:
    
    JpacSemiInc()=delete;
    //constructor giving jpac amplitude pointer (which we will now own)
    //and decay particles 
    JpacSemiInc( jpacInc_ptr inc, particle_ptrs parts,
		  const std::vector<int> pdgs  );

    constexpr double sixteen_pi_3() const{return 16.*TMath::Pi()*TMath::Pi()*TMath::Pi();}
    
   void PostInit(ReactionInfo* info) override;

    void CreateMassEnvelope();
    const TH1D* GetM2Distribution() {
      if( _createMassEnvelope==false)
	return nullptr;
      return &_massEnvelope->GetTH1();
    }
    
    void HistMaxXSection(TH1D& hist) final;
    double FindMaxOfIntensity() final;//{return get_max(); }
    double dsigma_dcosth() const final;
    double dsigma_dcosth(double W,double costheta) const final;
    double sigma(double W) const final{
      set_W(W);
      auto xs=  _inc->integrated_xsection(get_s());
    std::cout<<"JpacSemiInc::sigma @ s = "<<get_s() <<" = "<< xs <<" photon mass  = "<<GetPhoton()->M()<<"\n";
      return xs;
    }
    
    bool HasAngularDistribution() override{return true; } //I have an angular distribution
 
    double PhaseSpaceFactor() const final{
      //with jacobian for d2sigma/drdcth, now use PhaseSpaceFactoratQ2!=0
      //= Jacobian*K/16pi
      auto result =   _jacobian*sqrt(kine::Kallen(get_MB2(), get_t(), escat::M2_pr()) ) / (2*get_W()*PgammaCM() ) / sixteen_pi_3();
      //But mass phase space of B->pi nucleon already done...
      // multiply by Pbreakup = PDK(get_MB(),Mpi,Mnucleon)
      //  auto baryonPhaseSpace = TMath::Sqrt(dynamic_cast<const DecayingParticle*>(GetBaryon())->Model()->PhaseSpaceWeightSq(TMath::Sqrt(get_MB2())));
      double baryonPhaseSpace = 1.0;
      if(IsSampling()) baryonPhaseSpace = TMath::Sqrt(dynamic_cast<const DecayingParticle*>(GetBaryon())->Model()->PhaseSpaceWeightSq(TMath::Sqrt(get_MB2())));
      // std::cout<<"JpacSemiInc PhaseSpaceFactor() "<<TMath::Sqrt(get_MB2())<<" "<<baryonPhaseSpace<<std::endl;
      return result/baryonPhaseSpace;
    }
    // double PhaseSpaceWeightSq(double W) final{
    //   double result = 0;
    //   std::cout<<"PhaseSpaceWeightSq JpacSemiInc  " << (result = kine::PDK2(W,Product(0)->Mass(),Product(1)->Mass()))  << std::endl;
    //   for(auto* p:UnstableProducts()){
    // 	std::cout<<"PhaseSpaceWeightSq JpacSemiInc  " << p->Pdg()<<" "<<p->Mass()<<" "<<p->PhaseSpaceWeightSq()<<std::endl;
    // 	result*=p->PhaseSpaceWeightSq();
    //   }
    //   return result;
    // }
    double MatrixElementsSquared_T() const override {
      // if(get_W()<_inc->_kinematics->Wth()) return 0;
     //make sure mass is above threshold, as meson and baryon masses may have changed
      if( get_W() < (TMath::Sqrt(get_MB2())+GetMeson()->P4().M()) ) return 0;
       // Pass the total energy to the kinematics object
      _inc->_kinematics->_s = get_s();
      _inc->_kinematics->set_meson_mass( GetMeson()->Mass() );
      // std::cout<<"JpacSemiInc::MatrixElementsSquared_T meson mass "<<GetMeson()->P4().M()<<" "<< GetMeson()->Mass()<<" "<<_inc->_kinematics->get_meson_mass()<<std::endl;
      //std::cout<<"JpacSemiInc::MatrixElementsSquared_T W = "<<get_W()<<" "<<TMath::Sqrt(get_MB2())+GetMeson()->P4().M()<<" tmin "<< _inc->_kinematics->TMINfromM2(get_MB2() )<<std::endl;

      //check mb
      auto incXS=_inc->d3sigma_d3p(get_s(), get_t(),get_MB2());
      // std::cout<<"JpacSemiInc::MatrixElementsSquared_T check baryon mass "<<GetBaryon()->M2()<<" =? "<<M2FromR(get_r())<<" and r = "<<get_r()<<" and getMB2 "<<get_MB2()<<" ps "<<PhaseSpaceFactoratQ20()<<" inc "<<incXS<<std::endl;
      if(get_MB2()==0.) { std::cout<<"JpacSemiInc::MatrixElementsSquared_T zero baryon mass "<<std::endl;exit(0);}
      //convert cross section back to matrix element
      return incXS/PhaseSpaceFactoratQ20()*sixteen_pi_3();
      //return _inc->dsigma_drdcth(get_s(),get_r(),get_cosThCM());
 
    }
    double get_r() const{ return _r;}
    double get_MB2() const{ return _mB2;}
    
    void set_r(double val) {_r=val;}
    void set_MB2(double val) {_mB2=val;}

    void SetCreateM2Distribution(){_createMassEnvelope=true;}

private:
    
    double DiffXSFromRCosThatQ20() const{
      JacobianRCosTh();
      set_t( TFromWRCosThatQ20(get_W(),get_r(),get_cosThCM()));

      //std::cout<<"DiffXSFromRCosThatQ20() "<<_jacobian<<" "<<get_t()<<std::endl;
      return  DiffXS();
    }
    
    double M2FromR(double r) const
    {
      auto mX2 =GetMeson()->M2();
      auto qmax2 = kine::PDK2(get_W(), GetMeson()->Mass(),_minMassB);
      //eqn(A13) Mb^2= Mm^2 + s - 2WsEm
      return mX2 + get_s() - 2. * get_W()*sqrt(mX2  + qmax2  * r * r );
    };
    double RFromM2(double M2) const{
      //r^2 = (Em^2 - Mm^2)/qmax^2 = q2/qmax^2 ;  qmax2 =PDK2(W,Mm,Mb_min); q2 = PDK2(W,Mm,Mb); Em^2=Mm^2+q2
      auto Mm = GetMeson()->Mass();
      auto qmax2 = kine::PDK2(get_W(), Mm,_minMassB);
      auto q2  =  kine::PDK2(get_W(), Mm,TMath::Sqrt(M2));
      auto r2 = q2/qmax2;
      return TMath::Sqrt(r2);
    }
    double TFromWRCosThatQ20(double W, double r, double cosTh) const{
      auto mB2 = M2FromR(r);
      //std::cout<<"TFromWRCosThatQ20 "<<mB2<<" ct "<<cosTh<<" target "<<GetTarget()->M()<<" meson "<<GetMeson()->Mass()<<std::endl;
      
      return kine::tFromcosthW(cosTh,W,0.0,GetTarget()->M(),GetMeson()->Mass(),TMath::Sqrt(mB2));
    }

  
   void CalcKine() const final{
     //Requires valid Particle LorentzVectors 
     //2 body kinemtics
     TwoBodyProduction::CalcKine();

     //here we start with MB2 and must calculate r
     _mB2 = GetBaryon()->M2();
     //additional for this dsigma
      //jpac Eqn(A9)
      //meson momentum in CM
      auto qf = kine::PDK(get_W(), GetMeson()->Mass(), GetBaryon()->Mass());
      //maximimum possible meson momentum @ minimum baryon mass
      auto qmax = kine::PDK(get_W(), GetMeson()->Mass(),_minMassB);
      _r = qf/qmax;

      auto EmCM= TMath::Sqrt(qf*qf+GetMeson()->M2());

      _jacobian = 2*TMath::Pi()*_r*_r*qmax*qmax*qmax/EmCM;

//calculate mass weight here
if(_createMassEnvelope)SetMassWeight(_massEnvelope->GetCurrentWeight());
    }
    double JacobianRCosTh() const{
      auto qmax = kine::PDK(get_W(), GetMeson()->Mass(),_minMassB);
      auto qf = _r*qmax;
      auto EmCM= TMath::Sqrt(qf*qf+GetMeson()->M2());
      _mB2 = M2FromR(_r);
      
      //_jacobian = 2*TMath::Pi()*r*r*qmax*qmax*qmax/EmCM;
      return _jacobian = 2*TMath::Pi()*qf*qf*qmax/EmCM;
    }
   double PhaseSpaceFactoratQ20() const {
     //divide out phase space form jpac::inclusive
     return  sqrt(kine::Kallen(get_MB2(), get_t(), GetTarget()->M2())  / kine::Kallen(get_s(), 0.,  GetTarget()->M2()) );
    }



    jpacInc_ptr _inc={nullptr}; //I am not the owner
    double _minMassB=0;
    mutable double _mB2=0;
    mutable double _jacobian=0;
    mutable double _r=0;

    std::unique_ptr<DistTH1> _massEnvelope;
    bool _createMassEnvelope={false};
    ClassDefOverride(elSpectro::JpacSemiInc,1); //class JpacSemiInc
    
  };//class JpacSemiInc

}//namespace elSpectro
