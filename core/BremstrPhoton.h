//////////////////////////////////////////////////////////////
///
///Class:		BremstrPhoton
///Description:
///            Derive methods for 
///            Generating LorentzVectors for Nucleus->nucleon + spectator 
#pragma once

#include "DecayVectors.h"
#include "Distribution.h"
#include "DistTF1.h"
#include "DistTH1.h"
#include "LorentzVector.h"
#include <TRandom.h> //for gRandom
#include <Math/VectorUtil.h> //for boosts etc.

namespace elSpectro{

  using ROOT::Math::VectorUtil::boost;
  
  class BremstrPhoton : public DecayVectors {
    
  public:

    //default deuteron Fermi distribution, can overwrite when create object 
    BremstrPhoton(Double_t ebeam,Distribution* dist=new DistTF1{TF1("bremsstrahlungDist","1.0/x[0]*(4./3 - 4./3*x[0] + x[0]*x[0])",0.,1)}):_bremDist{dist}, _ebeam{ebeam}{
      auto tf1=dynamic_cast< DistTF1*>(_bremDist.get());
      if(tf1!=nullptr){
	//tf1->GetTF1().SetParameter(0,ebeam);
      }
     };
    BremstrPhoton(Double_t ebeam,Double_t emin,Double_t emax,Distribution* dist=new DistTF1{TF1("bremsstrahlungDist","1.0/x[0]*(4./3 - 4./3*x[0] + x[0]*x[0])",0,1)}):_bremDist{dist}, _ebeam{ebeam}{
      auto tf1=dynamic_cast< DistTF1*>(_bremDist.get());
      if(tf1!=nullptr){
	//	tf1->GetTF1().SetParameter(0,ebeam);
	tf1->GetTF1().SetRange(emin/ebeam,emax/ebeam);
      }
   
    };
    BremstrPhoton(Double_t ebeam,Double_t emin,Double_t emax,Distribution* poldist,double polplane, Distribution* dist=new DistTF1{TF1("bremsstrahlungDist","1.0/x[0]*(4./3 - 4./3*x[0] + x[0]*x[0])",0,1)}):_bremDist{dist}, _linPolDist{poldist},_polPlane{polplane},_ebeam{ebeam}{
      auto tf1=dynamic_cast< DistTF1*>(_bremDist.get());
      if(tf1!=nullptr){
	//	tf1->GetTF1().SetParameter(0,ebeam);
	tf1->GetTF1().SetRange(emin/ebeam,emax/ebeam);
      }
   
    };
   
    virtual ~BremstrPhoton()=default;
    BremstrPhoton(const BremstrPhoton& other); //need the virtual destructor...so rule of 5
    BremstrPhoton(BremstrPhoton&&)=default;
    BremstrPhoton& operator=(const BremstrPhoton& other);
    BremstrPhoton& operator=(BremstrPhoton&& other) = default;
 

    double Generate(const LorentzVector& parent,
		    const particle_ptrs& products)  final;

  
    double GetBeamEnergy() const {return _ebeam;}
    double GetMinEnergy() const {return _ebeam*_bremDist->GetMinX();}
    double GetMaxEnergy() const {return _ebeam*_bremDist->GetMaxX();}
    DistTF1* GetPhotonEdist() const { auto tf1=dynamic_cast< DistTF1*>(_bremDist.get()); return tf1;}
    
    void PostInit(ReactionInfo* info);

  protected :

 
  private:

    Double_t _ebeam=0;
    
    LorentzVector _photon;

    std::unique_ptr<Distribution> _bremDist;
    std::unique_ptr<Distribution> _linPolDist;

    PhotonPolarisationVector* _photonPol={nullptr};
    double _polPlane={0};
    
    ClassDef(elSpectro::BremstrPhoton,1); //class DecayVectors
 

  };

   
}
  
