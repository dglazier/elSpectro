//////////////////////////////////////////////////////////////
///
///Class:		DecayModel
///Description:
///             Control behaviour of Particle decay to Particle products
///             Defined by
///             1) list of Particle products
///             2) Intensity function given the produt LorentzVectors
///
///            Note drevied classes should include a constructor to initialise
///            DecayModel( particle_ptrs , const std::vector<int> pdgs );
#pragma once

#include "Particle.h"
#include "NBodyPhaseSpace.h"
#include "CurrentEventInfo.h"
#include "FunctionsForKinematics.h"
#include "ReactionInfo.h"
#include <TObject.h> //for ClassDef
#include <vector>
#include <string>
#include <algorithm>
#include <numeric> //for accumulate

namespace elSpectro{

  class DecayingParticle;
  class PhaseSpaceDecay;

  using particle_ptrs= std::vector<Particle*> ;
  using particle_const_ptrs= std::vector<const Particle*> ;
  using decaying_ptrs= std::vector<DecayingParticle*> ;
  using particle_objs= std::vector<Particle> ;
  using decaying_objs= std::vector<DecayingParticle> ;
 
  class DecayModel {

  public:
    
    DecayModel()=default;
    virtual ~DecayModel()=default;
    DecayModel(const DecayModel& other);// = default; //need the virtual destructor...so rule of 5
    DecayModel(DecayModel&& other);//=default;
    DecayModel& operator=(const DecayModel& other);// = default;
    DecayModel& operator=(DecayModel&& other);// = default;
//only declaring default constructor
    //so other 5 constructors also defaulted(rule of 5)
    void CopyOther(const DecayModel& other);
    
    //constructor to decay into particles
    //    DecayModel( particle_ptrs , const std::vector<int> pdgs );
    DecayModel( const decaying_objs& decs, const particle_objs& stables);
    
    // Each model must define its intensity
    //    virtual double Intensity() const=0;
    virtual double  Intensity() const=0;

    virtual const particle_ptrs& Products() const{ return _products;}
    const Particle*  Product(UInt_t i) const{ return (Products()[i]);}
    const LorentzVector&  ProductP4(UInt_t i) const{ return (Products()[i])->P4();}
    const decaying_objs& UnstableProducts() const { return _unstables;}
    const particle_objs& StableProducts() const { return _stables;}

    void  EventParticles(particle_ptrs& parts);

    // void AddUnstableProducts(){std::copy(std::begin(ps), std::end(ps), std::back_inserter(_products));}
    void SwapProducts(int index1, int index2){//swap pointers to particles in vector
      std::iter_swap(_products.begin() + index1, _products.begin() + index2);
    }
    
    virtual bool RegenerateOnFail() const noexcept =0;
    virtual bool HasAngularDistribution(){return true; }

    //Any preliminaries required
    virtual bool ReadyForDecay() {
      return CheckThreshold();
    }
    
    void ChooseDecay();
    
    bool CheckThreshold() const{
      SumAllProducts();
      //  std::cout<<"CheckThreshold() "<<_sumOfMasses<<" "<<_parent.M()<<std::endl;
      if(_sumOfMasses>_parent.M()) return false;
      else return true;
    }
    void SumAllProducts( ) const {
      _parent.SetXYZT(0,0,0,0);
      _sumOfMasses=0;
      
      for(auto* entry:_products){
	_parent+=entry->P4();
	_sumOfMasses+=entry->P4().M();
      }
    }

    const LorentzVector& ParentVector()const {SumAllProducts( );return _parent;}

    double SumOfProductMasses()const {return _sumOfMasses;}

    virtual double MinimumMassPossible() const;
    
    void  GetStableMasses( std::vector<double >& masses) const;
    
    virtual void Print() const;

    const std::string& GetName()const {return _name;}

    virtual double PhaseSpaceWeightSq(double W);

    void DetermineProductMasses();
    
    virtual void PostInit(ReactionInfo* info);

    virtual bool CanUseSDME()const noexcept{return false;}

    virtual void SetParent(DecayingParticle* pa){ _parentPtr=pa;}
    DecayingParticle* Parent()const noexcept{return _parentPtr;}

    virtual double dsigma() const {return 1;}

  protected:

    friend DecayingParticle;
    friend PhaseSpaceDecay;
    
    std::string _name;
    
    void ResetProducts(particle_ptrs ps);
    particle_ptrs& MutableProducts() { return _products;}
    decaying_objs& MutableUnstableProducts() { return _unstables;}
    particle_objs& MutableStableProducts()   { return _stables;}
    Particle*  MutableProduct(UInt_t i)      { return _products[i];}
 
    
  private:

    DecayingParticle* _parentPtr={nullptr};
    
    particle_ptrs _products;
    particle_objs _stables; //products which are stable
    decaying_objs _unstables; //products which decay
    
    std::vector<double> _unstableReservedMass; //mass reserved for other unstable products
    
    mutable LorentzVector _parent;
    mutable double _sumOfMasses=0;

    
    ClassDef(elSpectro::DecayModel,1); //class DecayModel
    
  };//class DecayModel

}//namespace elSpectro
