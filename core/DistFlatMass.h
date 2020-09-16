//////////////////////////////////////////////////////////////
///
///Class:		DistFlatMass
///Description:
///             wrapper for TF1 distributions
/// 

#pragma once

#include "Distribution.h"
//#include "Particle.h"
#include "DecayModel.h"
#include "DecayingParticle.h"
//
#include <TF1.h>
#include <string>

namespace elSpectro{

  class DistFlatMassMaster;
  
  class DistFlatMass : public Distribution {

    
  public :

    DistFlatMass(DistFlatMassMaster* master);
 
    double SampleSingle()   noexcept override;
    
    dist_pair SamplePair()   noexcept final {return dist_pair{0,0};} ;

    double CurrentValue() const noexcept final {return 1.;}
 
    double GetX() const noexcept ;
    
    double GetWeightFor(double valX) const {return 1;} //phase space ==1, handled by MassPhaseSpace
    
    double MaxValue() const noexcept final{return 1;}
    double MinValue()  const noexcept final{return 1;}

    double GetMinX() const noexcept final{return 0;}
    double GetMaxX() const noexcept final{return 1;}

  protected :
    void SetIndex(uint index){_index=index;}
    uint Index() const noexcept {return _index;}
  private:
    //no one should use default constructor
    DistFlatMass()=default;
    
    DistFlatMassMaster* _master={nullptr}; 
    uint _index = {0};

      
    ClassDef(elSpectro::DistFlatMass,1); //class Distribution
 

  };

  class DistFlatMassMaster : public DistFlatMass {

    
  public :

    DistFlatMassMaster(DecayingParticle* original, particle_ptrs ps);

    double GetMass(uint i)const noexcept{
      return _invMass[i];
    }

    double SampleSingle()   noexcept final ;

    uint AddClient(){
  

      //SetIndex(_size);//increment my own index
      _size++;
      auto newIndex=_size-2;
      if(newIndex==Index()){ //my index is _invMass.size()  -1
	std::cerr<<"AddClient::DistFlatMassMaster, can't have index greater than mine "<<Index()<<std::endl;
	exit(0);
      }
      
      return  newIndex; //index of client
      //by time all clients added _size = _products.size()-1
    }
    
  private:
    
    DecayingParticle* _parent={nullptr};
    std::vector<double> _invMass;
    std::vector<double> _prodMass;
    std::vector<double> _stableMass;

    particle_ptrs _products={nullptr};

    uint _size={0};
  };


  //inline functions which rely on forward declaration
  inline  double DistFlatMass::SampleSingle()   noexcept {
    return  _master->GetMass(_index);
  }
  
  inline double DistFlatMass::GetX() const noexcept { return _master->GetMass(_index);}

}
