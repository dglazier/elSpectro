#include "DistFlatMass.h"

namespace elSpectro{

  DistFlatMass::DistFlatMass(DistFlatMassMaster* master):
    _master{master}
  {
    if(this!=master)_index = master->AddClient();
   }
  //////////////////////////////////////////////////////////
  DistFlatMassMaster::DistFlatMassMaster(DecayingParticle* original, particle_ptrs ps):
    DistFlatMass(this),
    _parent{original},
    _products{ps}
  {
    _invMass.resize(_products.size()-2);
    _prodMass.resize(_products.size());
    SetIndex(_invMass.size()-1);//last entry in invMass
    _size++;
  }
  
    double DistFlatMassMaster::SampleSingle()   noexcept{
      //might need a while loop to make sure TCM>0 when
      //sample product dymamic mass
      auto Tcm=_parent->Mass();
      uint imass=0;
      for(const auto& p : _products ){
	//throw a value of child mass if it has a distribution
	//lock it here so it cannot be changed by anyone else
	//i.e. DecayModel::PhaseSpaceWeightSq
	p->UnlockMass();
	p->DetermineDynamicMass(-1,Tcm);//set a maximum to limit unphysical samples
	p->LockMass();
	_prodMass[imass++]=p->Mass();
     }
      
      //Note TCM is calculated using all stable masses
      //If we use an unstable particle as one of the products
      //Then we cannot just take its sampled mass as this does
      //not give the correct phase space distribution
      //subtle point (took a while to debug!)
      std::vector<double> masses;
      _parent->Model()->GetStableMasses(masses);
      Tcm=std::accumulate(masses.begin(),masses.end(), Tcm,  std::minus<double>());

      double sum = _prodMass[0];
      
      int nrand=_size;
      double randArray[nrand];
      gRandom->RndmArray(nrand,randArray);
      //Sorting gives factor 2 speed up (probably due to unphysical values being found earlier)
      if(nrand>1)std::sort(randArray,randArray + nrand);

      for (uint n=0; n< _size; ++n) {
	sum      += _prodMass[n+1];
	_invMass[n] = randArray[n]*Tcm + sum;
      }      
      
      return _invMass[Index()];
    }

}
