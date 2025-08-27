#include "DistFlatMass.h"

namespace elSpectro{

  DistFlatMass::DistFlatMass(DistFlatMassMaster* master):
    _master{master}
  {
    if(this!=master)_index = master->AddClient();
   }
   void DistFlatMass::Print(){
     std::cout<< "DistFlatMass master "<<(_master==this)<<" "<<_master<<" "<<_master->Parent()<<" this "<<this<<std::endl;
      // _master->Print();
     std::cout<< "index "<<_index<<" "<<_maxX<<std::endl;

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
  void DistFlatMassMaster::SetParticlePtrs(particle_ptrs ps){
    if(ps.size()!=_products.size()){
      std::cerr<<"DistFlatMassMaster::SetParticlePtrs different size to original " <<ps.size()<<" to "<<_products.size()<<std::endl;
    }
    _products=ps;
    //std::cout<<"DistFlatMassMaster::SetParticlePtrs "<<_parent->Pdg()<<std::endl;
    for(auto* p : _products ){
      //make sure I am the master!
      if( dynamic_cast<DistFlatMass*>( p->MassDistribution() )!=nullptr){
	dynamic_cast<DistFlatMass*>( p->MassDistribution() )->SetMaster(this);
     }
      //std::cout<<" "<<p->Pdg();
    }
    //std::cout<<std::endl;
  }
  double DistFlatMassMaster::SampleSingle()   noexcept{
      //might need a while loop to make sure TCM>0 when
      //sample product dymamic mass
      auto Tcm=_parent->Mass();
      uint imass=0;
      // std::cout<<"DistFlatMassMaster::SampleSingle "<<_parent->Pdg()<<" "<<Tcm<<std::endl;
      //Note TCM is calculated using all stable masses
      //If we use an unstable particle as one of the products
      //Then we cannot just take its sampled mass as this does
      //not give the correct phase space distribution
      //subtle point (took a while to debug!)
      std::vector<double> masses;
      _parent->Model()->GetStableMasses(masses);
      
      Tcm=std::accumulate(masses.begin(),masses.end(), Tcm,  std::minus<double>());
      //std::cout<<_parent->Mass()<<" DistFlatMassMaster::SampleSingle "<<Tcm<<std::endl;
      for(const auto& p : _products ){
	p->TakePdgMass();
	//throw a value of child mass if it has a distribution
	//lock it here so it cannot be changed by anyone else
	//i.e. DecayModel::PhaseSpaceWeightSq
	double sumMass=0.;
	if(p->IsDecay()!=DecayType::Stable){
	  p->TakePdgMass();
	  //add on this particle's stable masses
	  //so they can count towards its mass
	  std::vector<double> pmasses;
	  static_cast<DecayingParticle*>(p)->Model()->GetStableMasses(pmasses);

	  for(const auto &mm:pmasses){
	    sumMass+=mm;
	  }
	  Tcm+=sumMass;
	  
	  //  std::cout<<_parent->Mass()<<"DistFlatMassMaster::SampleSingle "<<sumMass<<" tcm "<<Tcm<<" size "<<pmasses.size()<<" sum "<<static_cast<DecayingParticle*>(p)->Model()->SumOfProductMasses()<<" nprods "<<static_cast<DecayingParticle*>(p)->Model()->Products().size()<<" "<<pmasses[0]<<" "<<pmasses[1]<<" "<<std::accumulate(masses.begin(),masses.end(), 0,  std::minus<double>())<<std::endl;
	}
	
	p->UnlockMass();
	//std::cout<<"DistFlatMassMaster::SampleSingle "<<_parent->Pdg()<<" "<<p->Pdg()<<" "<<Tcm<<std::endl;
	p->DetermineDynamicMass(-1,Tcm);//set a maximum to limit unphysical samples
	p->LockMass();
	_prodMass[imass++]=p->Mass();


	//if(p->IsDecay()!=DecayType::Stable){
	//sbtract product masses again
	Tcm-=sumMass;
	//	}
	//	std::cout<<"DistFlatMassMaster::SampleSingle "<<_parent->Pdg()<<" "<<p->Pdg()<<" "<<Tcm<<" "<<_prodMass[imass-1]<<" "<<_size<<std::endl;
      }
      
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
  
 //inline functions which rely on forward declaration
  double DistFlatMass::SampleSingle()   noexcept {
    if(_master->GetMass(_index)==0){
      std::cout<<"Error DistFlatMass::SampleSingle() "<<_index<<" "<< _master->GetMass(_index)<<" "<< _master->Size()<<" parent "<<_master->Parent()->Mass()<<" "<<_master->Parent()->Pdg()<<" prod0 "<<_master->Products()[0]->Pdg()<<" "<<_master->Products()[0]->Mass()<<" prod1 "<<_master->Products()[1]->Pdg()<<" "<<_master->Products()[1]->Mass()<<std::endl;
      exit(0);
    }
    return  _master->GetMass(_index);
  }
  
  double DistFlatMass::GetX() const noexcept { return _master->GetMass(_index);}

}
