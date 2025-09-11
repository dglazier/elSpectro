#include "Manager.h"
#include "DecayModel.h"
#include "DecayingParticle.h"
#include <algorithm>

namespace elSpectro{

  ///////////////////////////////////////////////////////////////////////
  /// construct with particles, derived classes need the same constructor
  //  DecayModel::DecayModel( particle_ptrs ps, const std::vector<int> pdgs){
  DecayModel::DecayModel( const decaying_objs& decs, const particle_objs& stables):_unstables{decs},_stables{stables}{
  
  
    //we can copy their objects to list of all products
    _products.reserve(_stables.size()+_unstables.size());
    for(auto& prod: _stables ){
       _products.push_back(&prod);
    }
    for(auto& prod: _unstables ){
      //std::cout<<"DecayModel add unstable product "<<prod.Pdg()<<" "<<&prod<<" "<<&_unstables<<" "<<prod.Model()->Product(0)->Pdg()<<" "<<prod.Model()->Product(0)<<std::endl;
     _products.push_back(&prod);
    }
    //   if(decs.size())std::cout<<"DecayModel dec "<<" "<<&decs[0]<<" "<<decs.size()<<" "<<decs[0].Pdg()<<" "<<decs[0].Model()->Product(0)->Pdg()<<std::endl;//<<decs[0].Model()->Product(0)->Pdg()
  }
  /////////////////////////////////////////////////////
  DecayModel::DecayModel(const DecayModel& other){
    CopyOther(other);
  }
  /////////////////////////////////////////////////////
  DecayModel::DecayModel(DecayModel&& other){
    CopyOther(other);
    _products.clear();
    _stables.clear();
    _unstables.clear();
    _parentPtr=nullptr;
  }
  //////////////////////////////////////////////////////
  DecayModel& DecayModel::operator=(const DecayModel& other) {
    CopyOther(other);
    return *this;
  }
  //////////////////////////////////////////////////////
  DecayModel& DecayModel::operator=(DecayModel&& other){
    CopyOther(other);
    _products.clear();
    _stables.clear();
    _unstables.clear();
    _parentPtr=nullptr;
    return *this;
  }
  ////////////////////////////////////////////////////////////
  void DecayModel::CopyOther(const DecayModel& other){

    //std::cout<<"DecayModel::CopyOther " <<this<<" "<<&other<<std::endl;
    _name=other._name;
    _parentPtr=other._parentPtr;
    _stables=other._stables;
    _unstables=other._unstables;
    _unstableReservedMass=other._unstableReservedMass;
    _parent=other._parent;
    _sumOfMasses=other._sumOfMasses;

     _products.reserve(_stables.size()+_unstables.size());
    for(auto& prod: _stables ){
      // std::cout<<"ASSIGNMENT DecayModel add stable product "<<prod.Pdg()<<" "<<&prod<<" "<<&_stables<<std::endl;
      _products.push_back(&prod);
    }
    for(auto& prod: _unstables ){
      //std::cout<<"ASSIGNMENT DecayModel add unstable product "<<prod.Pdg()<<" "<<&prod<<" "<<std::endl;
     _products.push_back(&prod);
    }
  }
  /////////////////////////////////////////////////////////////
  
  
  void  DecayModel::ResetProducts(particle_ptrs ps){
    _products.clear();
    
    //First create temp lists to copy ps particles
    //into stable and unstable vectors
    //this must be done before clearing current stable and unstable
    //as this may delete the ptr in ps
    particle_objs temp_stables;
    decaying_objs temp_unstables;
    for(auto prod: ps){
      auto dp=dynamic_cast<DecayingParticle*>(prod);
      if(dp!=nullptr){
	temp_unstables.push_back(*dp);//make a copy
      }
       else{
	 temp_stables.push_back(*prod);//make a copy
       }
    }
     if(_stables.empty()==false)_stables.clear();
     if(_unstables.empty()==false)_unstables.clear();
     _stables = std::move(temp_stables);  
     _unstables = std::move(temp_unstables);  
    //now vectors have all their entries
    //we can assign their ptrs to products
     _products.reserve(_stables.size()+_unstables.size());
     for(auto& prod: _unstables ){
      _products.push_back(&prod);
    }
     for(auto& prod: _stables ){
      _products.push_back(&prod);
    }
 
  }

  
  void DecayModel::PostInit(ReactionInfo* info){
    //std::cout<<"DecayModel::PostInit "<<_unstables.empty()<<std::endl;
    if(_unstables.empty()) return;//nothing to do
    
    std::vector<double> prodmasses;
    for(auto& p:_unstables){
        p.PostInit(info);
	prodmasses.push_back(p.MinimumMassPossible());
    }

    for(uint i=0;i<_unstables.size();++i){ //for each product
      //calculate mass reserved for subsequent unstable products
      float reserve=0;
      for(uint j=i+1;j<_unstables.size();++j){
	reserve+=prodmasses[j];
      }
      _unstableReservedMass.push_back(reserve);
    }
    // std::cout<<"DecayModel::PostInit done "<<std::endl; 
 
  }
 
  void  DecayModel::EventParticles(particle_ptrs& parts){
    auto nstables = _stables.size();
    for(size_t istab=0;istab<nstables;++istab){
      // std::cout<<"DecayModel add stable product "<<&prod<<" "<<&_stables<<std::endl;
      if(std::count(_notEventParticles.begin(), _notEventParticles.end(), istab)!=0) continue;
      parts.push_back(&_stables[istab]);
    }
      //for(auto& entry:_stables){
      //parts.push_back(&entry);
      //}
    //if unstable particle add its children
    for(auto& entry:_unstables){
      entry.EventParticles(parts);
    }
 
  }
  
  void  DecayModel::GetStableMasses( std::vector<double >& masses) const{

      //if stable particle add its mass
    for(const auto& entry:_stables){
      // std::cout<<"DecayModel::GetStableMasses stable "<<entry.Pdg()<<" "<<entry.PdgMass()<<std::endl;
      masses.push_back(entry.PdgMass()); 
    }
    //if unstable particla add its child masses
    for(const auto& entry:_unstables){
      //if mass is a delta function, count as a stable particle
      //or else phase space calcualtion very inefficent
      // std::cout<<"DecayModel::GetStableMasses unstable "<<Parent()->Pdg()<<" "<<entry.Pdg()<<" "<<entry.PdgMass()<<" no distribution "<<(entry.MassDistribution()==nullptr)<<std::endl;
      if(entry.MassDistribution()==nullptr){
      //if(entry.IsDecaying()==false){
	masses.push_back(entry.PdgMass()); 
   	continue;
      }
       //entry.MassDistribution()->Print();
      entry.Model()->GetStableMasses(masses);
    }
  }
  void DecayModel::DetermineProductMasses(){
    for(auto& p:_unstables){
      //This should be the only call to DetermineDynamicMass in the code.
      p.DetermineProductMasses(); //first get product masses for mass Minimum
      p.DetermineDynamicMass(); //now get the particles own mass
    }
 
  }
  bool DecayModel::SampleMasses(double W){
    //while(SampleMassesIteration(W)==false){}
    return SampleMassesIteration(W);
  }
  bool DecayModel::SampleMassesIteration(double W){
  ///std::cout<<GetName()<<" DecayModel::PhaseSpaceWeightSq "<<Parent()->Pdg()<<std::endl;
    // if(Parent()->Pdg()==-2211)std::cout<<GetName()<<" DecayModel::PhaseSpaceWeightSq start "<<MinimumMassPossible()<<" "<<W<<std::endl;
    
    // std::cout<<GetName()<<" DecayModel::PhaseSpaceWeightSq start "<<MinimumMassPossible()<<" "<<W<<" "<<MinimumMassPossible()-W<<" "<<_unstables.size()<<" "<<_stables.size()<<std::endl;
    //Note use weight squared to reduce sqrt calls
    
    if(_products.size()!=2){
      std::cerr<<"DecayModel::PhaseSpaceWeightSq must be 2-body decay"<<std::endl;
      exit(0);
    }

    double TCM=W;
 
    for(const auto& p:_stables){
      // p.TakePdgMass();
      // if(p.Pdg()==22&&p.Mass()>0)
     if(p.Pdg()==22) continue;
     //if(Parent()->Pdg()==-323)  std::cout<< "DecayModel::SampleMass stables "<<p.Mass()<<" "<<p.Pdg()<<" "<<TCM<<std::endl;
     TCM-=p.Mass();    
    }
    
    uint iu=0; //synch reserved mass vector
    for(auto& p:_unstables){
      //This should be the only call to DetermineDynamicMass in the code.
      //Unless somewhere else uses LockMass in which case
      //this call to DetermineDynamicMass will not change its value

      //Note in case there are additional unstable particle we
      //must subtract off their minimum masses
      // std::cout<<GetName()<<" DecayModel::PhaseSpaceWeightSq get Detemine dynamic mass "<<W<<" parent "<<Parent()->Pdg()<<" "<<Parent()->MinimumMassPossible()<<" particle "<<p.Pdg()<<" tcm "<<TCM<<" old mass "<<p.Mass()<<" diff "<<TCM-p.MinimumMassForChannel()<<" min mass "<<p.MinimumMassPossible()<<" min channel "<<p.MinimumMassForChannel()<<" "<<_unstableReservedMass[iu]<<" max mass "<<TCM-_unstableReservedMass[iu]<<std::endl;
      p.DetermineDynamicMass(-1,TCM-_unstableReservedMass[iu++]);
      ///p.DetermineDynamicMass(-1,TCM-p.MinimumMassForChannel());
      TCM -= p.Mass();
      //std::cout<<GetName()<<"DecayModel::SampleMass got Detemine dynamic mass "<<W <<" "<<Parent()->Pdg()<<" "<<p.Pdg()<<" tcm "<<TCM<<" pmass "<<p.Mass()<<" "<<p.MinimumMassForChannel()<<" channels "<<dynamic_cast<DecayingParticle*>(Parent())->Channels().CurrChannel()<<" "<<dynamic_cast<DecayingParticle*>(&p)->Channels().CurrChannel()<<" "<<_unstableReservedMass[iu-1]<<"\n";

      if(TCM<0){
	if(TCM>-1E-6){
	  //exit(0);
	  std::cout<<"DecayModel::PhaseSpaceWeightSq not good "<<W<<" "<<TCM<<" "<<Parent()->Pdg()<<" "<<Parent()->Mass()<<" "<<p.Pdg()<<" "<<p.Mass()<<"channel"<<std::endl;
	  TCM=0.;
	}
	// for(const auto& sp:_stables){
	//   std::cout<<" mass stable "<<sp.Mass()<<std::endl;;    
	// }
 	//return false;
      }//below threshold, start again
    }

    //Allow for rounding errors in check
    if(TCM<0){
      if(TCM>-1E-6){
       	std::cout<<"DecayModel::PhaseSpaceWeightSq not good "<<Parent()->Pdg()<<" "<<GetName()<<" "<< MinimumMassPossible()<<" W "<<W<<" T "<<TCM<<" after stables "<<std::endl;
      }
      else{
 	std::cout<<"DecayModel::PhaseSpaceWeightSq not good "<<Parent()->Pdg()<<" "<<GetName()<<" "<< MinimumMassPossible()<<" W "<<W<<" T "<<TCM<<" after stables "<<Parent()->MinimumMassForChannel()<<" "<<dynamic_cast<DecayingParticle*>(Parent())->Channels().CurrChannel()<<std::endl;
	return false;
      }
    }//below threshold (probably precission issue), start again

    //masses of decaying products
    // for(auto& p:_unstables){
    //   //if product is 
    //   //if(p.MassDistribution()!=nullptr){
    // 	p.SampleMasses(p.Mass());
    // 	//std::cout<<"DecayModel::PhaseSpaceWeightSq unstable "<<p.Pdg()<<" "<<result<<" "<<" "<<_unstables.size()<<" "<<Parent()->Pdg()<<" "<<W<<" "<<sqrt(temp)<<" "<<std::endl;
    //   }
    return true;
  }
  
  double DecayModel::PhaseSpaceWeightSq(double W,bool resample){
    if(resample)SampleMasses(W);//failed phase space, need to resample

    double result=1.;
    // if(Parent()->MassDistribution()!=nullptr){
    result  *= kine::PDK2(W,_products[0]->Mass(),_products[1]->Mass());
    // std::cout<<"DecayModel::PhaseSpaceWeightSq "<<result<<" "<<_products[0]->Mass()<<" "<<_products[1]->Mass()<<" "<<_unstables.size()<<" "<<Parent()->Pdg()<<" "<<W<<" "<<sqrt(kine::PDK2(W,_products[0]->Mass(),_products[1]->Mass()))<<std::endl;
      // }
  
    for(auto& p:_unstables){
      //if product is 
      if(p.MassDistribution()!=nullptr){
	auto temp = p.PhaseSpaceWeightSq(resample);
	result*=temp;
	//std::cout<<"DecayModel::PhaseSpaceWeightSq unstable "<<p.Pdg()<<" "<<result<<" "<<" "<<_unstables.size()<<" "<<Parent()->Pdg()<<" "<<W<<" "<<sqrt(temp)<<" "<<std::endl;
      }
    
      else{//call function to assign product masses, but we don't need the weight
	// std::cout<<"DecayModel::PhaseSpaceWeightSq just get masses"<<std::endl;
        p.PhaseSpaceWeightSq(resample);
      }
      
    }
    //  std::cout<<"DecayModel::PhaseSpaceWeightSq "<<result<<std::endl;
    return result;
    
  }
  double DecayModel::MinimumMassPossible() const {
    // std::cout<<"DEcayModel::MinimumMassPossible "<<_products.size()<<" "<<_stables.size()<<" "<<&_stables<<" "<<_unstables.size()<<std::endl;
    //  double minmass=0;
    // for(auto& p:_stables){
    //   //std::cout<<" stable "<<&p<<std::endl;
    // }
    // for(auto& p:_unstables){
    // 	std::cout<<" unstable "<<&p<<std::endl;
    // }
    if(_threshold>0.) return _threshold;
    _threshold=0.;
    for(auto* entry:_products){
      //std::cout<<"\tDEcayModel::MinimumMassPossible "<<entry<<" "<<entry->Pdg()<<std::endl;
      //     if(entry->Pdg()!=2212&&entry->Pdg()!=223) exit(0);
      _threshold+=entry->MinimumMassPossible();
       //std::cout<<"\tDEcayModel::MinimumMassPossible "<<minmass<<std::endl;
	
    }
    return _threshold;
  }
  void DecayModel::ChooseDecay() const{
      for(const auto& p:_unstables){
	p.ChooseDecay();
	//	std::cout<< "DecayModel::ChooseDecay()  " <<p.Pdg()<<" "<< p.Mass()<<std::endl;
	//p.Model()->SampleMasses(p.Mass());
      }
  }
  
  void DecayModel::SetParent(DecayingParticle* pa){
      _parentPtr=pa;
      //pass on to products in case then need
      //for example for DistFlatMass
      for(auto prod : _products){
	prod->SetParent(static_cast<Particle*>(pa));
      }
    }
  
  void DecayModel::Print() const{
    std::cout<<"DecayModel::Print() "<<" "<<GetName()<<std::endl;
    for(const auto& p:_products){
      p->Print();
    }
  }

}
