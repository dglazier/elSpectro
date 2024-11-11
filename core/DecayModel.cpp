#include "Manager.h"
#include "DecayModel.h"
#include "DecayingParticle.h"
#include <algorithm>

namespace elSpectro{

  ///////////////////////////////////////////////////////////////////////
  /// construct with particles, derived classes need the same constructor
  //  DecayModel::DecayModel( particle_ptrs ps, const std::vector<int> pdgs){
  DecayModel::DecayModel( const decaying_objs& decs, const particle_objs& stables):_unstables{decs},_stables{stables}{

    //first add pre-existing particles
    //std::copy(std::begin(ps), std::end(ps), std::back_inserter(_products));
    
    //store list of unstable particles which decay
    //make our own copy
    // for(auto* prod: ps ){ 
    //   auto dp=dynamic_cast<DecayingParticle*>(prod);
    //   if(dp!=nullptr){
    // 	_unstables.push_back(*dp);//make a copy
    //  }
    //   else{
    // 	_stables.push_back(*prod);//make a copy
    //   }
    // }
    
    // //now the non decaying particles
    // auto& pman = Manager::Instance().Particles();
    // for(const auto& pdg : pdgs){
    //   _stables.push_back( Particle{pdg} );
    //   // _stables.push_back( pman.Take( new Particle{pdg} ) );
    // }

    //now vectors have all their entries
    //we can assign their ptrs to products
    _products.reserve(_stables.size()+_unstables.size());
    for(auto& prod: _stables ){
      std::cout<<"DecayModel add stable product "<<&prod<<" "<<&_stables<<std::endl;
      _products.push_back(&prod);
    }
    for(auto& prod: _unstables ){
      _products.push_back(&prod);
    }
  	    
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
    _name=other._name;
    _parentPtr=other._parentPtr;
    _stables=other._stables;
    _unstables=other._unstables;
    _unstableReservedMass=other._unstableReservedMass;
    _parent=other._parent;
    _sumOfMasses=other._sumOfMasses;

     _products.reserve(_stables.size()+_unstables.size());
    for(auto& prod: _stables ){
      std::cout<<"ASSIGNMENT DecayModel add stable product "<<&prod<<" "<<&_stables<<std::endl;
      _products.push_back(&prod);
    }
    for(auto& prod: _unstables ){
      _products.push_back(&prod);
    }
  }
  /////////////////////////////////////////////////////////////
  
  
  void  DecayModel::ResetProducts(particle_ptrs ps){
    std::cout<<"DecayModel::ResetProducts "<<ps.size()<<" "<<_products.size()<<" "<<_unstables.size()<<" "<<_stables.size()<<std::endl;
    for(auto& p: ps)std::cout<<" pdg code "<<p->Pdg()<<std::endl;
    _products.clear();
     std::cout<<"DecayModel::ResetProducts "<<ps.size()<<" "<<_products.size()<<" "<<_unstables.size()<<" "<<_stables.size()<<std::endl;
     for(auto& p: _unstables){p.Channels().SetChannel(0);std::cout<<"a pdg code "<<p.Pdg()<<" "<<p.MassDistribution()<<" model "<<p.Channels().N()<<" "<<p.Model()->Products().size()<<std::endl;}
     //_unstables.clear();
     std::cout<<"DecayModel::ResetProducts "<<ps.size()<<" "<<_products.size()<<" "<<_unstables.size()<<" "<<_stables.size()<<std::endl;
     if(_unstables.empty()==false)_unstables.erase(_unstables.begin());
      for(auto& p: _unstables){std::cout<<"a pdg code "<<p.Pdg()<<" "<<p.MassDistribution()<<std::endl;}
      _unstables.clear();
 //_stables.erase(_stables.begin());
       //for(auto& p: _stables){std::cout<<"b pdg code "<<p.Pdg()<<" "<<p.MassDistribution()<<std::endl;}
       //_stables.erase(_stables.begin());
     //     _stables.pop_back();
     //for(auto& p: _stables){std::cout<<"c pdg code "<<p.Pdg()<<" "<<p.MassDistribution()<<std::endl;}
      //_stables.erase(_stables.begin());
     //  _stables.pop_back();
     //for(auto& p: _stables){std::cout<<"d pdg code "<<p.Pdg()<<" "<<p.MassDistribution()<<std::endl;}
     //std::cout<<" pdg done "<<std::endl;
    _stables.clear();
    std::cout<<"DecayModel::ResetProducts "<<ps.size()<<" "<<_products.size()<<" "<<_unstables.size()<<" "<<_stables.size()<<std::endl;

    //store list of unstable particles which decay
    std::cout<<"DecayModel::ResetProducts "<<ps.size()<<std::endl;
    for(auto prod: ps){
      auto dp=dynamic_cast<DecayingParticle*>(prod);
      std::cout<<"DecayModel::ResetProducts "<<dp<<std::endl;
       if(dp!=nullptr)
	_unstables.push_back(*dp);//make a copy
      else
	_stables.push_back(*prod);//make a copy
    }
    //now vectors have all their entries
    //we can assign their ptrs to products
     std::cout<<"DecayModel::ResetProducts new unstables"<<_unstables.size()<<std::endl;
    for(auto& prod: _unstables ){
      _products.push_back(&prod);
    }
   std::cout<<"DecayModel::ResetProducts new stables"<<_stables.size()<<std::endl;
    for(auto& prod: _stables ){
      _products.push_back(&prod);
    }

  }

  
  void DecayModel::PostInit(ReactionInfo* info){
    if(_unstables.empty()) return;//nothing to do
    
    std::vector<double> prodmasses;
    for(auto& p:_unstables){
        p.PostInit(info);
	prodmasses.push_back(p.MinimumMassPossible());
    }

    for(uint i=0;i<_unstables.size();++i){ //for each product
      //calculate mass reserved for subsequent unstable products
      float reserve=0;
      for(uint j=i+1;j<_unstables.size();++j)
	reserve+=prodmasses[j];
      
      _unstableReservedMass.push_back(reserve);
    }

  }
 
  void  DecayModel::EventParticles(particle_ptrs& parts){
    for(auto& entry:_stables){
      parts.push_back(&entry);
    }
    //if unstable particle add its children
    for(auto& entry:_unstables){
      entry.EventParticles(parts);
    }
 
  }
  
  void  DecayModel::GetStableMasses( std::vector<double >& masses) const{
    
    //if stable particle add its mass
    for(const auto& entry:_stables)
      masses.push_back(entry.PdgMass()); 
    
    //if unstable particla add its child masses
    for(const auto& entry:_unstables)
      entry.Model()->GetStableMasses(masses);
    
  }
  void DecayModel::DetermineProductMasses(){
    for(auto& p:_unstables){
      //This should be the only call to DetermineDynamicMass in the code.
      p.DetermineProductMasses(); //first get product masses for mass Minimum
      p.DetermineDynamicMass(); //now get the particles own mass
    }
 
  }
  double DecayModel::PhaseSpaceWeightSq(double W){
    //std::cout<<GetName()<<" DecayModel::PhaseSpaceWeightSq "<<Parent()<<std::endl;
    // if(Parent()->Pdg()==-2211)std::cout<<GetName()<<" DecayModel::PhaseSpaceWeightSq start "<<MinimumMassPossible()<<" "<<W<<std::endl;
    
    //std::cout<<GetName()<<" DecayModel::PhaseSpaceWeightSq start "<<MinimumMassPossible()<<" "<<W<<" "<<_unstables.size()<<" "<<_stables.size()<<std::endl;
    //Note use weight squared to reduce sqrt calls
    
    if(_products.size()!=2){
      std::cerr<<"DecayModel::PhaseSpaceWeightSq must be 2-body decay"<<std::endl;
      exit(0);
    }
    double result=1;
    double TCM=W;
 
    for(const auto& p:_stables){
      TCM-=p.Mass();    
    }
    
    uint iu=0; //synch reserved mass vector
    for(auto& p:_unstables){
      //This should be the only call to DetermineDynamicMass in the code.
      //Unless somewhere else uses LockMass in which case
      //this call to DetermineDynamicMass will not change its value

      //Note in case there are additional unstable particle we
      //must subtract off their minimum masses
      //std::cout<<GetName()<<" DecayModel::PhaseSpaceWeightSq got Detemine dynamic mass "<<Parent()->Pdg()<<" "<<p->Pdg()<<" "<<TCM<<" "<<p->Mass()<<" "<<_unstableReservedMass[iu]<<std::endl;
      p.DetermineDynamicMass(-1,TCM-_unstableReservedMass[iu++]);
      TCM -= p.Mass();
      //std::cout<<GetName()<<"2 DecayModel::PhaseSpaceWeightSq got Detemine dynamic mass "<<Parent()->Pdg()<<" "<<p->Pdg()<<" "<<TCM<<" "<<p->Mass()<<" "<<_unstableReservedMass[iu]<<std::endl;

      if(TCM<0){
	
	return 0.;
      }//below threshold, start again
    }

    //Allow for rounding errors in check
    if(TCM<0){
      if(TCM<-1E-5){
	std::cout<<"DecayModel::PhaseSpaceWeightSq "<<Parent()->Pdg()<<" "<<GetName()<<" "<< MinimumMassPossible()<<" W "<<W<<" T "<<TCM<<" after stables "<<std::endl;}
      return 0.;
    }//below threshold (probably precission issue), start again
      
 
    result  *= kine::PDK2(W,_products[0]->Mass(),_products[1]->Mass());
    
    for(auto& p:_unstables){
      result*=p.PhaseSpaceWeightSq();
    }
    
    return result;
    
  }
  double DecayModel::MinimumMassPossible() const {
    //std::cout<<"DEcayModel::MinimumMassPossible"<<_products.size()<<" "<<_stables.size()<<" "<<&_stables<<" "<<_unstables.size()<<std::endl;
    double minmass=0;
    for(auto& p:_stables){
      //std::cout<<" stable "<<&p<<std::endl;
    }
    // for(auto& p:_unstables){
    // 	std::cout<<" unstable "<<&p<<std::endl;
    // }
    for(auto* entry:_products){
      //std::cout<<"\tDEcayModel::MinimumMassPossible "<<entry<<" "<<entry->Pdg()<<std::endl;
      minmass+=entry->MinimumMassPossible();
      //std::cout<<"\tDEcayModel::MinimumMassPossible "<<minmass<<std::endl;
	
    }
    return minmass;
  }
  void DecayModel::ChooseDecay(){
      for(auto& p:_unstables){
	p.ChooseDecay();
      }
  }

  void DecayModel::Print() const{
    std::cout<<"DecayModel::Print() "<<" "<<GetName()<<std::endl;
    for(const auto& p:_products){
      p->Print();
    }
  }

}
