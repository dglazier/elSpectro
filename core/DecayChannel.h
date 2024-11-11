//////////////////////////////////////////////////////////////
///
///Class:		DecayChannel
///Description:
///            

#pragma once

#include "DecayModel.h"
#include "DecayVectors.h"

namespace elSpectro{

  using decaymodel_ptr = std::shared_ptr<DecayModel>;
  using decayer_ptr = std::shared_ptr<DecayVectors>;
  
  template <class T>
     decaymodel_ptr CloneModel(T t){
     decaymodel_ptr m{new T(t)};
    return m;
  }
  template <class T>
    decayer_ptr CloneDecayer(T t){
    decayer_ptr v{new T(t)};
    return v;
  }
 
  /* template <class T> */
  /*   decayer_ptr CloneDecayer(T* t){ */
  /*   decayer_ptr v{new T(*t)}; */
  /*   return v; */
  /* } */
  
  class DecayChannel {
  public:
    
    void AddDecay(DecayingParticle* parent,double bratio,decaymodel_ptr mod,decayer_ptr dec){
      _brRatioSum.push_back(_brRatioSum.back()+bratio);
      _models.push_back(std::move(mod));
      _models.back()->SetParent(parent);
      _decayers.push_back(std::move(dec));
       _idecay++; //so we can edit new decay
    }

    void SetParent(DecayingParticle* parent){
      for(uint i = 0; i<_models.size();++i){
	SetChannel(i);
	CurrModel()->SetParent(parent);
      }
      SetChannel(0);
    }
    
    void PostInit(ReactionInfo* info){
     std::cout<<"DecayChannel::PostInit "<<dynamic_cast<ReactionElectroProd*>(info) <<std::endl; 
     for(uint i = 0; i<_models.size();++i){
	SetChannel(i);
	CurrModel()->PostInit(info);
 	CurrDecayer()->PostInit(info);
      }
      SetChannel(0);
    }
    uint ChooseDecay() const{
      //choose random number between 0 and sum of all branch ratios
      //find index corresponding to that value
      return _brRatioSum.size() == 1 ? _idecay :
	(_idecay = (std::lower_bound(_brRatioSum.begin(),_brRatioSum.end(),gRandom->Uniform(0,_brRatioSum.back()))) - _brRatioSum.begin() -1) ; //return 0 if 1 decay, if not choose.
    }
    uint CurrChannel() const {return _idecay>=N() ? 0 : _idecay;}

    void SetChannel(uint val) const {_idecay=val;};

    DecayModel* CurrModel()  const {return _models[CurrChannel()].get();}

    DecayVectors* CurrDecayer() const {return _decayers[CurrChannel()].get();}

    uint N() const {return _models.size();}
    
    void SetDecayer(uint index,decayer_ptr   dec){
      _idecay=index;
      _decayers[index]=std::move(dec);
    }

    void SetModel(uint index,decaymodel_ptr  mod){
      _idecay=index;
      _models[index]=std::move(mod);
    }

    double Threshold() const{
      double threshold = 1E30;
      //take the minimum possible mass of all models
      for(const auto model : _models){
	auto mm = model->MinimumMassPossible();
	if(mm<threshold) threshold = mm;
      }
      return threshold;
    }
    
  private:
    
    std::vector<decaymodel_ptr > _models;
    std::vector<decayer_ptr  > _decayers;
    std::vector<double> _brRatioSum={0.0};
    mutable uint _idecay=0;
  };
}
