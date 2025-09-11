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
      _brRatios.push_back(bratio);
      _models.push_back(std::move(mod));
      _models.back()->SetParent(parent);
      _decayers.push_back(std::move(dec));
      _idecay++; //so we can edit new decay
      std::cout<< "DecayChannel AddDecay " << _models.size()<<" "<<_models.back()<<std::endl;
     }

    void SetParent(DecayingParticle* parent){
      for(uint i = 0; i<_models.size();++i){
	SetCurrChannel(i);
	CurrModel()->SetParent(parent);
      }
      SetCurrChannel(0);
    }
    
    void PostInit(ReactionInfo* info){
      //    std::cout<<"DecayChannel::PostInit "<<dynamic_cast<ReactionElectroProd*>(info) <<std::endl; 
     for(uint i = 0; i<_models.size();++i){
	SetCurrChannel(i);
	CurrModel()->PostInit(info);
 	CurrDecayer()->PostInit(info);
      }
     MaxThreshold();
     SetCurrChannel(0);
    }
    uint ChooseDecay(double W) const;
    
    uint CurrChannel() const {return _idecay>=N() ? 0 : _idecay;}

    void SetCurrChannel(uint val) const {
      if(val>=N()){
	throw std::runtime_error(Form("DecayChannel::Threshold, asked for model %d, but only have %d",val,N()));
      }
      _idecay=val;
    };

    DecayModel* CurrModel()  const { return _models[CurrChannel()].get();}

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

    double CurrThreshold() const{
     return _models[_idecay]->MinimumMassPossible();
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
    
    double MaxThreshold() const{
      if(_maxThreshold>0.) return _maxThreshold;
      _maxThreshold=0.;
      //take the minimum possible mass of all models
      for(const auto model : _models){
	auto mm = model->MinimumMassPossible();
	if(mm>_maxThreshold) _maxThreshold = mm;
      }
      return _maxThreshold;
    }
    
  private:

    friend ProductionProcess;
    
    void SetBranchRatios(const std::vector<double>& brs) const{
      if(_brRatioSum.size()!=brs.size()+1){
	std::cerr<<"DecayChannel::SetBranchRatios must have same number of branching ratios as channels = " <<_brRatioSum.size()<< " not "<< brs.size()<<" "<<_models.size()<<std::endl;exit(0);
      }

      _brRatioSum.clear();
      _brRatios.clear();
      _brRatioSum.push_back(0.0); //sum vector must start at 0.
      for(auto br:brs){
	_brRatioSum.push_back(_brRatioSum.back()+br);
	_brRatios.push_back(br);
      }
      
    }
    std::vector<decaymodel_ptr > _models;
    std::vector<decayer_ptr  > _decayers;
    mutable std::vector<double> _brRatioSum={0.0};
    mutable std::vector<double> _brRatios={0.0};
    mutable double _maxThreshold=0.;
    mutable uint _idecay=0;
  };

 
}
