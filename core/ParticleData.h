#pragma once
#include "DecayChannel.h"
#include "DistTF1.h"
#include "TwoBodyFlat.h"
#include <TNamed.h>
#include <TDatabasePDG.h>


namespace elSpectro{
  namespace particles{
  
  
    class ParticleData : public TNamed {
    
    
    public:
      ParticleData() = default;
    
      ParticleData(const TString& name, const TString& type="");
      ParticleData(int pdg,const TString& type=""):_pdg(pdg){}

      void AddDecay(double branch_ratio,std::vector<ParticleData> hadrons);
    
      int PDGCode()const {return _pdg;}
 
      bool ShouldItDecay() const {
	//particle may decay if short lived ( <1ns)
	auto rootPdg = TDatabasePDG::Instance()->GetParticle(GetName());
	//note pi0 lifetime 7.81e-09ns
	return  rootPdg->Lifetime()<1E-9&&rootPdg->Lifetime()>0;
      }
      bool IsWide() const{
	return Width()>0.001; //minimum width 1MeV
      }
   
      void SetMass(double m){_mass = m;}
      double Mass()const {return _mass;}

      void SetWidth(double w){_width=w;};
      double Width()const {return _width;}

      bool IsDecaying(){
	return _decayChannels.empty()!=true ;
      }
      uint NDecays() const{return _decayChannels.size();}
    
      decaymodel_ptr ChannelModel(uint ichn) const;
 
      double BranchRatio(uint ichn) const{
	return _branchRatio[ichn];
      }
    
      decayer_ptr ChannelDecayer(uint ichn) const{
	//eventually should allow polarised/non-flat things here...
	auto res = elSpectro::CloneDecayer( TwoBodyFlat()  );
	std::cout<< "ParticleData::ChannelDecayer done "<<res.get()<<std::endl;
	return res;
      }
    
      double BranchRatio() const;
      
      elSpectro::DistTF1 BreitWignerDistribution(){
	return elSpectro::DistTF1{TF1("Mass",Form("TMath::BreitWigner(x,%lf,%lf)",Mass(),Width()),Mass()-5*Width(),Mass()+5*Width())};
      }

      void Print(Option_t *option="")const override;
      
    private:
    
      std::vector< std::vector<ParticleData> > _decayChannels;//
      mutable elSpectro::particle_ptrs _decayPtrs;//
    
      std::vector<double> _branchRatio;//
      double _mass=0;//
      double _width=0;//
      int _pdg=0;//
 

      //ClassDef(elSpectro::particles::ParticleData,1);
    };
  
 
 

  }
}
