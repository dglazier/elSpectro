#include "ExcitationSpectra.h"
#include "DecayingParticle.h"

namespace elSpectro{
  //////////////////////////////////////////////////////

  void ExcitationSpectra::CreateSpectra(const DecayChannel& channels,double Wmax){

    _minW = channels.Threshold();
    _maxW = Wmax;

    
    for (uint ich = 0 ; ich<channels.N() ; ++ich ){
      channels.SetCurrChannel(ich);
      _minW = channels.CurrThreshold();
      SpectraFromModel(dynamic_cast<ProductionModel*>( channels.CurrModel() ) );
    }
    // _total = SumSpectra();

    //take the absolute minimum for now
    _minW = channels.Threshold();

  }
  
  DistTH1 ExcitationSpectra::SumSpectra(){
    if(_spectra.empty()==true){
      return DistTH1(TH1D{});
    }
    TH1D sum = _spectra[0].GetTH1();
    bool missFirst=false;
    for(auto& spec:_spectra){
      if(missFirst) sum.Add( &(spec.GetTH1()) );
      missFirst=true;
    }
    return DistTH1( sum );
  }
  void ExcitationSpectra::CreateTemplateHist(){
    //W bins want focussed on threshold
    std::vector<double > WBins;
    int NW=20;//safe at 60 but much faster with 20!
    double WRange = _maxW - _minW;
    double deltaW = WRange/NW;
    for(int iW=0;iW<NW+1;++iW){
      int Nsteps = NW-iW;
      if(Nsteps==0)Nsteps==1;
      if(iW==0) Nsteps = 50; //extra at threshold
      for(int iWi=0;iWi<Nsteps;++iWi){
	WBins.push_back(_minW + iW*deltaW+static_cast<double>(iWi*deltaW)/Nsteps );
      }
    }
    WBins.push_back(_maxW);
    //increasing order
    std::sort(WBins.begin(),WBins.end());
    std::cout<<"ExcitationSpectra::CreateTemplateHist "<<WBins.size()<<" "<<WBins[0]<<" "<<WBins[1]<<std::endl;
    _template = { "XS_W","XS_W",(WBins.size()-1),WBins.data() };
   
  }
  void ExcitationSpectra::SpectraFromModel(ProductionModel* model){
    if(model==nullptr){
      std::cerr<< "ExcitationSpectra::SpectraFromModel invalid model, must inherit from ProductionModel. Exiting..." <<std::endl;
      exit(0);
      
    }
    //Note threshold changes for each channel
    //so need to redefine template
    CreateTemplateHist();
  
    // Need to generate a 1D distribution dependent on W
    // integrating overall all other variables (masses, angles)
    _spectra.push_back( DistTH1( model->CrossSectionW(_template) ) );

  }


}
