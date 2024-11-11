//////////////////////////////////////////////////////////////
///
///Class:	ExcitationSpectra
///Description:
///             For a given Formation, store all the W dependent
///             excitation spectrums given by TwoBodyProduction
///             Given W select 2-body channel for this event
#pragma once

#include "DecayChannel.h"
#include "ProductionModel.h"
#include "DistTH1.h"


namespace elSpectro{

  class ExcitationSpectra {

  public :
    
    ExcitationSpectra()=default;

    uint ChooseChannel(double W){return 0;}

    void CreateSpectra(const DecayChannel& channels,double Wmax);
    void CreateTemplateHist();
    void SpectraFromModel(ProductionModel* model);
    DistTH1 SumSpectra();
    const DistTH1& TotalCrossSection()const {return _total;}

  private:
    
    std::vector<DistTH1> _spectra;
    TH1D _template;
    DistTH1 _total;
    double _minW=0;
    double _maxW=0;
    
  };

}
