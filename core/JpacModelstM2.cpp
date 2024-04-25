#include "JpacModelstM2.h"
#include "FunctionsForJpac.h"
#include "FunctionsForGenvector.h"

namespace elSpectro{
  ///////////////////////////////////////////////////////
  ///constructor includes subseqent decay of Ngamma* system
  JpacModelstM2::JpacModelstM2(jpacPhoto::inclusive_production* inc ,
			      particle_ptrs parts, const std::vector<int> pdgs) :
    _inc{inc},
    DecayModelstM2{ parts, pdgs }
  {
    _name={"JpacModelstM2"};

 
    std::cout<<"JpacModelstM2::JpacModelstM2 "<<_amp<<std::endl;
  }

 void DecayModelstM2::HistIntegratedXSection(TH1D& hist){

    auto M1 = 0;//assume real photon for calculation
    auto M2 = _target->M();
    auto M3 = _meson->Mass(); //should be pdg value here
    auto M4 = _baryon->Mass();
    auto Wmin = M3+M4;
 
      for(int ih=1;ih<=hist.GetNbinsX();ih++){
	_W=hist.GetXaxis()->GetBinCenter(ih);
	if( _W < Wmin )
	  hist.SetBinContent(ih, 0);
	else// S=W*W ds = dW.2W
	  hist.SetBinContent(ih, _inc->integrated_xsection(W*W) );
	
    }

  }
  
  
}
