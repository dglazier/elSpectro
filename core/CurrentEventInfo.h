//////////////////////////////////////////////////////////////
///
///Class:		CurrentEventInfo
///Description:
///             Interface to any information that decay models may depend on
///             In particular prouction kinematics and polarisations
#pragma once

namespace elSpectro{
  
  class CurrentEventInfo{


  public:

    double _weight={1};


    //variables for weighting phasespace
    /* double _emMinPhaseSpace={0}; */
    /* double _emMaxPhaseSpace={0}; */
    /* double _wgtMaxPhaseSpace={1}; */
    /* double _wgtPhaseSpace={1}; */
    
  };

  class PhotoProdInfo : public CurrentEventInfo {


  public:

    double _W={0};
    double _Q2={0};

    // double _polPhoto;
    // double _polNuc;
    
  };


}
