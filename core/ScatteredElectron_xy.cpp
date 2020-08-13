#include "ScatteredElectron_xy.h"
#include "FunctionsForElectronScattering.h"
#include "FunctionsForGenvector.h"
#include "Manager.h"
#include <TRandom.h>
#include <TMath.h>

namespace elSpectro{
  
  ScatteredElectron_xy::ScatteredElectron_xy(double eb, double mion, double Wmin):
    _random_xy{eb,mion,Wmin}
  {
   }
  
  ////////////////////////////////////////////////////////////////////
  ///Caclulate electron scattering kinematics from
  ///parent electron beam in nucleon target rest frame
  ///Sample a and y which are bound [0,1] so independent of eenrgy scale
  ///Return a weight that gives phase-space distribution
  double ScatteredElectron_xy::Generate(const LorentzVector& parent, const particle_ptrs& products)  {

     double Ee = escat::E_el(parent.Z()); //parent in rest frame of ion, momentum=e momentum
    double Mion= parent.T()-Ee; // energy of parent = Mion + E(e-)

    double xx,yy;

    std::tie(xx,yy) = _random_xy.SamplePair();

    double Egamma = Ee * yy;
 
    double W = sqrt( Mion*(Mion + 2*Egamma ) -  escat::Q2_xy( Ee,xx,yy));
    /*
    Manager::Instance().AcceptPhaseSpace(W);
    while(Manager::Instance().AcceptPhaseSpace(W) == false ){
      //try again
      _gStarNmodel->DetermineProductMasses();
      _random_xy.SetWThreshold(_gStarNmodel->MinimumMassPossible());
       std::tie(xx,yy) = _random_xy.SamplePair();
      Egamma = Ee * yy;
      W = sqrt( Mion*(Mion + 2*Egamma ) -  escat::Q2_xy( Ee,xx,yy));
      
      }
    */
    
    double Esc = Ee - Egamma;
    histy.Fill(yy);
    histW.Fill(W);

    //calculate cos(theta) from e,x,y (via Q2 and Mass proton)
    double costh = escat::CosTh_xy(Ee,xx,yy);
    costh = costh>1 ? 1 : costh; //protect <=1
     
    auto sinth=TMath::Sqrt(1-costh*costh);

    //random phi
    double phi=gRandom->Uniform()*2*TMath::Pi();
    auto cosphi=TMath::Cos(phi);
    auto sinphi=TMath::Sin(phi);

    double Psc=escat::P_el(Esc);

 									 
    auto x_sc = Psc * sinth * cosphi;
    auto y_sc = Psc * sinth * sinphi;
    auto z_sc = Psc * costh;

    _scattered.SetXYZT(x_sc,y_sc,z_sc,Esc);//scattered electron
    
    //Must make sure scattered e- is in the same frame as the parent
    //still in rest system of nucl, just need rotation
    //Please note this needs checked for correct rotation
    RotateZaxisToCMDirection(parent);
    
   _gamma_ion= parent - _scattered; //residual gamma* + ion system
   
   if(products[0]->Pdg()==11){
     products[0]->SetXYZT(_scattered.X(),_scattered.Y(),_scattered.Z(),_scattered.T());
     products[1]->SetXYZT(_gamma_ion.X(),_gamma_ion.Y(),_gamma_ion.Z(),_gamma_ion.T());
   }
   else{
     products[1]->SetXYZT(_scattered.X(),_scattered.Y(),_scattered.Z(),_scattered.T());
     products[0]->SetXYZT(_gamma_ion.X(),_gamma_ion.Y(),_gamma_ion.Z(),_gamma_ion.T());
   }
 
   return 1.;//in this case distribition already accounts for virtual photon flux
  }

  
}
