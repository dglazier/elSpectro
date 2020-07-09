#include "ScatteredElectron_xy.h"
#include "FunctionsForElectronScattering.h"
#include "FunctionsForGenvector.h"
#include <TRandom.h>
#include <TMath.h>

namespace elSpectro{
  
   ScatteredElectron_xy::ScatteredElectron_xy(Distribution* dist):
    _random_xy{dist}
  {
    
  }
  
  ////////////////////////////////////////////////////////////////////
  ///Caclulate electron scattering kinematics from
  ///parent electron beam in nucleon target rest frame
  ///Sample a and y which are bound [0,1] so independent of eenrgy scale
  ///Return a weight that gives phase-space distribution
  double ScatteredElectron_xy::Generate(const LorentzVector& parent, const particle_ptrs& products)  {

    //std::cout<<"ScatteredElectron_xy::Generate "<<parent.T()<<" "<<products.size()<<std::endl;
    double Ee = escat::E_el(parent.Z()); //parent in rest frame of ion, momentum=e momentum
    double Mion= parent.T()-Ee; // energy of parent = Mion + E(e-)

    double xx,yy;
    std::tie(xx,yy) = _random_xy->SamplePair();

    // std::cout<<" x = "<<xx*1E6<<" y "<<yy*1E6<<std::endl;
    double Egamma = Ee * yy;
    double Esc = Ee - Egamma;

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
    // std::cout<<"Scattered electrons "<<_scattered.M()<<" "<<Esc*Esc-Psc*Psc<<std::endl;
    // std::cout<<(Ee==Esc)<<" Ee "<<Ee*1E2<<"Esc "<<Esc*1E2<<" P "<<Psc*1E2<<" phi "<<phi <<" costh "<<costh<<std::endl;
    //std::cout<<_scattered<<" "<<_scattered.Theta()<<std::endl;
    
    //Must make sure scattered e- is in the same frame as the parent
    //still in rest system of nucl, just need rotation
    //Please note this needs checked for correct rotation
    RotateZaxisToCMDirection(parent);
    
   _gamma_ion= parent - _scattered; //residual gamma* + ion system

   //   std::cout<<_scattered<<" "<<TMath::Pi()*1E6-_scattered.Theta()*1E6<<std::endl;
  
   if(products[0]->Pdg()==11){
     products[0]->SetXYZT(_scattered.X(),_scattered.Y(),_scattered.Z(),_scattered.T());
     products[1]->SetXYZT(_gamma_ion.X(),_gamma_ion.Y(),_gamma_ion.Z(),_gamma_ion.T());
   }
   else{
     products[1]->SetXYZT(_scattered.X(),_scattered.Y(),_scattered.Z(),_scattered.T());
     products[0]->SetXYZT(_gamma_ion.X(),_gamma_ion.Y(),_gamma_ion.Z(),_gamma_ion.T());
   }
    
    return 1.;
    //    return _random_xy->CurrentValue();
  }

  
}
