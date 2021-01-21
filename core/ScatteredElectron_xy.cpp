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

    double xx,yy;
    std::tie(xx,yy) = _random_xy.SamplePair();
    CompleteGivenXandY(parent, products, xx, yy);
    double Ee = escat::E_el(parent.Z()); //parent in rest frame of ion, momentum= 
    double Mion= parent.T()-Ee; // energy of parent = Mion + E(e-)

    double Egamma = Ee * yy;
 
    double W = sqrt( Mion*(Mion + 2*Egamma ) -  escat::Q2_xy( Ee,xx,yy));
    histy.Fill(yy);
    histyQ2.Fill(escat::Q2_xy( Ee,xx,yy),yy);
    histyx.Fill(xx,yy);
    histW.Fill(W);
 
    return 1;//in this case distribution already accounts for virtual photon flux
  
    /*
    double Ee = escat::E_el(parent.Z()); //parent in rest frame of ion, momentum=e momentum
    double Mion= parent.T()-Ee; // energy of parent = Mion + E(e-)

    double xx,yy;

    std::tie(xx,yy) = _random_xy.SamplePair();

    double Egamma = Ee * yy;
 
    double W = sqrt( Mion*(Mion + 2*Egamma ) -  escat::Q2_xy( Ee,xx,yy));
   
    double Esc = Ee - Egamma;
    // if(escat::Q2_xy( Ee,xx,yy)<0.3&&escat::Q2_xy( Ee,xx,yy)>0.01) histy.Fill(yy);
    histy.Fill(yy);
    histyQ2.Fill(escat::Q2_xy( Ee,xx,yy),yy);
    histyx.Fill(xx,yy);
    histW.Fill(W);
   
    //calculate cos(theta) from e,x,y (via Q2 and Mass proton)
    double costh = escat::CosTh_xy(Ee,xx,yy);
    costh = costh>1 ? 1 : costh; //protect <=1

    histyCosTh.Fill(costh,yy);

    auto sinth=TMath::Sqrt(1-costh*costh);

    //random phi done in RotatetoParent
    double Psc=escat::P_el(Esc);
    
    auto x_sc = Psc * sinth ;
    auto y_sc = 0;
    auto z_sc = Psc * costh;

    _scattered.SetXYZT(x_sc,y_sc,z_sc,Esc);//scattered electron
    
    //Must make sure scattered e- is in the same frame as the parent
    //still in rest system of nucl, just need rotation
    RotateToParent(parent,_scattered);
    auto phot=parent -_scattered;

    if(products[0]->Pdg()==11){
      products[0]->SetP4(_scattered);
      products[1]->SetP4(parent -_scattered);
    }
    else{
      products[1]->SetP4(_scattered);
      products[0]->SetP4(parent -_scattered);
    }
    return 1;//in this case distribution already accounts for virtual photon flux
    
    */
  }

  double ScatteredElectron_xy::GenerateGivenXandY(const LorentzVector& parent, const particle_ptrs& products,double xx, double yy)  {


 
    CompleteGivenXandY(parent, products, xx, yy);
  
    return 1;//in this case distribution already accounts for virtual photon flux
    
    
  }


  double ScatteredElectron_xy::CompleteGivenXandY(const LorentzVector& parent, const particle_ptrs& products,double xx, double yy){

    double Ee = escat::E_el(parent.Z()); //parent in rest frame of ion, momentum=e momentum
    double Mion= parent.T()-Ee; // energy of parent = Mion + E(e-)

    double Egamma = Ee * yy;
 
    double W = sqrt( Mion*(Mion + 2*Egamma ) -  escat::Q2_xy( Ee,xx,yy));
   
    double Esc = Ee - Egamma;
    // if(escat::Q2_xy( Ee,xx,yy)<0.3&&escat::Q2_xy( Ee,xx,yy)>0.01) histy.Fill(yy);
   
    //calculate cos(theta) from e,x,y (via Q2 and Mass proton)
    double costh = escat::CosTh_xy(Ee,xx,yy);
    costh = costh>1 ? 1 : costh; //protect <=1

    histyCosTh.Fill(costh,yy);

    auto sinth=TMath::Sqrt(1-costh*costh);

    //random phi done in RotatetoParent
    double Psc=escat::P_el(Esc);
    
    auto x_sc = Psc * sinth ;
    auto y_sc = 0;
    auto z_sc = Psc * costh;

    _scattered.SetXYZT(x_sc,y_sc,z_sc,Esc);//scattered electron
    //Must make sure scattered e- is in the same frame as the parent
    //still in rest system of nucl, just need rotation
    RotateToParent(parent,_scattered);
    if(products[0]->Pdg()==11){
      products[0]->SetP4(_scattered);
      products[1]->SetP4(parent -_scattered);
    }
    else{
      products[1]->SetP4(_scattered);
      products[0]->SetP4(parent -_scattered);
    }
    return 1.;
  }
}
