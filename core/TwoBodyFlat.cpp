#include "TwoBodyFlat.h"
#include "FunctionsForKinematics.h"
#include <TMath.h>

namespace elSpectro{


  ////////////////////////////////////////////////////////////////////
  ///Caclulate two body decay from masses and random costh and phi
  ///Return a weight that gives phase-space distribution
  double TwoBodyFlat::Generate(const LorentzVector& parent, const particle_ptrs& products)  {

    _W=parent.M(); //may also be needed by derived classes

    _weight=1;//reset weight
    
    //std::cout<<"TwoBodyFlat::Generate "<<W<<" "<<parent.E()<<" "<<products.size()<<std::endl;
    //sample mass of products in case they are from a distribution
    auto m2_a =products[0]->M2();
    auto m2_b =products[1]->M2();
    
    //std::cout<<"TwoBodyFlat "<<_W<<" "<<m2_a<<" "<<m2_b<<std::endl;
    //  auto pdk=kine::PDK(_W,products[0]->Mass(),products[1]->Mass());
    //auto edk=TMath::Sqrt(pdk*pdk+m2_a);
    
    if((_W - TMath::Sqrt(m2_a) - TMath::Sqrt(m2_b) ) < 0 ) return 0;//non physical
    //std::cout<<"TwoBodyFlat::Generate "<<_W<<" "<<m2_a<<" "<<m2_b<<" "<<" "<<((_W*_W + m2_a - m2_b))/(2.0*_W)<<" "<<((_W*_W + m2_a - m2_b))/(4.0*_W*_W) - m2_a<<" "<<pdk<<" "<<edk<<std::endl;
 
    auto e_a = (_W*_W + m2_a - m2_b)/(2.0*_W); // E decay product a
    auto p_a = TMath::Sqrt(e_a*e_a - m2_a); // p for both
    auto e_b = TMath::Sqrt(p_a*p_a + m2_b); // E for decay product b

    //std::cout<<"CM "<<e_a<<" "<<p_a<<" "<<e_b<<std::endl;
    auto costh = RandomCosTh();
    auto sinth=TMath::Sqrt(1-costh*costh);
    
    // std::cout<<"TwoBodyFlat CosTh "<<TMath::ACos(costh)<<" "<<costh<<std::endl;
    
    auto phi =  RandomPhi();   
    auto sinphi=TMath::Sin(phi);
    auto cosphi=TMath::Sqrt(1-sinphi*sinphi);

    //momentum components in CM frame
    auto x_a = p_a * sinth * cosphi;
    auto y_a = p_a * sinth * sinphi;
    auto z_a = p_a * costh;

    
    //set the local child 4-vectors
    _a.SetXYZT( x_a, y_a, z_a, e_a);
    _b.SetXYZT(-x_a,-y_a,-z_a, e_b);


    //Must make sure scattered products are in the same frame as the parent
    //still in rest system of nucl, just need rotation
    //Please note this needs checked for correct rotation
    RotateZaxisToCMDirection(parent);
    
    //boost from parent rest frame 
    auto boostFromParent=-parent.BoostToCM();
    
    _a=boost(_a,boostFromParent);//ROOT::Math::VectorUtil::boost;
    _b=boost(_b,boostFromParent);
 
     //set the lorentz vectors of the decay children
    products[0]->SetXYZT(_a.X(),_a.Y(),_a.Z(),_a.T());
    products[1]->SetXYZT(_b.X(),_b.Y(),_b.Z(),_b.T());
   
    //  std::cout<<"Masses "<< products[0]->P4().M()<<" "<< products[1]->P4().M()<<" "<<std::endl;
    //std::cout<<"Check theta "<<products[0]->P4().Theta()*TMath::RadToDeg()<<" "<<products[1]->P4().Theta()*TMath::RadToDeg()<<std::endl;
    //check
    // auto sum =_a+_b;
    //std::cout<<sum.X()<<"="<<parent.X()<<"   "<<sum.Y()<<"="<<parent.Y()<<"   "<<sum.Z()<<"="<<parent.Z()<<"   "<<sum.T()<<"="<<parent.T()<<"   "<<std::endl;
    return _weight; //this is already a pure phase space distribuition
  }
}
