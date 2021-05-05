#include "QuasiFreeNucleon.h"
#include <TMath.h>

namespace elSpectro{


  ////////////////////////////////////////////////////////////////////
  ///Caclulate two body decay from masses and random costh and phi
  ///Return a weight that gives phase-space distribution
  double QuasiFreeNucleon::Generate(const LorentzVector& parent, const particle_ptrs& products)  {
    _weight=1;
    //Get Fermi momentum sampled from distribution given in input file
    auto ptar = _fermiDist->SampleSingle();
    //Get random flat cos theta
    auto costhtar   = RandomCosTh();
    auto sinthtar      = TMath::Sqrt( 1-costhtar*costhtar );
    //Get random flat phi
    //  done in boost
    //  phtar      = RandomPhi();
    //calculate momentum components
    auto pxtar      = ptar * sinthtar;
    auto pytar      = 0;
    auto pztar      = ptar * costhtar;
  
    // Force spectator on mass shell
    //get its mass from input data and PDG database
    auto smass2=products[1]->M2(); //Note interacting particle is [0]
    auto Espec = TMath::Sqrt(ptar*ptar + smass2);
    //set spectator 4 momentum (opposite p to quasi target)
    _spectator.SetXYZT(-pxtar,-pytar,-pztar, Espec);
  
    //calculate quasi target energy from energy conservation
    //This will give offshell mass
    auto Etar  = parent.M() - Espec;
    _nucleon.SetXYZT(pxtar,pytar,pztar, Etar);


    //Must make sure scattered products are in the same frame as the parent
    //theta=0 => moving along parent direction
    //i.e. boost vector should only have z component
    //Please note this needs checked for correct rotation
    //Also add the random phi angle...
    //std::cout<<"spectator "<<_spectator<<" "<<_spectator.M()<<std::endl;
    //std::cout<<"nucleon "<<_nucleon<<" "<<_nucleon.M()<<std::endl;
    //_a.SetXYZT( x_a, y_a, z_a, e_a);
    BoostToParent(parent,_spectator);
    // std::cout<<"spectator boosted "<<_spectator<<" "<<_spectator.M()<<std::endl;
    products[1]->SetP4(_spectator);
    products[0]->SetP4( _nucleon= parent - _spectator );
    // std::cout<<"nucleon boosted"<<_nucleon<<" "<<_nucleon.M()<<std::endl;
    if(_nucleon.M()<0.01) _weight=0; //fix minimum offshell mass to 10MeV...
    return _weight; 
  }


  namespace QuasiFree{
    ///////////////////////////////////////////////////
    //Calculation of nucleon momentum distribution
    //based on CDBonn potential
    //code supplied by  Mikhail Bashkanov
    DistTH1* CDBonnMomentum(const Int_t Nbins,const Double_t pmin,const Double_t pmax){

      Double_t CC;
      Double_t c[11],m[11];

      c[0] = 0.88472985;
      c[1] = -0.26408759;
      c[2] = -0.044114404;
      c[3] = -14.397512;
      c[4] = 85.591256;
      c[5] = -318.76761;
      c[6] = 703.36701;
      c[7] = -900.49586;
      c[8] = +661.45441;
      c[9] = -259.58894;
      c[10] = 42.2607181;

      CC=sqrt(197.326);
      for(Int_t i=0;i<11;i++){
	m[i]=0.231538+i*0.9;
      }
      auto spart = [&CC,&c,&m](double p){
		     p*=1000;
		     double u=0;
		     for(Int_t i=0;i<11;i++){
		       u=u+c[i]*CC/(p*p+m[i]*m[i]*CC*CC*CC*CC);
		     }
		     u=u*sqrt(4*TMath::Pi());
		     return u;
		   };

      Double_t D[11];
      D[0] = 0.022623762;
      D[1] = -0.50471056;
      D[2] = 0.56278897;
      D[3] = -16.079764;
      D[4] = 111.26803;
      D[5] = -446.6749;
      D[6] = 1098.5907;
      D[7] = -1611.4995;
      D[8] = 1374.064275;
      D[9] = -630.4371973;
      D[10] = 120.6876674;

  
      auto dpart = [&CC,&D,&m](double p){
		     p*=1000;

		     double w=0;
		     for(Int_t i=0;i<11;i++){
		       w=w+D[i]*CC/(p*p+m[i]*m[i]*CC*CC*CC*CC);
		     }

		     w=w*sqrt(4*TMath::Pi());
		     return w;

		   };

      //Now fill histogram
      Double_t binOffset = (pmax-pmin)/Nbins/2;
      TH1D hist{"Dmomentum","Dmomentum",Nbins,pmin-binOffset,pmax-binOffset};
      for(Int_t ibin=1;ibin<=Nbins;ibin++){
	Double_t p=hist.GetBinCenter(ibin);
	Double_t s2=spart(p);s2*=s2;
	Double_t d2=dpart(p);d2*=d2;
	Double_t val=p*p*(s2+d2)*(5.16137052874752581e+08);
	hist.SetBinContent(ibin,val);
      }

      return new DistTH1{hist};
    }
  }

  ///name space QuasiFree
}
