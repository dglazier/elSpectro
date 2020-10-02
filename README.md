# elSpectro 

Event Generator framework for incorporating Spectroscopy into electro / photoproduction reactions.

# Prerequisites

ROOT 6
HepMC3

# Installation

     git clone https://github.com/dglazier/elSpectro

     cd elSpectro

 # Set Environment

     setenv ELSPECTRO /path/to/elSpectro (or $PWD)
     setenv JPACPHOTO /path/to/jpacPhoto
     setenv HEPMC3 /path/to/HepMC3/

     mkdir build

     cd build

     cmake ../

     make install

     cmake ../

     make install

Note the double cmake ../; make install; is required for ROOT pcm files 

## Running examples

     cd examples

### Examples comparing to simple weighting of TGenPhaseSpace

1) Decay a rho meson to 2 pi

       root Load.C  ComparePhaseSpaceto2.C

2) Decay g+p -> rho(pi+,pi-) p

       root Load.C ComparePhaseSpaceto3Rho.C

3) Decay g+p -> X(rho(pi+,pi-) , phi(K+,K-) ) p

       root  Load.C ComparePhaseSpaceto5RhoPhi.C

### Examples of ElectroProduction

1) EIC e + p -> e' rho(pi+,pi-) p

        root Load.C RunRhoProton.C
        hW.Draw()
        hQ2.Draw()
        hRhoM.Draw()

### Examples of ElectroProduction of Jpac amplitudes

1) e + p -> e' Y (Jpsi (e+e-)rho(pi+,pi-)) p

       root Load.C JpacAmpVectorJpsiPiPi_hepmc3.C

### Examples of MesonEx Quasi-real PhotoProduction

Note forward tagger acceptance can be included with lines like

     production->SetLimitTarRest_eThmin(1.5*TMath::DegToRad());
     production->SetLimitTarRest_eThmax(5.5*TMath::DegToRad());
     production->SetLimitTarRest_ePmin(0.4);
     production->SetLimitTarRest_ePmax(6);


1) e + p -> e' X (pi+pi-) p

      root Load.C MesonEx_p2pi.C

2) e + p -> e' X (pi+pi+pi-) n

      root Load.C MesonEx_n3pi.C
 
3)  e + p -> e' P_c -> Jpsi(e+e-) p

      root Load.C MesonEx_JpsiPenta.C

Also just does phase space Jpsi, which does not need jpacPhoto, see code for details