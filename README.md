# elSpectro 

Event Generator framework for incorporating Spectroscopy into electro / photoproduction reactions.

# Prerequisites
 
ROOT 6

with MathMore RooFit GenVector EG

We require jpacPhot but this is included as a submodule. You may link to your own version of you prefer.

# Installation

     git clone --recurse-submodules https://github.com/dglazier/elSpectro

     cd elSpectro

 # Set Environment

     setenv ELSPECTRO /path/to/elSpectro (or $PWD)
     setenv PATH ${PATH}:${ELSPECTRO}/bin

Note, if you want to use jpacPhoto with jpacBox you need to add the C++ boost library location to your path

      setenv PATH ${PATH}:/where/is/boost

# Build with cmake
 
     mkdir build; cd build; cmake ../

     cmake --build . --target install


## Running examples

     cd examples

Note the --i option means retain interactive root session, if not included elsepctro will exit once script is complete.

### Examples comparing to simple weighting of TGenPhaseSpace

1) Decay a rho meson to 2 pi

      elspectro --i ComparePhaseSpaceto2.C

2) Decay g+p -> rho(pi+,pi-) p

      elspectro --i  ComparePhaseSpaceto3Rho.C

3) Decay g+p -> X(rho(pi+,pi-) , phi(K+,K-) ) p

       elspectro --i ComparePhaseSpaceto5RhoPhi.C

4)  Decay g+p -> X(rho(pi+,pi-) , 4pi ) p

       elspectro --i ComparePhaseSpaceto4PiRho.C


### Examples of ElectroProduction of Jpac amplitudes

1) e + p -> e' X (Jpsi (e+e-)rho(pi+,pi-)) p

      elspectro 'EIC_JPAC_X3872.C("high",5,41,1E33,10)'

Which will run with 5GeV e- energy, 41 GeV proton, Luminosity=10^33 for 10 days

Or with diagnostic histgrams

      elspectro 'EIC_JPAC_X3872_Hists.C("high",5,41,1E33,10)'

The first argument can be "high" or "low" giving different parameterisations.

To set luminosity and days change last 2 arguments, e.g. for luminosoty 10^33 and 25 days, e- energy 100 and p energy 100 with high energy paramterisation :

      elspectro 'EIC_JPAC_X3872.C("high",100,100,1E33,25)'

To just run a fixed number of events leave last argument 0 and nLumi=number of events

     elspectro 'EIC_JPAC_X3872.C("high",100,100,1E4)'

2) e + p -> e' Z(3900) (Jpsi (e+e-) pi+) n

To run with luminosity 10^33 for 25 days

      elspectro 'EIC_JPAC_nZc_Hists.C("low",5,41,1E33,25)'

or to just run 1000 events

      elspectro 'EIC_JPAC_nZc_Hists.C("low",5,41,1000)'

### Examples of MesonEx Quasi-real PhotoProduction

Note forward tagger acceptance can be included with lines like

     production->SetLimitTarRest_eThmin(1.5*TMath::DegToRad());
     production->SetLimitTarRest_eThmax(5.5*TMath::DegToRad());
     production->SetLimitTarRest_ePmin(0.4);
     production->SetLimitTarRest_ePmax(6);


1) e + p -> e' X (pi+pi-) p

     elspectro MesonEx_p2pi.C

2) e + p -> e' X (pi+pi+pi-) n

     elspectro MesonEx_n3pi.C
 
3)  e + p -> e' P_c -> Jpsi(e+e-) p

     elspectro MesonEx_JpsiPenta.C

Also just does phase space Jpsi, which does not need jpacPhoto, see code for details