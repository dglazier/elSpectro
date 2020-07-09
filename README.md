# elSpectro 

Event Generator framework for incorporating Spectroscopy into electro / photoproduction reactions.

# Prerequisites

ROOT 6

# Installation

     git clone https://github.com/dglazier/elSpectro

     mkdir build

     cd build

     cmake ../

     make install

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

        root Load.C
        hW.Draw()
        hQ2.Draw()
        hRhoM.Draw()