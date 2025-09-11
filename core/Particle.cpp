#include "Particle.h"
#include "DistFlatMass.h"
#include <TDatabasePDG.h>
#include <iostream>

namespace elSpectro{

  Particle::Particle() = default;
  Particle::~Particle() = default;
  Particle::Particle(const Particle& other) = default;
  Particle::Particle(Particle&&) = default;
  Particle& Particle::operator=(const Particle& other) = default;
  Particle& Particle::operator=(Particle&& other) = default;

  Particle::Particle(int pdg)
    : _pdg{pdg}
  {
    TDatabasePDG *pdgDB = TDatabasePDG::Instance();
    auto particle = pdgDB->GetParticle(pdg);

    if (!particle) {
      std::cerr << "Particle::Particle pdg " << pdg << " does not exist in PDG table" << std::endl;
      exit(1);
    }

    _pdgMass = particle->Mass();
    _dynamicMass = _pdgMass;
    SetXYZT(0, 0, 0, _pdgMass); // initialise at rest
  }

  void Particle::Print() const
  {
    std::cout << "Particle: PDG = " << Pdg()
              << ", minimum mass = " << MinimumMassPossible()
              << ", mass = " << _dynamicMass
              << ", PDG mass = " << _pdgMass << "\n"
              << " P4 = ("
              << P4().X() << ", "
              << P4().Y() << ", "
              << P4().Z() << ", "
              << P4().T() << ")\n";
  }

  void Particle::SetXYZT(double xx, double yy, double zz, double tt)
  {
    _vec.SetXYZT(xx, yy, zz, tt);
    _dynamicMass = _vec.M();
  }

  void Particle::SetXYZ(double xx, double yy, double zz)
  {
    auto m2 = _vec.M2();
    auto P2 = xx * xx + yy * yy + zz * zz;
    _vec.SetXYZT(xx, yy, zz, TMath::Sqrt(P2 + m2));
  }

  void Particle::SetP4(const LorentzVector& p4)
  {
    _vec = p4;
    _dynamicMass = _vec.M();
  }

  void Particle::SetP4M(double mm)
  {
    auto P2 = _vec.P2();
    _vec.SetXYZT(_vec.X(), _vec.Y(), _vec.Z(), TMath::Sqrt(P2 + mm * mm));
    _dynamicMass = mm;
  }

  void Particle::DetermineDynamicMass(double xmin, double xmax)
  {
    if(_massDist == nullptr ){
      TakePdgMass();
      return;
    }
    if(_massLocked) return;

    _dynamicMass = -1;
    _massWeight = 0;
    auto minposs = MinimumMassPossible();

    auto minRange = (xmin == -1) ? minposs : xmin;
    if(minRange > minposs) minRange = minposs;
    auto maxRange = (xmax == -1) ? _massDist->GetMaxX() : xmax;
    // Handle cases where the range is numerically zero or slightly negative due to precision issues.
    if((maxRange-minRange)<0) {
      if((maxRange-minRange)>-1E-6) {
        _dynamicMass = minRange;
        return;
      }
    }
    
    if(minRange > maxRange){// Unphysical situation where the calculated minimum is greater than the maximum allowed.
      std::cout << "Warning  Particle::DetermineDynamicMass min " << minRange
                << " greater than max " << maxRange << " for " << _pdg
                << " minposs " << minposs << " " << xmax << " "
                << _massDist->GetMaxX()
                << " equal " << (minposs == _massDist->GetMaxX()) << std::endl;
      _dynamicMass = minRange;
      return;
    }

    while(_dynamicMass < minposs){
      _dynamicMass = _massDist->SampleSingle(minRange, maxRange);
      _massWeight = _massDist->GetCurrentWeight();

      if(_dynamicMass == 0) {
        std::cout << "Error  Particle::DetermineDynamicMass zero mass" << std::endl;
        std::cout << _pdg << "  DetermineDynamicMass( " << MinimumMassPossible()
                  << " " << _dynamicMass << " " << _massWeight << " "
                  << minRange << " " << maxRange << std::endl;
        exit(1);
      }
    }
    SetP4M(_dynamicMass);
  }

  void Particle::SetParent(Particle* parent){
    //need to update the mass distribution parent pointer
    if(dynamic_cast<DistFlatMassMaster*>(MassDistribution()) ){
      auto mdist = dynamic_cast<DistFlatMassMaster*>(MassDistribution());
      mdist->SetParentPtr(dynamic_cast<DecayingParticle*>(parent));
    }
    
  }
}
