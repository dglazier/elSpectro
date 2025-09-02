/**
 * @file Particle.h
 * @author D. Glazier
 * @brief Defines the base Particle class for the elSpectro framework.
 * @date 2025-09-02
 */

#pragma once

#include "LorentzVector.h"
#include "Distribution.h"
#include "SDME.h"
#include <TObject.h>
#include <TMath.h>
#include <TRandom.h>
#include <vector>
#include <memory>
#include <iostream>

namespace elSpectro{

  class DecayModel;         
  class DistFlatMassMaster; 
 
  /**
   * @enum DistType
   * @brief Type of distribution (mass or mass squared).
   */
  enum class DistType {kMass, kMassSquared};

  /**
   * @enum DecayType
   * @brief Type to indicate particle decay status.
   */
  enum class DecayType{ Stable, Attached, Detached, Production };

  /**
   * @class Particle
   * @brief Controls the behaviour and properties of particles in elSpectro.
   *
   * Particle is defined by:
   *   1) Its instantaneous LorentzVector
   *   2) Any subsequent Decays
   *
   * Acts as the base class for all particle types within the framework.
   */
  class Particle {

  public:
    Particle();
    virtual ~Particle();
    Particle(const Particle& other);
    Particle(Particle&&);
    Particle& operator=(const Particle& other);
    Particle& operator=(Particle&& other);

    /**
     * @brief Construct a Particle from a PDG code.
     * @param[in] pdg The PDG code.
     */
    Particle(int pdg);
    
    /**
     * @brief Returns whether this particle is decaying.
     * @return False by default.
     */
    virtual bool IsDecaying() const {return false;}

    /**
     * @brief Get the particle's Lorentz four-vector.
     * @return Constant reference to the LorentzVector.
     */
    LorentzVector const& P4() const {return _vec;}
    /**
     * @brief Get a pointer to the particle's Lorentz four-vector.
     * @return Pointer to LorentzVector.
     */
    LorentzVector* P4ptr() {return &_vec;}
    
    /**
     * @brief Get the PDG code for this particle.
     * @return PDG integer code.
     */
    int Pdg() const {return _pdg;}

    /**
     * @brief Set the four-vector components (x, y, z, t).
     * @param[in] xx X component.
     * @param[in] yy Y component.
     * @param[in] zz Z component.
     * @param[in] tt Time component.
     */
    void SetXYZT(double xx, double yy, double zz, double tt);

    /**
     * @brief Set the spatial three-vector; time is calculated by mass.
     * @param[in] xx X component.
     * @param[in] yy Y component.
     * @param[in] zz Z component.
     */
    void SetXYZ(double xx, double yy, double zz);

    /**
     * @brief Set the four-vector from a LorentzVector.
     * @param[in] p4 The LorentzVector to copy.
     */
    void SetP4(const LorentzVector& p4);

    /**
     * @brief Set the four-vector to the maximum possible mass.
     */
    void TakeMaximumMass();

    /**
     * @brief Set the four-vector to the minimum possible mass.
     */
    void TakeMinimumMass();

    /**
     * @brief Set the four-vector to the PDG mass.
     */
    void TakePdgMass(){
      SetP4M( PdgMass() );
    }
    /**
     * @brief Boost the particle by a given velocity vector.
     * @param[in] vboost BetaVector (velocity/c).
     */
    void Boost(const  elSpectro::BetaVector& vboost ){
      _vec=ROOT::Math::VectorUtil::boost(_vec,vboost);
    }

    /**
     * @brief Get mass squared (M^2) of the particle.
     * @return Mass squared in GeV^2.
     */
    double M2() const {
      if(Pdg()==22) return 0.;
      return _dynamicMass*_dynamicMass;
    }
    /**
     * @brief Get mass of the particle.
     * @return Mass in GeV.
     */
    double Mass() const {
      if(Pdg()==22) return 0.;
      return _dynamicMass;
    }
    
    /**
     * @brief Set the mass distribution for this particle.
     * @param[in] dist Shared pointer to Distribution.
     */
    void SetMassDist(std::shared_ptr<Distribution> dist){
      _massDist = dist;
    }
   
    /**
     * @brief Get the mass distribution for this particle, if any.
     * @return Pointer to Distribution, or nullptr.
     */
    Distribution* MassDistribution() const {return _massDist.get();}
    
    /**
     * @brief Set the parent particle pointer.
     * @param[in] parent Pointer to parent Particle.
     */
    void SetParent(Particle* parent);
 
    /**
     * @brief Set the PDG mass for this particle and update four-vector mass.
     * @param[in] val The new PDG mass.
     */
    void SetPdgMass(double val){ _pdgMass=val; SetP4M(val); }
    
    /**
     * @brief Get the PDG mass for this particle.
     * @return PDG mass value.
     */
    double PdgMass() const noexcept { return _pdgMass; }

    /**
     * @brief Get the minimum physically possible mass for this particle.
     * @return Minimum mass value.
     */
    virtual double MinimumMassPossible()const noexcept {
      return  PdgMass();
    }
    /**
     * @brief Get the maximum physically possible mass for this particle.
     * @return Maximum mass value.
     */
    virtual double MaximumMassPossible()const noexcept {
      return  PdgMass();
    }
    /**
     * @brief Get the minimum mass allowed by the decay channel.
     * @return Minimum channel mass.
     */
    virtual double MinimumMassForChannel() const noexcept {
      return MinimumMassPossible();
    }
 
    /**
     * @brief Get the mass weight for this particle (used in mass distributions).
     * @return Mass weight (probability).
     */
    double MassWeight() const noexcept { return _massWeight; }

    /**
     * @brief Print the particle state to the console.
     */
    virtual void Print() const;

    /**
     * @brief Set the vertex ID associated with this particle.
     * @param[in] vertexID The vertex ID.
     */
    void SetVertexID(int vertexID){
      _vertexID=vertexID;
    }
    /**
     * @brief Set the vertex position for this particle.
     * @param[in] v The LorentzVector for the vertex.
     */
    void SetVertexPosition(const LorentzVector& v){
      _vertex=v;
    }
    /**
     * @brief Get the vertex position.
     * @return Constant reference to LorentzVector.
     */
    const LorentzVector& VertexPosition()const noexcept{return _vertex;}
    /**
     * @brief Get the vertex ID.
     * @return Vertex ID.
     */
    int VertexID()const noexcept{return _vertexID;}

    /**
     * @brief Returns decay type (Stable by default).
     * @return DecayType status.
     */
    virtual DecayType IsDecay() const noexcept {return DecayType::Stable;}
  
    /**
     * @brief Initialise the SDME for this particle.
     * @param[in] J Spin.
     * @param[in] alphaMax Maximum allowed value for alpha.
     * @return Pointer to initialised SDME.
     */
    SDME* InitSDME(uint J,uint alphaMax){
      _sdme=SDME(J,alphaMax);
      return &_sdme;
    }
    /**
     * @brief Get the SDME for this particle.
     * @return Const pointer to SDME.
     */
    const SDME* GetSDME() const noexcept{ return &_sdme; }

    /**
     * @brief Lock the mass, preventing it from being changed by mass distributions.
     */
    void LockMass(){_massLocked=true;}
    /**
     * @brief Unlock the mass, allowing mass distributions to update it.
     */
    void UnlockMass(){_massLocked=false;}
    
    /**
     * @brief Set the four-vector mass component, adjusting energy accordingly.
     * @param[in] mm The new mass.
     */
    void SetP4M(double mm);

  private:

    friend DecayModel;        
    friend DistFlatMassMaster;

    /**
     * @brief If mass comes from a distribution, sample it.
     * @param[in] xmin Minimum mass range (default -1 for auto).
     * @param[in] xmax Maximum mass range (default -1 for auto).
     * @note This function may update _dynamicMass and _massWeight.
     */
    void DetermineDynamicMass(double xmin=-1, double xmax=-1);

    LorentzVector _vec;               ///< Particle's Lorentz four-vector.
    SDME _sdme;                       ///< Spin Density Matrix Elements for the particle.
    double _pdgMass {0};              ///< PDG mass value for the particle (GeV).
    double _dynamicMass {0};          ///< Mass that may be updated by distributions (GeV).
    double _massWeight {1};           ///< Weight of the sampled mass (for distributions).
    int _pdg {0};                     ///< PDG integer code.
    int _vertexID {0};                ///< ID of the production vertex.
    LorentzVector _vertex;            ///< Position of the production vertex.
    std::shared_ptr<Distribution> _massDist {nullptr}; ///< Mass distribution, if applicable.
    bool _massLocked {false};         ///< Is the mass locked (cannot be changed by distributions)?

    ClassDef(elSpectro::Particle,1);  ///< ROOT class definition macro.
    
  }; // class Particle

  /**
   * @brief Unique pointer type to Particle.
   */
  using particle_uptr = std::unique_ptr<Particle>;

} // namespace elSpectro
