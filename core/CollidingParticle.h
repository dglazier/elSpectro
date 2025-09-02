/**
 * @file CollidingParticle.h
 * @author D. Glazier
 * @brief Defines the CollidingParticle class.
 * @date 2025-09-02
 */

#ifndef COLLIDINGPARTICLE_H
#define COLLIDINGPARTICLE_H

#include "Particle.h"
#include "DecayChannel.h" // Include the DecayChannel header
#include "DecayModel.h"   // For accessing product particles
#include "FunctionsForGenvector.h"

namespace elSpectro {

class DecayVectors; // Forward declaration

/**
 * @struct BeamProperties
 * @brief A structure to hold properties of a beam particle.
 */
struct BeamProperties {
    Double_t _dirTheta = 0.0;       ///< Beam polar angle in radians.
    Double_t _dirPhi = 0.0;         ///< Beam azimuthal angle in radians.
    Double_t _sizeHor = 0.0;        ///< Horizontal beam spot size in cm.
    Double_t _sizeVer = 0.0;        ///< Vertical beam spot size in cm.
    Double_t _divHor = 0.0;         ///< Horizontal beam divergence in radians.
    Double_t _divVer = 0.0;         ///< Vertical beam divergence in radians.
};


/**
 * @class CollidingParticle
 * @brief Represents a particle that participates in an initial collision, like a beam or target.
 */
class CollidingParticle : public Particle {
public:
    /**
     * @brief Default constructor.
     */
    CollidingParticle();

    /**
     * @brief Constructor for a simple beam particle with a defined momentum.
     * @param[in] pdg The PDG code of the interacting particle.
     * @param[in] momentum The momentum magnitude along the z-axis in GeV/c.
     */
    CollidingParticle(Int_t pdg, Double_t momentum);
    
    /**
     * @brief Constructor for a particle that is a product of another decay.
     * @param[in] pdg The PDG code of the interacting particle from the parent's decay.
     * @param[in] momentum The momentum of the parent particle in GeV/c.
     * @param[in] parentpdg The PDG code of the parent particle.
     * @param[in] model The DecayModel describing the parent's decay.
     * @param[in] decayer The DecayVectors describing the parent's decay kinematics.
     */
    CollidingParticle(Int_t pdg, Double_t momentum, Int_t parentpdg, decaymodel_ptr  mod,decayer_ptr  dec);


    /**
     * @brief Sets the decay channel for this particle by copying from another.
     * @param[in] channel A constant reference to the DecayChannel object to copy.
     */
    void SetDecayChannel(const DecayChannel& channel);

    /**
     * @brief Gets a mutable reference to the internal DecayChannel object.
     * @return A reference to the DecayChannel object.
     */
    DecayChannel& GetDecayChannel() { return _decayChannel; }

    /**
     * @brief Gets the DecayModel from the associated channel.
     * @return A pointer to the DecayModel, or nullptr if not set.
     */
  DecayModel* Model()  const {
    // Get the primary DecayModel from the owned DecayChannel.
    return (_decayChannel.N() == 0) ? nullptr:_decayChannel.CurrModel();
  }

  /**
   * @brief Gets the DecayVectors from the associated channel.
   * @return A pointer to the DecayVectors, or nullptr if not set.
   */
  DecayVectors* Decayer() {
    // Get the primary DecayVectors from the owned DecayChannel.
    return  (_decayChannel.N() == 0) ? nullptr:_decayChannel.CurrDecayer();
  }  
  /**
     * @brief Generates the final state kinematics for the event.
     * @details This is the core generation step, delegating to the DecayVectors object.
     * @return The physics weight of the generated event.
     */
    double Generate();

    /**
     * @brief Generates a complete event, including intensity and beam effects.
     * @details This method calls Generate() and then applies the model's intensity
     * and any configured beam smearing effects.
     * @return The final, weighted result for the event.
     */
    Double_t GenerateComponents();

  /**
     * @brief Post-initialization routine for the particle.
     * @details This function is called after the main object construction to perform
     * any setup that requires external information, such as linking decay products
     * to the correct vertex and initializing sub-models.
     * @param[in] info Pointer to the global ReactionInfo object.
     */
    void PostInit(ReactionInfo* info);

    // --- Setters for beam/target properties ---

    /**
     * @brief Sets the nominal four-momentum of the beam particle.
     * @param[in] p4 The nominal four-momentum.
     */
    void SetNominalP4(const LorentzVector& p4) { _nominalP4 = p4; }

    /**
     * @brief Sets the PDG code of the particle this particle interacts with.
     * @param[in] pdg The PDG code of the interaction partner.
     */
    void SetInteractingPdg(Int_t pdg) { _interactingPdg = pdg; }
    
    /**
     * @brief Sets the index of the interacting product particle from the DecayModel.
     * @param[in] idx The index in the product list. Use -1 for the nominal beam.
     */
    void SetInteractingIndex(Int_t idx) { _interactIdx = idx; }

    /**
     * @brief Sets the central direction of the beam.
     * @param[in] theta The polar angle in radians.
     * @param[in] phi The azimuthal angle in radians.
     */
    void SetAngleThetaPhi(Double_t theta, Double_t phi) {
      _beamProps._dirTheta = theta; _beamProps._dirPhi = phi;
      auto p4=_nominalP4; //copy 4-vector
      //and rotate it
      genvector::LorentzRotateY(p4,theta);
      genvector::LorentzRotateZ(p4,phi);
      //Set the rotated vector
      SetP4(p4);
      _nominalP4=p4;
    }

    /**
     * @brief Sets the beam spot size at the interaction point.
     * @param[in] hor The horizontal size (e.g., sigma) in cm.
     * @param[in] ver The vertical size (e.g., sigma) in cm.
     */
    void SetBeamSpotSize(Double_t hor, Double_t ver) { _beamProps._sizeHor = hor; _beamProps._sizeVer = ver; }

    /**
     * @brief Sets the beam divergence.
     * @param[in] hor The horizontal divergence (e.g., sigma) in radians.
     * @param[in] ver The vertical divergence (e.g., sigma) in radians.
     */
    void SetBeamDivergence(Double_t hor, Double_t ver) { _beamProps._divHor = hor; _beamProps._divVer = ver; }

    // --- Getters for beam/target properties ---
    
    const LorentzVector& GetNominal4Vector() const { return _nominalP4; }
    Int_t GetInteractingPdg() const { return _interactingPdg; }
    Int_t GetInteractingIndex() const { return _interactIdx; }
    const BeamProperties& GetBeamProperties() const { return _beamProps; }

    /**
     * @brief Gets the four-vector of the particle to be used in the interaction.
     * @details If the interacting index is -1, this returns the nominal beam P4.
     * Otherwise, it returns the P4 of the corresponding product particle from the decay model.
     * @return A constant reference to the interacting four-vector.
     */
    const LorentzVector& GetInteracting4Vector() const {
        return _interactIdx == -1 ? GetNominal4Vector() : Model()->Product(_interactIdx)->P4();
    }

  /**
   * @brief Add additional decay channel
   */
   void AddDecay(double bratio,decaymodel_ptr  mod,decayer_ptr  dec){
     _decayChannel.AddDecay(nullptr,bratio,std::move(mod),std::move(dec));
    }

  /**
   * @brief Define Decay type as production
   */
  DecayType IsDecay() const noexcept final {return DecayType::Production;}

  /**
   * @brief Choose a decay if required
   */
    void ChooseDecay() const {
      if(_decayChannel.N() == 0) return;
      
      //Randomly select a decay channel based on branching ratio
      bool good = false;
      while(good==false){
	auto idec = _decayChannel.ChooseDecay(Mass());
	//recurse daughter particles
	good = Model()->SampleMasses(Mass());
	Model()->ChooseDecay();
      }
    }
  
private:
    DecayChannel _decayChannel; ///< The single decay channel owned by this particle.

    // Beam and target properties
    LorentzVector _nominalP4;      ///< Nominal beam four-momentum.
    Int_t _interactingPdg = 0;      ///< PDG code of the particle this one interacts with.
    Int_t _interactIdx = -1;        ///< Index of the interacting product particle.
    BeamProperties _beamProps;      ///< Struct containing all beam-specific properties.

    ClassDef(CollidingParticle, 1);
};
} // namespace elSpectro

#endif

