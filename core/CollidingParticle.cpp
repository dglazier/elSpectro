#include "CollidingParticle.h"
#include "DecayModel.h"
#include "DecayVectors.h"
#include "TRandom3.h"
#include "TSystem.h"
#include "TMath.h"
#include <iostream>

namespace elSpectro {

  CollidingParticle::CollidingParticle() : Particle() {
    // Default constructor
  }

  CollidingParticle::CollidingParticle(Int_t pdg, Double_t momentum) : Particle(pdg) {
    // Constructor for a simple beam particle with a defined momentum.
    _interactingPdg = pdg;
    auto mass = PdgMass(); // Assumes Particle base class can look up mass from PDG
    LorentzVector lv(0, 0, momentum, TMath::Sqrt(momentum * momentum + mass * mass));
    SetP4(lv);
    _nominalP4 = P4();
    _interactIdx = -1;
  }

  CollidingParticle::CollidingParticle(Int_t pdg, Double_t momentum, Int_t parentpdg, decaymodel_ptr  model, decayer_ptr  decayer)
    : Particle(parentpdg) {
    
    // Set the model and decayer on the internal channel
    _decayChannel.AddDecay(nullptr,1.,std::move(model),std::move(decayer));

    // Set interacting particle pdg
    _interactingPdg = pdg;

    // Set the parent LorentzVector
    auto mass = PdgMass();
    LorentzVector lv(0, 0, momentum, TMath::Sqrt(momentum * momentum + mass * mass));
    SetP4(lv);
    _nominalP4 = P4();

    // Find the relevant particle pointer in the model
    // This is the particle which is used in the production process
    UInt_t position = 0;
    for (auto& p : model->Products()) {
      if (p->Pdg() == pdg) {
	if (_interactIdx != -1) {
	  std::cerr << "CollidingParticle::CollidingParticle, multiple particles with pdg = " << pdg << std::endl;
	  exit(0);
	} else {
	  _interactIdx = position;
	}
      }
      position++;
    }

    // We need a "nominal" 4-momentum for our interacting particle for integrations etc.
    // To do this we boost it from rest into the lab frame of the parent.
    if (_interactIdx != -1) {
      // Get the interacting particle (at rest in the parent frame)
      _nominalP4 = GetInteracting4Vector();
      // Boost it to the lab frame
      decayer->BoostToParentWithRandPhi(P4(), _nominalP4);
    }
  }


  void CollidingParticle::SetDecayChannel(const DecayChannel& channel) {
    // Use the copy-assignment operator to copy the channel's state
    _decayChannel = channel;
  }


  double CollidingParticle::Generate() {
    // Use the decayer to generate the event kinematics for the products defined by the model.
    if(_decayChannel.N() == 0) return 1.0;
    auto model = Model();
    return Decayer()->Generate(P4(), model->Products());
  }
  
  void CollidingParticle::PostInit(ReactionInfo* info) {
    // This function handles setup after the particle is fully constructed.

    // First, delegate the PostInit call to the DecayChannel.
   // This allows the channel to initialize its internal components (model and decayer)
   // and potentially make choices about the decay.
    if(_decayChannel.N() == 0) return;
    _decayChannel.PostInit(info);
   
   
   auto model = Model();
   if (model) {
      // Get the list of product particles from the decay model.
      auto& products = model->Products();

      // Assign the same vertex ID from this particle (the parent)
      // to all of its daughter particles.
      for (auto* prod : products) {
	prod->SetVertexID(this->VertexID());
      }

      // Call the PostInit function on the sub-model to continue initialization.
      model->PostInit(info);
    }
  }
    
  Double_t CollidingParticle::GenerateComponents() {
    auto model = Model();
    if (model == nullptr) {
      return 1.0;
    }

    Double_t weight = 0.0;
    // Loop to ensure a physical event is generated (weight is not zero)
    while (weight == 0) {
      weight = this->Generate();
      if (weight > 0) { // Only multiply if the generation was successful
	weight *= model->Intensity();
      }
    }

    // --- Apply beam divergence smearing ---
    /*
      if (_beamProps._divHor > 0 || _beamProps._divVer > 0) {
      LorentzVector p4 = P4();

      // Use TRandom3 for Gaussian smearing. gRandom is a ROOT global pointer.
      if (!gRandom) gRandom = new TRandom3(0);

      // Generate random smearing angles based on the beam divergence
      double dTheta = gRandom->Gaus(0, _beamProps._divVer); // Vertical divergence affects theta
      double dPhi = gRandom->Gaus(0, _beamProps._divHor);   // Horizontal divergence affects phi
        
      // Apply the smearing by rotating the momentum vector
      p4.RotateY(dTheta);
      p4.RotateZ(dPhi);

      // Set the new, smeared four-momentum for the particle
      SetP4(p4);
      }
    */
    return weight;
  }

} // namespace elSpectro
