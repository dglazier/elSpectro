#ifdef __CINT__

#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;
#pragma link C++ nestedclasses;

#pragma link C++ namespace elSpectro;


//#pragma link C++ class elSpectro::LorentzVector+;
#pragma link C++ class elSpectro::Particle+;

#pragma link C++ class elSpectro::DecayingParticle+;

#pragma link C++ class elSpectro::DecayModel+;
#pragma link C++ class elSpectro::DecayModelQ2W+;
#pragma link C++ class elSpectro::DecayModelst+;
#pragma link C++ class elSpectro::JpacModelst+;
#pragma link C++ class elSpectro::GenericModelst+;
#pragma link C++ class elSpectro::DecayGammaN_Test+;
#pragma link C++ class elSpectro::PhaseSpaceDecay+;
#pragma link C++ class elSpectro::VectorSDMEDecay+;

#pragma link C++ class elSpectro::DecayVectors+;
#pragma link C++ class elSpectro::TwoBodyFlat+;
#pragma link C++ class elSpectro::TwoBody_stu+;
#pragma link C++ class elSpectro::ScatteredElectron_xy+;
#pragma link C++ class elSpectro::QuasiFreeNucleon+;

#pragma link C++ class elSpectro::ProductionProcess+;
#pragma link C++ class elSpectro::ElectronScattering+;

#pragma link C++ class elSpectro::Distribution+;
#pragma link C++ class elSpectro::DistConst+;
#pragma link C++ class elSpectro::DistUniform+;
#pragma link C++ class elSpectro::DistTF1+;
#pragma link C++ class elSpectro::DistTH1+;
#pragma link C++ class elSpectro::DistTH2+;
#pragma link C++ class elSpectro::DistFlatMass+;
#pragma link C++ class elSpectro::DistFlatMassMaster+;
#pragma link C++ class elSpectro::DistVirtPhotFlux_xy+;

#pragma link C++ class elSpectro::Writer+;
#pragma link C++ class elSpectro::LundWriter+;
#pragma link C++ class elSpectro::HepMC3Writer+;


#pragma link C++ class elSpectro::ParticleManager+;
#pragma link C++ class elSpectro::DecayManager+;
#pragma link C++ class elSpectro::MassPhaseSpace+;
#pragma link C++ class elSpectro::Manager+;

#pragma link C++ defined_in "Interface.h";
#pragma link C++ defined_in "FunctionsForElectronScattering.h";
#pragma link C++ defined_in "FunctionsForKinematics.h";
#pragma link C++ defined_in "FunctionsForGenvector.h";


#endif
