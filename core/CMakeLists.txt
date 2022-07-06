set(ELSPECTRO elSpectro)
set(DIRNAME core)
set(BASE $ENV{ELSPECTRO})

include_directories($ENV{JPACPHOTO}/include)
include_directories($ENV{JPACPHOTO}/include/amplitudes)

ROOT_GENERATE_DICTIONARY(G__${ELSPECTRO}
  Particle.h
  DecayingParticle.h
  CollidingParticle.h
  DecayModel.h
  PhaseSpaceDecay.h
  NuclearBreakup.h
  Bremsstrahlung.h
  SDMEDecay.h
  TensorSDMEDecay.h
  VectorSDMEDecay.h
  DecayModelW.h
  DecayModelQ2W.h
  DecayModelst0.h
  DecayModelst.h
  JpacModelst.h
  GenericModelst.h
  DecayVectors.h
  DontDecay.h
  TwoBodyFlat.h
  TwoBody_stu.h
  ScatteredElectron_xy.h
  QuasiFreeNucleon.h
  BremstrPhoton.h
  ProductionProcess.h
  ElectronScattering.h
  PhotoProduction.h
  Distribution.h
  DistTF1.h
  DistConst.h
  DistUniform.h
  DistTH1.h
  DistTH2.h
  DistFlatMass.h
  SDME.h
  PhotonPolarisationVector.h
  DistVirtPhotFlux_xy.h
  ParticleManager.h
  DecayManager.h
  MassPhaseSpace.h
  Writer.h
  HepMC3Writer.h
  LundWriter.h
  GlueXWriter.h
  EICSimpleWriter.h
  FunctionsForJpac.h
  Manager.h
  Interface.h
  amplitude_blend.hpp
  LINKDEF ElSpectroLinkDef.h
  )

add_library(${ELSPECTRO} SHARED
  Particle.cpp
  DecayingParticle.cpp
  CollidingParticle.cpp
  DecayModel.cpp
  PhaseSpaceDecay.cpp
  NuclearBreakup.cpp
  Bremsstrahlung.cpp
  SDMEDecay.cpp
  TensorSDMEDecay.cpp
  VectorSDMEDecay.cpp
  DecayModelW.cpp
  DecayModelQ2W.cpp
  DecayModelst.cpp
  DecayModelst0.cpp
  JpacModelst.cpp
  GenericModelst.cpp
  DecayVectors.cpp
  DontDecay.cpp
  TwoBodyFlat.cpp
  TwoBody_stu.cpp
  ScatteredElectron_xy.cpp
  QuasiFreeNucleon.cpp
  BremstrPhoton.cpp
  ProductionProcess.cpp
  ElectronScattering.cpp
  PhotoProduction.cpp
  Distribution.cpp
  DistTF1.cpp
  DistTH1.cpp
  DistTH2.cpp
  DistFlatMass.cpp
  DistVirtPhotFlux_xy.cpp
  SDME.cpp
  PhotonPolarisationVector.cpp
  ParticleManager.cpp
  DecayManager.cpp
  MassPhaseSpace.cpp
  Writer.cpp
  HepMC3Writer.cpp
  LundWriter.cpp
  GlueXWriter.cpp
  EICSimpleWriter.cpp
  Manager.cpp
  FunctionsForJpac.cpp
  amplitude_blend.cpp
  G__${ELSPECTRO}.cxx
  )

target_link_libraries(${ELSPECTRO}   ROOT::Core ROOT::Rint ROOT::RIO ROOT::RooFit ROOT::MathMore ROOT::EG ROOT::GenVector )

install(TARGETS ${ELSPECTRO} 
  LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}")

install (FILES  ${CMAKE_CURRENT_BINARY_DIR}/libelSpectro_rdict.pcm    
  DESTINATION "${CMAKE_INSTALL_LIBDIR}")
