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
  SDMEDecay.h
  TensorSDMEDecay.h
  VectorSDMEDecay.h
  DecayModelQ2W.h
  DecayModelst.h
  DecayGammaN_Test.h
  JpacModelst.h
  GenericModelst.h
  DecayVectors.h
  TwoBodyFlat.h
  TwoBody_stu.h
  ScatteredElectron_xy.h
  QuasiFreeNucleon.h
  ProductionProcess.h
  ElectronScattering.h
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
#  CurrentEventInfo.h
  ParticleManager.h
  DecayManager.h
  MassPhaseSpace.h
  Writer.h
  HepMC3Writer.h
  LundWriter.h
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
  SDMEDecay.cpp
  TensorSDMEDecay.cpp
  VectorSDMEDecay.cpp
  DecayModelQ2W.cpp
  DecayModelst.cpp
  JpacModelst.cpp
  GenericModelst.cpp
  DecayGammaN_Test.cpp
  DecayVectors.cpp
  TwoBodyFlat.cpp
  TwoBody_stu.cpp
  ScatteredElectron_xy.cpp
  QuasiFreeNucleon.cpp
  ProductionProcess.cpp
  ElectronScattering.cpp
  Distribution.cpp
  DistTF1.cpp
  DistTH1.cpp
  DistTH2.cpp
  DistFlatMass.cpp
  DistVirtPhotFlux_xy.cpp
  SDME.cpp
  PhotonPolarisationVector.cpp
#  CurrentEventInfo.cpp
  ParticleManager.cpp
  DecayManager.cpp
  MassPhaseSpace.cpp
  Writer.cpp
  HepMC3Writer.cpp
  LundWriter.cpp
  EICSimpleWriter.cpp
  Manager.cpp
  FunctionsForJpac.cpp
  amplitude_blend.cpp
  G__${ELSPECTRO}.cxx
  )

target_link_libraries(${ELSPECTRO} ${ROOT_LIBRARIES} )

install(TARGETS ${ELSPECTRO} 
  LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}")

install (FILES  ${CMAKE_CURRENT_BINARY_DIR}/libelSpectro_rdict.pcm    
  DESTINATION "${CMAKE_INSTALL_LIBDIR}")
