set(ELSPECTRO elSpectro)
set(DIRNAME core)
set(BASE $ENV{ELSPECTRO})

include_directories($ENV{JPACPHOTO}/include)


ROOT_GENERATE_DICTIONARY(G__${ELSPECTRO}
  Particle.h
  DecayingParticle.h
  DecayModel.h
  PhaseSpaceDecay.h
  VectorSDMEDecay.h
  DecayModelQ2W.h
  DecayModelst.h
  DecayGammaN_Test.h
  JpacModelst.h
  DecayVectors.h
  TwoBodyFlat.h
  TwoBody_stu.h
  ScatteredElectron_xy.h
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
  CurrentEventInfo.h
  ParticleManager.h
  DecayManager.h
  MassPhaseSpace.h
  Writer.h
  HepMC3Writer.h
  LundWriter.h
  FunctionsForJpac.h
  Manager.h
  Interface.h
  LINKDEF ElSpectroLinkDef.h
  )

add_library(${ELSPECTRO} SHARED
  Particle.cpp
  DecayingParticle.cpp
  DecayModel.cpp
  PhaseSpaceDecay.cpp
  VectorSDMEDecay.cpp
  DecayModelQ2W.cpp
  DecayModelst.cpp
  JpacModelst.cpp
  DecayGammaN_Test.cpp
  DecayVectors.cpp
  TwoBodyFlat.cpp
  TwoBody_stu.cpp
  ScatteredElectron_xy.cpp
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
  CurrentEventInfo.cpp
  ParticleManager.cpp
  DecayManager.cpp
  MassPhaseSpace.cpp
  Writer.cpp
  HepMC3Writer.cpp
  LundWriter.cpp
  Manager.cpp
  FunctionsForJpac.cpp
  G__${ELSPECTRO}.cxx
  )

target_link_libraries(${ELSPECTRO} ${ROOT_LIBRARIES} )

install(TARGETS ${ELSPECTRO} 
  LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}")

FILE(GLOB ROOT_ELSPECTRO_PCM  ${CMAKE_BINARY_DIR}/${DIRNAME}/libelSpectro_rdict.pcm)

install (FILES  ${ROOT_ELSPECTRO_PCM}    
  DESTINATION "${CMAKE_INSTALL_LIBDIR}")

#include_directories(${CMAKE_CURRENT_SOURCE_DIR}/src)
#add_subdirectory (${CMAKE_CURRENT_SOURCE_DIR}/src)
