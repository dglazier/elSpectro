set(ELSPECTRO elSpectro)
set(DIRNAME core)
set(BASE $ENV{ELSPECTRO})



ROOT_GENERATE_DICTIONARY(G__${ELSPECTRO}
  Particle.h
  DecayingParticle.h
  DecayModel.h
  PhaseSpaceDecay.h
  DecayModelQ2W.h
  DecayVectors.h
  TwoBodyFlat.h
  ScatteredElectron_xy.h
  ProductionProcess.h
  ElectronScattering.h
  Distribution.h
  DistTF1.h
  DistVirtPhotFlux_xy.h
  CurrentEventInfo.h
  ParticleManager.h
  DecayManager.h
  MassPhaseSpace.h
  Manager.h
  Interface.h
  LINKDEF ElSpectroLinkDef.h)

add_library(${ELSPECTRO} SHARED
  Particle.cpp
  DecayingParticle.cpp
  DecayModel.cpp
  PhaseSpaceDecay.cpp
  DecayModelQ2W.cpp
  DecayVectors.cpp
  TwoBodyFlat.cpp
  ScatteredElectron_xy.cpp
  ProductionProcess.cpp
  ElectronScattering.cpp
  Distribution.cpp
  DistTF1.cpp
  DistVirtPhotFlux_xy.cpp
  CurrentEventInfo.cpp
  ParticleManager.cpp
  DecayManager.cpp
  MassPhaseSpace.cpp
  Manager.cpp
  G__${ELSPECTRO}.cxx)

target_link_libraries(${ELSPECTRO} ${ROOT_LIBRARIES})

install(TARGETS ${ELSPECTRO} 
  LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}")

FILE(GLOB ROOT_ELSPECTRO_PCM ${CMAKE_BINARY_DIR}/*pcm ${CMAKE_BINARY_DIR}/${DIRNAME}/*pcm)

install (FILES  ${ROOT_ELSPECTRO_PCM}    
  DESTINATION "${CMAKE_INSTALL_LIBDIR}")
