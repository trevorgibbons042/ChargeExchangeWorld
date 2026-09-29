#setup, project name: g4outline
cmake_minimum_required(VERSION 3.5 FATAL_ERROR)
project(g4outline)

#finds Geant4 packages, ui_all and vis_all make GUI stuff, gdml for cad files
find_package(Geant4 REQUIRED ui_all vis_all)
include(${Geant4_USE_FILE})
include_directories(${PROJECT_SOURCE_DIR}/include)

#finds custom cc and hh, is a folder tho
file(GLOB sources ${PROJECT_SOURCE_DIR}/*.cc) #finds .cc
file(GLOB sources ${PROJECT_SOURCE_DIR}/1cc/*.cc) # finds .cc in 1cc folder
file(GLOB headers ${PROJECT_SOURCE_DIR}/2hh/*.hh) # finds .hh in 2hh folder

file(GLOB MACRO_FILES ${PROJECT_SOURCE_DIR}/3mac/*.mac)
file(COPY ${MACRO_FILES} DESTINATION ${CMAKE_BINARY_DIR}) #copies macros to build folder
file(GLOB DATA_FILES ${PROJECT_SOURCE_DIR}/4dat/*.dat)
file(COPY ${DATA_FILES} DESTINATION ${CMAKE_BINARY_DIR}/4dat) #copies data to build folder

#executes the main cc, name of exe: sim.exe, name of cc: 1sim.cc
add_executable(sim 1sim.cc ${sources} ${headers})
target_include_directories(sim PRIVATE ${PROJECT_SOURCE_DIR}/2hh)
target_link_libraries(sim PRIVATE ${Geant4_LIBRARIES})
add_custom_target(g4outline DEPENDS sim)