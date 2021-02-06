#----------------------------------------------------------------------------
# Setup the project
cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(The-MicroBlackHole)

#----------------------------------------------------------------------------
# Find Geant4 package, activating all available UI and Vis drivers by default
# You can set WITH_GEANT4_UIVIS to OFF via the command line or ccmake/cmake-gui
# to build a batch mode only executable
#
option(WITH_GEANT4_UIVIS "Build example with Geant4 UI and Vis drivers" ON)
if(WITH_GEANT4_UIVIS)
  find_package(Geant4 REQUIRED ui_all vis_all)
else()
  find_package(Geant4 REQUIRED)
endif()

#----------------------------------------------------------------------------
# Setup Geant4 include directories and compile definitions
# Setup include directory for this project
#
include(${Geant4_USE_FILE})
include_directories(${PROJECT_SOURCE_DIR}/include)


#----------------------------------------------------------------------------
# Locate sources and headers for this project
# NB: headers are included so they will show up in IDEs
#
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)

#----------------------------------------------------------------------------
# Add the executable, and link it to the Geant4 libraries
#
#set(CRY_LIB "/home/pablom/Documents/cry/cry_v1.7/lib/libCRY.a")
#set(CRY_LIB "/home/pablom/Documents/geant4/work/MTP/cry/cry_v1.7/lib/libCRY.a")
#"/home/carlos/geant4/MPWork/MTP/cry/cry_v1.7/lib/libCRY.a"
#set(CRY_LIB $ENV{G4CRYLIBDIR})
set(CRY_LIB "$ENV{G4WORKDIR}/MTP/cry/cry_v1.7/lib/libCRY.a")
find_library(CRY_LIB CRY $ENV{CRY_LIB_PATH})
include_directories("$ENV{G4WORKDIR}/MTP/cry/cry_v1.7/src")


#Leibniz
#include_directories("/home/mpablo/Documents/jsoncpp/include/")
#set(JSONCPP_LIB "/home/mpablo/Documents/jsoncpp/lib/libjsoncpp.a")
#include_directories("/home/pablom/Documents/geant4/work/MTP/cry/cry_v1.7/src/")

#server
#include_directories("/home/pablom/Documentos/work/jsoncpp/include/")
#set(JSONCPP_LIB "/home/pablom/Documentos/work/jsoncpp/lib/libjsoncpp.a")
#find_library(JSONCPP_LIB JSONCPPY $ENV{JSONCPP_LIB_PATH})
#include_directories("/home/pablom/Documentos/work/MTP/cry/cry_v1.7/src/")
include_directories("/home/pablom/jsoncpp/include/")
set(JSONCPP_LIB "/home/pablom/jsoncpp/lib/libjsoncpp.a")
find_library(JSONCPP_LIB JSONCPPY $ENV{JSONCPP_LIB_PATH})

# CADMesh
find_package(cadmesh)

# Find JSONCPP library on local system.
#find_package(JSONCPP)

#if(NOT APPLE)
#  include_directories(SYSTEM AFTER "/home/mpablo/Documents/jsoncpp/include/")
#else( )
#  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -isystem \"/home/mpablo/Documents/jsoncpp/include/\"")
#endif( )
##############ADDED FOR JSONCPP#######################

add_executable(MicroBlackHole MicroBlackHole.cc ${sources} ${headers})
#target_link_libraries(MicroBlackHole ${Geant4_LIBRARIES} ${CRY_LIB})
#target_link_libraries(MicroBlackHole ${Geant4_LIBRARIES} jsoncpp ${CRY_LIB})
#target_link_libraries(MicroBlackHole ${Geant4_LIBRARIES} ${JSONCPP_LIB} ${CRY_LIB})
#CAD
target_link_libraries(MicroBlackHole ${Geant4_LIBRARIES} ${JSONCPP_LIB} ${CRY_LIB} ${cadmesh_LIBRARIES})



#----------------------------------------------------------------------------
# Copy all scripts to the build directory, i.e. the directory in which we
# build B1. This is so that we can run the executable directly because it
# relies on these scripts being in the current working directory.
#
set(TheWorldFunction_SCRIPTS
  )

foreach(_script ${TheWorldFunction_SCRIPTS})
  configure_file(
    ${PROJECT_SOURCE_DIR}/${_script}
    ${PROJECT_BINARY_DIR}/${_script}
    COPYONLY
    )
endforeach()

#----------------------------------------------------------------------------
# For internal Geant4 use - but has no effect if you build this
# example standalone
#
add_custom_target(The-MicroBlackHole DEPENDS MicroBlackHole)

#----------------------------------------------------------------------------
# Install the executable to 'bin' directory under CMAKE_INSTALL_PREFIX
#
install(TARGETS MicroBlackHole DESTINATION bin)


