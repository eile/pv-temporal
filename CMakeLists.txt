project(pv_temporal)
cmake_minimum_required(VERSION 2.8)

# --------------------------------------------------------------
# CMake setup, modules, policy, paths
# --------------------------------------------------------------
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR})

if(COMMAND cmake_policy)
  cmake_policy(SET CMP0003 OLD)
  cmake_policy(SET CMP0017 NEW)
  cmake_policy(SET CMP0023 OLD)
endif(COMMAND cmake_policy)
      
# CMake 2.6+ specific
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")

set (${PROJECT_NAME}_PACKAGE_VERSION 0)
set (${PROJECT_NAME}_PACKAGE_VERSION_MAJOR 1)
set (${PROJECT_NAME}_PACKAGE_VERSION_MINOR 0)

# --------------------------------------------------------------
# Plugin path
# --------------------------------------------------------------
get_filename_component(PLUGIN_PATH ${CMAKE_CURRENT_SOURCE_DIR} ABSOLUTE)

#--------------------------------------------------
# Location of cmake files
#--------------------------------------------------
SET(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

# --------------------------------------------------------------
# Find ParaView
# --------------------------------------------------------------
FIND_PACKAGE(ParaView REQUIRED)
INCLUDE(${PARAVIEW_USE_FILE})
# recent VTK/ParaView build uses special module files
include(vtkModuleAPI)
include(vtkModuleMacros)
include(ParaViewMacros)

# --------------------------------------------------
# Plugin config
# --------------------------------------------------
set(XML_PATH   "${CMAKE_CURRENT_SOURCE_DIR}")
set(SRC_PATH   "${CMAKE_CURRENT_SOURCE_DIR}")

#--------------------------------------------------
# Includes
#--------------------------------------------------
include_directories(
  ${PLUGIN_PATH}
)

#--------------------------------------------------
# Define plugin sources
#--------------------------------------------------

SET(PLUGIN_NAME pv_temporal)
ADD_PARAVIEW_PLUGIN(
  ${PLUGIN_NAME}
  "1.0" 

  SERVER_MANAGER_XML     
    ${XML_PATH}/vtkTemporalDifferenceFilter.xml
    ${XML_PATH}/vtkExponentialDecayFilter.xml

  SERVER_MANAGER_SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/vtkTemporalDifferenceFilter.cxx
    ${CMAKE_CURRENT_SOURCE_DIR}/vtkExponentialDecayFilter.cxx

  SERVER_SOURCES

  GUI_INTERFACES 
    ${BBP_IFACE}

  GUI_SOURCES
)

#--------------------------------------------------
# setup plugin linking
#--------------------------------------------------

TARGET_LINK_LIBRARIES(${PLUGIN_NAME} 
  PUBLIC 
# todo
)

set(pv_temporal_INCLUDES_BUILD_TIME ${CMAKE_CURRENT_SOURCE_DIR})

include(CMakeInstallation)

#--------------------------------------------------
# Install
#--------------------------------------------------
SET(INSTALL_PATH 
  "${CMAKE_INSTALL_PREFIX}/lib/paraview-${PARAVIEW_VERSION_MAJOR}.${PARAVIEW_VERSION_MINOR}"
)

INSTALL(
  TARGETS 
    ${PLUGIN_NAME}
  DESTINATION ${INSTALL_PATH}
)
