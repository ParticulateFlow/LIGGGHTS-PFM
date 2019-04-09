# Usage:
#
#
# Contribution :
#   Author:
#     Daniel Queteschiner, daniel.queteschiner@jku.at
#
#   Affiliation:
#     Department of Particulate Flow Modelling
#     Institute of Fluid Mechanics and Heat Transfer
#     Faculty of Engineering & Natural Sciences
#     Johannes Kepler University Linz

CMAKE_MINIMUM_REQUIRED(VERSION 2.6)

#--------------------------------------------------
# Find and Use ParaView
#--------------------------------------------------
IF (ParaView_SOURCE_DIR)
  INCLUDE_DIRECTORIES(
    ${VTK_INCLUDE_DIRS}
  )
ELSE (ParaView_SOURCE_DIR)
  FIND_PACKAGE(ParaView REQUIRED)
  INCLUDE(${PARAVIEW_USE_FILE})
ENDIF (ParaView_SOURCE_DIR)

IF(PARAVIEW_BUILD_QT_GUI)
  # Extend the auto-generated panel
  QT4_WRAP_CPP(MOC_SRCS pqSuperquadricGlyphPanel.h)
  ADD_PARAVIEW_OBJECT_PANEL(IFACES IFACE_SRCS
    CLASS_NAME pqSuperquadricGlyphPanel
    XML_NAME SuperquadricGlyphFilter # name of SourceProxy in server manager XML
    XML_GROUP filters
    )

  ADD_PARAVIEW_PLUGIN(SuperquadricGlyphFilter "1.0.0"
    REQUIRED_ON_SERVER
    REQUIRED_ON_CLIENT
    SERVER_MANAGER_XML vtkSuperquadricGlyphFilter_SM.xml
    SERVER_MANAGER_SOURCES
        vtkSuperquadricGlyphFilter.cxx
    GUI_INTERFACES ${IFACES}
    GUI_SOURCES ${MOC_SRCS} ${UI_SRCS} ${IFACE_SRCS}
        pqSuperquadricGlyphPanel.cxx
    )
ELSE(PARAVIEW_BUILD_QT_GUI)
  ADD_PARAVIEW_PLUGIN(SuperquadricGlyphFilter "1.0.0"
    SERVER_MANAGER_SOURCES
        vtkSuperquadricGlyphFilter.cxx
    SERVER_MANAGER_XML vtkSuperquadricGlyphFilter_SM.xml
    )
ENDIF(PARAVIEW_BUILD_QT_GUI)
