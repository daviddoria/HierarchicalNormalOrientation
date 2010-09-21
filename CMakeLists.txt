cmake_minimum_required(VERSION 2.6)

PROJECT(HierarchicalNormalOrientation)

FIND_PACKAGE(VTK REQUIRED)
INCLUDE(${VTK_USE_FILE})

ADD_EXECUTABLE(HierarchicalNormalOrientation Example.cxx vtkHierarchicalNormalOrientation vtkSuperPoints.cxx)
TARGET_LINK_LIBRARIES(HierarchicalNormalOrientation vtkHybrid vtkInfovis)
