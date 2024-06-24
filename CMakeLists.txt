cmake_minimum_required(VERSION 3.13.0)

# Define Project Name and Version
project(autoGR VERSION 0.7.5 LANGUAGES NONE)

if (NOT CMAKE_Fortran_COMPILER AND NOT DEV)
  message(FATAL_ERROR
    "\n !!! Error. Could not find the Fortran compiler specification \
    (the CMAKE_Fortran_COMPILER variable). !!! \
    \n\n !!! A Fortran compiler must be specified explicitly before proceeding. !!! \
    \n\n !!! Also, check that you used 'cmake -C ...', not 'cmake -c ...' \
    - this could be another reason why this error message triggers. !!!")
endif()

enable_language(Fortran)

mark_as_advanced(CLEAR
  CMAKE_Fortran_COMPILER
  CMAKE_Fortran_FLAGS
)

SET(CMAKE_Fortran_MODULE_DIRECTORY "${PROJECT_BINARY_DIR}/include")
MESSAGE(${CMAKE_Fortran_MODULE_DIRECTORY})


add_library(autoGR)

set_target_properties(autoGR
  PROPERTIES
  ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}
  LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}
  RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR})

set_target_properties(autoGR PROPERTIES LINKER_LANGUAGE Fortran)

set_target_properties(autoGR PROPERTIES OUTPUT_NAME autoGR)

TARGET_INCLUDE_DIRECTORIES(autoGR PRIVATE
  ${CMAKE_Fortran_MODULE_DIRECTORY})

install(TARGETS autoGR
  EXPORT autoGRConfig
  INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})

INSTALL(DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY}/
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
  PATTERN "*.mod")

install(EXPORT autoGRConfig NAMESPACE autoGR:: DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/autoGR)
include(CMakePackageConfigHelpers)
write_basic_package_version_file(
  ${CMAKE_CURRENT_BINARY_DIR}/autoGRConfigVersion.cmake
  VERSION ${PROJECT_VERSION}
  COMPATIBILITY AnyNewerVersion)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/autoGRConfigVersion.cmake
  DESTINATION lib/cmake/autoGR)
target_include_directories(autoGR INTERFACE
  $<BUILD_INTERFACE:${CMAKE_Fortran_MODULE_DIRECTORY}>)
export(EXPORT autoGRConfig NAMESPACE autoGR::
  FILE ${PROJECT_BINARY_DIR}/autoGRConfig.cmake)

target_sources(autoGR PRIVATE
    symlib/src/num_types.f90
    symlib/src/utilities.f90
    symlib/src/numerical_utilities.f90
    symlib/src/vector_matrix_utilities.f90
    symlib/src/symmetry.f90
    symlib/src/compare_structures.f90
    symlib/src/combinatorics.f90
    symlib/src/rational_mathematics.f90
    symlib/src/group_theory.f90
    symlib/src/itertools.f90
    symlib/src/classes.f90
    symlib/src/spacegroup.f90
)

target_sources(autoGR PRIVATE
    kgridGen/src/kpointgeneration.f90
    kgridGen/src/wrap_kpoints.f90
)

target_sources(autoGR PRIVATE
    src/input_structure.f90
    src/grid_utils.f90
    src/sp_hnfs.f90
    src/niggli.f90
    src/find_kgrids.f90
    src/control_file.f90
)