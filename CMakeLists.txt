cmake_minimum_required (VERSION 3.0)

project (gamil-dycore Fortran)

find_package (MPI)
if (NOT MPI_Fortran_FOUND)
  set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -DNO_MPI")
else ()
  include_directories (${MPI_Fortran_INCLUDE_PATH})
endif ()

if (CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
  set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fdefault-real-8 -ffpe-summary=none -ffree-line-length-none -Ofast -flto")
elseif (CMAKE_Fortran_COMPILER_ID STREQUAL "PGI")
  set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -r8 -traceback -g")
endif ()

if (NOT $ENV{NETCDF_ROOT} STREQUAL "")
  include_directories ("$ENV{NETCDF_ROOT}/include")
  link_directories ("$ENV{NETCDF_ROOT}/lib")
else ()
  message (FATAL_ERROR "CHECK")
  find_package (PkgConfig)
  if (PKG_CONFIG_FOUND)
    pkg_search_module (NETCDF REQUIRED netcdf)
    if (NETCDF_FOUND)
      include_directories (${NETCDF_INCLUDEDIR})
      link_directories (${NETCDF_LIBDIR})
    else ()
      message (FATAL_ERROR "Pkg-config could not find netcdf library!")
    endif ()
  else ()
    message (FATAL_ERROR "Unable to find pkg-config library!")
  endif ()
endif ()

add_subdirectory (${PROJECT_SOURCE_DIR}/lib/datetime)

set (source_files
  "src/string_mod.F90"
  "src/time_mod.F90"
  "src/io_mod.F90"
  "src/log_mod.F90"
  "src/params_mod.F90"
  "src/map_mod.F90"
  "src/mesh_mod.F90"
  "src/parallel_mod.F90"
  "src/types_mod.F90"
  "src/diag_mod.F90"
  "src/history_mod.F90"
  "src/restart_mod.F90"
  "src/data_mod.F90"
  "src/reduce_mod.F90"
  "src/diffusion_mod.F90"
  "src/dycore_mod.F90"
  "src/dycore_test.F90"
  "src/test_cases/barotropic/rossby_haurwitz_wave_test_mod.F90"
  "src/test_cases/barotropic/steady_geostrophic_flow_test_mod.F90"
  "src/test_cases/barotropic/mountain_zonal_flow_test_mod.F90"
  "src/test_cases/barotropic/jet_zonal_flow_test_mod.F90"
  "lib/quadpack.f90"
)

set (CMAKE_Fortran_MODULE_DIRECTORY "${CMAKE_BINARY_DIR}/.mods")

add_executable (dycore_test.exe ${source_files})
target_link_libraries (dycore_test.exe ${MPI_Fortran_LIBRARIES} fortran_datetime netcdf netcdff)
if (MPI_Fortran_FOUND)
  set_target_properties (dycore_test.exe PROPERTIES
    LINK_FLAGS "${MPI_Fortran_LINK_FLAGS}"
  )
endif ()
