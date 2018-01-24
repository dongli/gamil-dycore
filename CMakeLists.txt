cmake_minimum_required (VERSION 2.8)

project (gamil-dycore Fortran)

find_package (MPI)
if (NOT ${MPI_Fortran_FOUND})
  set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -DNO_MPI")
else ()
  set (CMAKE_Fortran_COMPILER ${MPI_Fortran_COMPILER})
endif ()

set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fdefault-real-8 -ffpe-summary=none -ffree-line-length-none -Ofast -flto")

find_package (PkgConfig)
if (${PKG_CONFIG_FOUND})
  pkg_search_module (NETCDF REQUIRED netcdf)
  if (${NETCDF_FOUND})
    include_directories (${NETCDF_INCLUDE_DIRS})
    link_directories (${NETCDF_LIBRARY_DIRS})
  else ()
    message (FATAL_ERROR "Pkg-config could not find netcdf library!")
  endif ()
else ()
  message (FATAL_ERROR "Unable to find pkg-config library!")
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
  "src/history_mod.F90"
  "src/restart_mod.F90"
  "src/dycore_mod.F90"
  "src/dycore_main.F90"
  "src/test_cases/barotropic/rossby_haurwitz_test_mod.F90"
)

set (CMAKE_Fortran_MODULE_DIRECTORY "${CMAKE_BINARY_DIR}/.mods")

add_executable (gamil_dycore.exe ${source_files})
target_link_libraries (gamil_dycore.exe fortran_datetime netcdf netcdff)
