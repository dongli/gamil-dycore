# Introduction

This is a redesign and reimplementation of GAMIL dynamical core (shortly as
dycore) developed in IAP/LASG. New dycore should be devised with the following
considerations ordered by priority from high to low:

- Switch between hydrostatic and nonhydrostatic (e.g. hydrostatic-pressure coordinate)
- Better vertical coordinate (e.g. SLEVE)
- Nonhydrostatic (e.g. filter vertical acoustic wave or implicit solver)
- Deep atmosphere
- Good extendability

**Currently, only a barotropic dycore is implemented, and we are working on a reduced tendency scheme to improve numerical stability around Poles.**

# Prerequistes

You should install the following software with proper version yourself:

- CMake
- NetCDF
- NCL

# Usages

Obtain the codes by using `git` command:

```
$ git clone https://github.com/dongli/gamil-dycore --recursive
```

Build the model:

```
$ mkdir build
$ cd build
$ FC=gfortran cmake ..
$ make
$ cd ..
```

There should be `dycore_test.exe` file generated.

Run the Rossby-Haurwitz wave test:

```
$ cd run
$ ../build/dycore_test.exe namelist.rh_test
$ ncl -Q ../src/test_cases/barotropic/plot_rossby_haurwitz_wave_test.ncl file_prefix=\"rh_test.360x181.dt240\"
```

Run the mountain zonal flow test:

```
$ cd run
$ ../build/dycore_test.exe namelist.mz_test
$ ncl -Q ../src/test_cases/barotropic/plot_mountain_zonal_flow_test.ncl file_prefix=\"mz_test.360x181.dt240\"
```

# Authors

- Bin Wang
- Li Dong
- ...
