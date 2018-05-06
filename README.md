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

# Usages

```
$ cd build
$ FC=gfortran cmake ..
$ make
$ cd ../run
$ ../build/dycore_test.exe namelist.rh_test
```

# Authors

- Bin Wang
- Li Dong
- ...
