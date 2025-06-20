project(DSARRFIT Fortran)
cmake_minimum_required(VERSION 2.8)

add_executable(dsarrfit.x ${PROJECT_SOURCE_DIR}/src/dsarrfit.f)

install(TARGETS dsarrfit.x)
