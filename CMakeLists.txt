cmake_minimum_required(VERSION 3.16)
project(MagFieldCalculator C)

set(CMAKE_C_STANDARD 99)

add_executable(MagFieldCalculator main.c geomagcalc.c geomagcalc.h)