cmake_minimum_required(VERSION 3.11)
project(SALEcVtsReader C)

SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O0 -fopenmp -Wunused-result")
SET(CMAKE_C_STANDARD 99)
LINK_LIBRARIES(z m)

add_executable(SALEcVtsReader main.c VtsTag.c VtsReader.c InputParser.c InputParser.h Utility.c Utility.h)
