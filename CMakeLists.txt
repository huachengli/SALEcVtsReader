cmake_minimum_required(VERSION 3.11)

project(SALEcVtsReader C)
SET(CMAKE_C_STANDARD 99)
SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O2 -fopenmp -Wno-unused-result")

LINK_LIBRARIES(z m)
SET(SRCS src/VtsTag.c
        src/VtsReader.c
        src/InputParser.c
        src/Utility.c
        src/Statistic.c)

add_executable(SALEcVtsReader src/main.c  ${SRCS})
add_executable(Remnant src/remnant.c ${SRCS})
