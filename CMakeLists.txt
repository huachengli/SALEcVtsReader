cmake_minimum_required(VERSION 3.11)

project(SALEcVtsReader C)
SET(CMAKE_C_STANDARD 99)
SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fopenmp -Wno-unused-result -O3")

LINK_LIBRARIES(z m)
SET(SRCS src/VtsTag.c
        src/VtsReader.c
        src/InputParser.c
        src/Utility.c
        src/Statistic.c)

SET(VTP_SRCS src/VtpReader.c
        src/VtpTag.c
        src/VtkWriter.c)

SET(CITCOMS_SRCS src/Utility2d.c)
add_executable(SALEcVtsReader src/main.c  ${SRCS})
add_executable(Remnant src/remnant.c ${SRCS})
add_executable(TracerView src/Tracer.c ${SRCS} src/TracerReader.c src/TracerReader.h src/AVLTree.c src/AVLTree.h)
add_executable(VtpTracer src/VtpTracer.c ${SRCS} ${VTP_SRCS} )
add_executable(CitcomsExport src/ToCitcoms2d.c ${SRCS} ${CITCOMS_SRCS})
target_compile_options(CitcomsExport PUBLIC -O0)