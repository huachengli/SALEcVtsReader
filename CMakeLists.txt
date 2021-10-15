cmake_minimum_required(VERSION 3.11)
project(SALEcVtsReader C)

set(CMAKE_C_STANDARD 99)
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O0")
set(CMAKE_C_COMPILER gcc)
LINK_LIBRARIES(m z)

add_executable(SALEcVtsReader main.c VtsReader.c VtsReader.h)
