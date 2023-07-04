cmake_minimum_required(VERSION 3.6 FATAL_ERROR)

set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)
include(PreventInSourceBuilds)

project(kww VERSION 3.7 LANGUAGES C)

## Options.

if(NOT DEFINED BUILD_SHARED_LIBS)
    option(BUILD_SHARED_LIBS "Build as shared library" ON)
endif()
if(NOT DEFINED LIB_MAN)
    option(LIB_MAN "Build library man pages" ON)
endif()
option(PEDANTIC "Compile with pedantic warnings" ON)
option(WERROR "Treat warnings as errors" OFF)
option(USE_FLOAT128 "Use float128. Required if long double is shorter than 80 bits" OFF)
option(PORTABLE "Under gcc, build a portable binary without host-specific optimization" OFF)

## Compiler settings.

set(CMAKE_C_STANDARD 11)
if(MSVC)
    message(FATAL_ERROR "kww must not be compiled with MSVC because it needs 80 bits long doubles")
endif()
if(WIN32)
    set(WINDOWS_EXPORT_ALL_SYMBOLS ON)
    message(STATUS "compiling for Windows")
endif()
add_compile_options(-fno-omit-frame-pointer)
if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU" OR CMAKE_C_COMPILER_ID STREQUAL "GNU" )
    if(NOT PORTABLE)
        add_compile_options(-march=native)
    endif()
    if(USE_FLOAT128)
        add_compile_options(-DUSE_FLOAT128)
    endif()
else()
    if(USE_FLOAT128)
        message(FATAL "float128 is only available under gcc")
    endif()
endif()
if(PEDANTIC)
    add_compile_options(-pedantic -Wall)
endif()
if(WERROR)
    add_compile_options(-Werror)
endif()
if (CMAKE_BUILD_TYPE MATCHES Debug)
    add_compile_options(-g)
else()
    add_compile_options(-O3)
endif()

## Install settings.

set(destination ${CMAKE_INSTALL_PREFIX})
configure_file("kww.pc.in" "kww.pc" @ONLY)
install(FILES "${CMAKE_CURRENT_BINARY_DIR}/kww.pc"
    DESTINATION "${destination}/lib/pkgconfig/")

include(CTest)

add_subdirectory(lib)
add_subdirectory(demo)
add_subdirectory(test)
if (LIB_MAN)
    add_subdirectory(man)
endif()
