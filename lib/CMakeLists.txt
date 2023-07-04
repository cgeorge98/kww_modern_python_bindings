set(lib kww)
set(${lib}_LIBRARY ${lib} PARENT_SCOPE)

set(src_files kww.c kww_lowlevel.c)
set(inc_files kww.h kww_lowlevel.h)

add_library(${lib} ${src_files})

set_target_properties(
    ${lib} PROPERTIES
    OUTPUT_NAME ${lib}
    VERSION ${kww_VERSION}
    SOVERSION ${kww_VERSION_MAJOR})

if(USE_FLOAT128)
    target_link_libraries(${lib} quadmath)
endif()

include(LinkLibMath)
link_libm(${lib})

install(
    TARGETS ${lib}
    LIBRARY DESTINATION ${destination}/lib
    RUNTIME DESTINATION ${destination}/lib
    ARCHIVE DESTINATION ${destination}/lib
    COMPONENT Libraries)
install(
    FILES ${inc_files}
    DESTINATION ${destination}/include
    COMPONENT Headers)
