set(demos
    kww_countterms
    runkww
    )

foreach(app ${demos})
    add_executable(${app} ${app}.c)
    target_include_directories(${app} PRIVATE ${CMAKE_CURRENT_BINARY_DIR} ${kww_SOURCE_DIR}/lib)
    target_link_libraries(${app} ${kww_LIBRARY})
endforeach()
