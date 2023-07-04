# test whether demo programs run
add_test(NAME countterms WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/lib COMMAND kww_countterms 10 10)
add_test(NAME runkww     WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/lib COMMAND runkww 0 c a .5 .5)

# test whether numeric results agree with our reference

add_executable(kwwtest kwwtest.c)
target_include_directories(kwwtest PRIVATE ${CMAKE_SOURCE_DIR}/lib)
target_link_libraries(kwwtest ${kww_LIBRARY})
add_test(NAME kwwtest COMMAND kwwtest WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
