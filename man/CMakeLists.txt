function(one_page pname section)
    add_custom_command(
        OUTPUT "${pname}.${section}"
        COMMAND pod2man -s ${section} -c "kww manual"
            "${CMAKE_CURRENT_SOURCE_DIR}/${pname}.pod"
            "${CMAKE_CURRENT_BINARY_DIR}/${pname}.${section}"
        DEPENDS ${pname}.pod
    )
    add_custom_command(
        OUTPUT ${pname}.html
        COMMAND pod2html --title="kww manual" --noindex
            ${CMAKE_CURRENT_SOURCE_DIR}/${pname}.pod
            > ${CMAKE_CURRENT_BINARY_DIR}/${pname}.html
        DEPENDS ${pname}.pod
    )
    install(
        FILES ${CMAKE_CURRENT_BINARY_DIR}/${pname}.${section}
        DESTINATION "${CMAKE_INSTALL_PREFIX}/man/man${section}"
    )
endfunction()

add_custom_target(
    man ALL
    DEPENDS kww.3
)

add_custom_target(
    html ALL
    DEPENDS kww.html
)

one_page(kww 3)
