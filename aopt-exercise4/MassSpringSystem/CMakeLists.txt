
add_library (MassSpringSystem INTERFACE
#        MassSpringSystemT_impl.hh
#        MassSpringSystemT.hh
        )

add_library (AOPT::MassSpringSystem ALIAS MassSpringSystem)

target_link_libraries(MassSpringSystem INTERFACE AOPT::AOPT)

#set_target_properties(MassSpringSystem PROPERTIES
#        CXX_STANDARD 14
#        CXX_STANDARD_REQUIRED YES
#        )

target_include_directories(MassSpringSystem
        INTERFACE
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
        )

target_compile_definitions(MassSpringSystem
        INTERFACE
        "INCLUDE_TEMPLATES"
        )