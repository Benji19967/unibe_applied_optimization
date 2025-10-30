get_filename_component(PROJECT_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)

add_executable(${PROJECT_NAME}
        main.cc
        )
target_link_libraries(${PROJECT_NAME} AOPT::AOPT)

set_target_properties(${PROJECT_NAME} PROPERTIES
        CXX_STANDARD 14
        CXX_STANDARD_REQUIRED YES
        )



add_executable(${PROJECT_NAME}-test
        unit_tests.cc
        )

target_link_libraries(${PROJECT_NAME}-test
        AOPT::AOPT
        gtest gtest_main
        )


gtest_add_tests(TARGET ${PROJECT_NAME}-test
        EXTRA_ARGS "--gtest_color=yes"
        )